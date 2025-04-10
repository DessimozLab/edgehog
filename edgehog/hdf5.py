import numpy
import pandas as pd
import sys
import networkx as nx

try:
    import tables
    from pyoma.browser.tablefmt import AncestralSyntenyRels, ExtantSyntenyRels
    from pyoma.browser import db
    from pyoma.browser.models import Genome
    from pyoma import version as pyoma_version
except ImportError as e:
    print(f"pyoma and pytables libraries are required to work on HDF5 files. "
          f"Please install edgehog with the `oma` extras added, i.e. `pip install edgehog[oma]`.")
    print(e)
    sys.exit(2)


class HDF5Writer:
    def __init__(self, fname, oma_db_fn, date_edges=False, orient_edges=False):
        self.fname = fname
        self.oma_db = oma_db_fn
        self.date_edges = date_edges
        self.orient_edges = orient_edges

    def __enter__(self):
        self.oma_h5 = tables.open_file(self.oma_db, 'r')
        self.h5 = tables.open_file(self.fname, 'w', filters=tables.Filters(complevel=6, complib="blosc"))
        self.tax2taxid = self._load_taxname_2_taxid()
        self.genome_lookup = self._load_genome_tab()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.h5.close()
        self.oma_h5.close()

    def _load_taxname_2_taxid(self):
        tab = self.oma_h5.get_node("/Taxonomy")
        return {row['Name'].decode(): int(row['NCBITaxonId']) for row in tab.read()}

    def _load_genome_tab(self):
        tab = self.oma_h5.get_node("/Genome")
        return {int(row['NCBITaxonId']): row for row in tab.read()}

    def _load_hog_at_level(self, taxid):
        try:
            tab = self.oma_h5.get_node('/AncestralGenomes/tax{}/Hogs'.format(taxid))
        except tables.NoSuchNodeError:
            tab = self.oma_h5.get_node('/Hogs_per_Level/tax{}'.format(taxid))
        return {row['ID'].decode(): row_nr for row_nr, row in enumerate(tab.read())}

    def add_graph_at_level(self, taxid, tree_node):
        hogid_lookup = self._load_hog_at_level(taxid)
        dtype = tables.dtype_from_descr(AncestralSyntenyRels)
        if self.orient_edges:
            try:
                orient_enum = AncestralSyntenyRels.columns['Orientation'].enum
            except KeyError:
                print(f"[WARNING] The installed pyoma library ({pyoma_version()}) does not allow to store the Orientation results. Please update pyoma if needed")
                self.orient_edges = False
        dfs = []
        for evidence, graph in zip(
                ("linearized", "parsimonious", "any"),
                (tree_node.linear_synteny, tree_node.top_down_synteny, tree_node.bottom_up_synteny)):
            data = []
            ev = AncestralSyntenyRels.columns['Evidence'].enum[evidence]
            print(f"process level {taxid} - graph {evidence} - |N|,|V| = {len(graph.nodes)},{len(graph.edges)}")

            for u, v, edge_data in graph.edges.data():
                w = edge_data["weight"]
                h1 = u.hog_id.rsplit('_')[0]
                h2 = v.hog_id.rsplit('_')[0]
                lca, orient, orient_score = -1, -1, 0
                if self.date_edges:
                    try:
                        lca_clade = edge_data['lca']
                        lca = self.tax2taxid[lca_clade]
                    except KeyError:
                        pass
                if self.orient_edges:
                    keys = ['unidirectional', 'divergent', 'convergent']
                    try:
                        orient_weights = numpy.array(list(map(lambda o: edge_data[o], keys)))
                    except KeyError:
                        pass
                    else:
                        orient = orient_enum[keys[numpy.argmax(orient_weights)]]
                        orient_score = orient_weights.max()
                rec = (hogid_lookup[h1], hogid_lookup[h2], w, ev, lca, orient, orient_score)
                data.append(rec)
            df = pd.DataFrame.from_records(numpy.array(data, dtype=dtype))
            dfs.append(df)
        df = pd.concat(dfs, ignore_index=True)
        sumdf = df.groupby(by=["HogRow1", "HogRow2"], as_index=False).agg({
            "Weight": "max",
            "Evidence": "min",
            # -1 indicates 'n/a', all other values should be consistent.
            "LCA_taxid": lambda x: x[x != -1].max() if not x[x != -1].empty else -1,
            "Orientation": lambda x: x[x != -1].max() if not x[x != -1].empty else -1,
            "OrientationScore": "max",

        })
        as_array = sumdf.astype({col: dtype.fields[col][0] for col in dtype.names}).to_records(index=False)
        tab = self.h5.create_table("/AncestralGenomes/tax{}".format(taxid),
                                   "Synteny", description=AncestralSyntenyRels,
                                   obj=as_array,
                                   expectedrows=len(as_array),
                                   createparents=True)
        for col in ("HogRow1", "HogRow2", "Weight", "Evidence", "LCA_taxid", ):
            tab.colinstances[col].create_csindex()
        print(f" hdf5 level {taxid} written: {len(tab)} rows.")

    def add_extant_graph(self, taxid, tree_node):
        if ExtantSyntenyRels is None:
            print(f"[WARNING] The installed pyoma library ({pyoma_version}) does not allow to store extant synteny results. Please update pyoma if needed")
            return

        dtype = tables.dtype_from_descr(ExtantSyntenyRels)
        if self.orient_edges:
            orient_enum = ExtantSyntenyRels.columns['Orientation'].enum

        gs = self.genome_lookup[taxid]
        os_code = gs['UniProtSpeciesCode'].decode()
        try:
            os = tree_node.genome_code
            if os != os_code:
                print(f"[WARNING] non-matching species codes: {os} vs {os_code}")
                return
        except AttributeError:
            pass

        # graph is input synteny graph stored in bottom_up_synteny
        graph = tree_node.bottom_up_synteny
        ev = ExtantSyntenyRels.columns['Evidence'].enum["linearized"]
        orient_enum = ExtantSyntenyRels.columns['Orientation'].enum

        print(f"process species {taxid} - graph |N|,|V| = {len(graph.nodes)},{len(graph.edges)}")
        data = []
        for u, v, edge_data in graph.edges.data():
            w = edge_data["weight"]
            if u.startswith(os_code):
                enr1 = int(u[len(os_code):]) + int(gs['EntryOff'])
                enr2 = int(v[len(os_code):]) + int(gs['EntryOff'])
            else:
                print(f"[WARNING] cannot parse protein ids: {u} / {v}")
                continue
            lca, orient, orient_score = -1, -1, 0
            if self.date_edges:
                try:
                    lca_clade = edge_data['lca']
                    lca = self.tax2taxid[lca_clade]
                except KeyError:
                    pass
            if self.orient_edges:
                keys = ['unidirectional', 'divergent', 'convergent']
                try:
                    orient_weights = numpy.array(list(map(lambda o: edge_data[o], keys)))
                except KeyError:
                    pass
                else:
                    orient = orient_enum[keys[numpy.argmax(orient_weights)]]
                    orient_score = orient_weights.max()
            rec = (enr1, enr2, w, ev, lca, orient, orient_score)
            data.append(rec)
        df = pd.DataFrame.from_records(numpy.array(data, dtype=dtype))
        sumdf = df.groupby(by=["EntryNr1", "EntryNr2"], as_index=False).agg({
            "Weight": "max",
            "Evidence": "min",
            # -1 indicates 'n/a', all other values should be consistent.
            "LCA_taxid": lambda x: x[x != -1].max() if not x[x != -1].empty else -1,
            "Orientation": lambda x: x[x != -1].max() if not x[x != -1].empty else -1,
            "OrientationScore": "max",
        })
        as_array = sumdf.astype({col: dtype.fields[col][0] for col in dtype.names}).to_records(index=False)

        tab = self.h5.create_table("/ExtantGenomes/{}".format(os_code),
                                   "Synteny", description=ExtantSyntenyRels,
                                   obj=as_array,
                                   expectedrows=len(as_array),
                                   createparents=True)
        for col in ("EntryNr1", "EntryNr2", "Weight", "Evidence", "LCA_taxid",):
            tab.colinstances[col].create_csindex()
        print(f" hdf5 extend {os_code} written: {len(tab)} rows.")

    def get_taxid_from_hog_names(self, graph):
        taxids = set([])
        for u in graph.nodes:
            try:
                hogid, taxid = u.hog_id.rsplit('_')
            except ValueError:
                continue
            taxids.add(taxid)
        if len(taxids) > 1:
            raise Exception("multiple taxids on this level: {}".format(taxids))
        elif len(taxids) == 0:
            return 0
        return taxids.pop()


def init_extant_graphs_from_hdf5(ham, hdf5_file, orient_edges):
    print('###################################')
    print('Initializing synteny graphs of extant genomes ...')
    h5 = db.Database(hdf5_file)
    ext_genomes = [Genome(h5, g) for g in h5.get_hdf5_handle().root.Genome.read()]
    for genome in ext_genomes:
        proteins = h5.main_isoforms(genome.uniprot_species_code)
        proteins.sort(order=['Chromosome', 'LocusStart'])
        graph = nx.Graph()
        contiguity_dict = dict()
        old_gene = ham.get_gene_by_id(proteins[0]['EntryNr'])
        old_gene_strand = proteins[0]['LocusStrand']
        old_contig = proteins[0]['Chromosome']
        graph.add_node(old_gene, contig = old_contig)
        for i in range(1, len(proteins)):
            gene = ham.get_gene_by_id(proteins[i]['EntryNr'])
            contig = proteins[i]['Chromosome']
            graph.add_node(gene, contig = contig)
            if contig == old_contig:
                graph.add_edge(gene, old_gene, weight=1, unidirectional=0, convergent=0, divergent=0, children = [None], extant_descendants = [None])
                if orient_edges:
                    old_gene_strand = proteins[i-1]['LocusStrand']
                    gene_strand = proteins[i]['LocusStrand']
                    try:
                        if old_gene_strand == gene_strand:
                            graph[old_gene][gene]['unidirectional'] = 1
                        elif old_gene_strand == 1:
                            graph[old_gene][gene]['convergent'] = 1
                        elif old_gene_strand == -1:
                            graph[old_gene][gene]['divergent'] = 1
                    except:
                        pass
            old_gene, old_contig = gene, contig
        gene.genome.taxon.add_feature('bottom_up_synteny', graph)
        gene.genome.taxon.add_feature('contiguity_dict', contiguity_dict)
        gene.genome.taxon.add_feature('genome_code', genome.uniprot_species_code)
    h5.close()
    return ham