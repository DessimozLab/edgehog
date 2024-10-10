import numpy
import pandas as pd
import sys
import networkx as nx

try:
    import tables
    from pyoma.browser.tablefmt import AncestralSyntenyRels
    from pyoma.browser import db
    from pyoma.browser.models import Genome
    from pyoma import version as pyoma_version
except ImportError:
    print(f"pyoma and pytables libraries are required to work on HDF5 files. "
          f"Please install edgehog with the `oma` extras added, i.e. `pip install edgehog[oma]`.")
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
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.h5.close()
        self.oma_h5.close()

    def _load_taxname_2_taxid(self):
        tab = self.oma_h5.get_node("/Taxonomy")
        return {row['Name'].decode(): int(row['NCBITaxonId']) for row in tab.read()}

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
        as_array = sumdf.to_records(index=False, column_dtypes={"LCA_taxid": numpy.int32})
        tab = self.h5.create_table("/AncestralGenomes/tax{}".format(taxid),
                                   "Synteny", description=AncestralSyntenyRels,
                                   obj=as_array,
                                   expectedrows=len(as_array),
                                   createparents=True)
        for col in ("HogRow1", "HogRow2", "Weight", "Evidence", "LCA_taxid", ):
            tab.colinstances[col].create_csindex()
        print(f" hdf5 level {taxid} written: {len(tab)} rows.")

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
    h5.close()
    return ham