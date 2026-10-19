"""Helpers to write small GFF files for the parser and extant graph tests."""


def gene(contig, start, end, strand, gene_id, *proteins, type="gene"):
    """A gene line followed by one CDS line per protein id."""
    rows = [(contig, type, start, end, strand, "ID=%s" % gene_id)]
    for i, prot in enumerate(proteins):
        rows.append((contig, "CDS", start, end, strand,
                     "ID=cds-%s-%d;Parent=%s;protein_id=%s" % (gene_id, i, gene_id, prot)))
    return rows


def write_gff(path, *blocks, header="##gff-version 3\n"):
    with open(path, "w") as fh:
        fh.write(header)
        for block in blocks:
            for contig, type_, start, end, strand, attrs in block:
                fh.write("\t".join([contig, "test", type_, str(start), str(end), ".", strand, ".", attrs]) + "\n")
    return str(path)
