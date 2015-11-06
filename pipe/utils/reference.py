

def chunk_genome(chunk_size, reference):
    """ 
    Parses bwa .ann file to retrieve chromosome sizes
    for chunking purposes
    """
    ann = open(reference + ".ann").read()
    # Parsing .ann files
    contigs = [x.split(" ")[1] for x in ann.split("\n")[1:-1:1]][::2]
    contig_sizes = map(int,[x.split(" ")[1] for x in ann.split("\n")[1:-1:1]][1::2])
    for chrom, size in zip(contigs, contig_sizes):
        for chunk in xrange(1,size, chunk_size):
            if chunk + chunk_size > size:
                chunk_end = size
            else:
                chunk_end = chunk + chunk_size-1
            yield "{chrom}:{chunk}-{chunk_end}".format(**locals())

