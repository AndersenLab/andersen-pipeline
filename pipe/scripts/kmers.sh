# fastq-kmers
# CHECK: {results_path}/kmers/{ID}.kmer.tsv
zcat {FQ1} {FQ2} | fastq-kmers {command_kmers.params} | \
awk -v ID={ID} 'NR == 1 {{ print "ID\t" $0}} NR > 1 {{print ID "\t" $0}}' > {results_path}/kmers/{ID}.kmer.tsv