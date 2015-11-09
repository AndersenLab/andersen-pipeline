# Telseq
# CHECK: {results_path}/telseq/{SM}.telseq.txt
samtools view -h {bam_merged} | \
awk '/^@/  {{ print }} {{if (length($10) > 80) {{ print }}}}'  | \
samtools view -b - | telseq {command_telseq.params} - | \
awk -v s={SM} '{{ print s "\t" $0}}' | \
egrep -v '\[-\]|BAMs' > {results_path}/telseq/{SM}.telseq.txt
