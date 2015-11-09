# Fastq
# CHECK: {results_path}/eav/{ID}-1.tsv
# CHECK: {results_path}/eav/{ID}-2.tsv
bam fastq {FQ1} > {results_path}/eav/{ID}-1.tsv
bam fastq {FQ2} > {results_path}/eav/{ID}-2.tsv
