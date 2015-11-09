# Coverage
# CHECK: {results_path}/eav/{SM}.global.coverage.tsv
# CHECK: {results_path}/eav/{SM}.100000.coverage.tsv
# CHECK: {results_path}/eav/{SM}.10000.coverage.tsv
# CHECK: {results_path}/eav/{SM}.talt.coverage.tsv
bam coverage {bam_merged} > {results_path}/eav/{SM}.global.coverage.tsv
bam coverage {bam_merged} --window=100000 > {results_path}/eav/{SM}.100000.coverage.tsv
bam coverage {bam_merged} --window=10000 > {results_path}/eav/{SM}.10000.coverage.tsv
bam coverage {bam_merged} --regions={analysis_path}/extra/talt_regions.WS245.bed > {results_path}/eav/{SM}.talt.coverage.tsv
