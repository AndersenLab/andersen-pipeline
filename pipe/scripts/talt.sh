# Coverage
# CHECK: {results_path}/talt/{SM}.talt.coverage.tsv
bam coverage {bam_merged} --tsv --regions={analysis_path}/extra/talt_regions.WS245.control.bed > {results_path}/talt/{SM}.talt.coverage.tsv
bam coverage {bam_merged} --tsv --regions={analysis_path}/extra/talt_regions.WS245.bed > {results_path}/talt/{SM}.talt.coverage.tsv
