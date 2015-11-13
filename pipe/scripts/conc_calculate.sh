# run concordance
# CHECK: {vcf_path}/snps_concordance/union_merged.bcftools.vcf.gz
# CHECK: {vcf_path}/snps_concordance/union_merged.bcftools.vcf.gz.csi
# CHECK: {results_path}/concordance/concordance.tsv
bcftools gtcheck -G 1 -H {vcf_path}/snps_concordance/union_merged.bcftools.vcf.gz > {results_path}/concordance/concordance.tsv
