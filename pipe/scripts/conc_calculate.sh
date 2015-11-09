# run concordance
# CHECK: {vcf_path}/snps_concordance/union_merged.bcftools.vcf.gz
# CHECK: {results_path}/concordance/concordance.tsv
bcftools merge -O v `ls {vcf_path}/snps_concordance/*.union.bcftools.vcf.gz` | tb geno het-polarization - | bcftools view -O z > {vcf_path}/snps_concordance/union_merged.bcftools.vcf.gz
bcftools index {vcf_path}/snps_concordance/union_merged.bcftools.vcf.gz
bcftools gtcheck -G 1 -H {vcf_path}/snps_concordance/union_merged.bcftools.vcf.gz > {results_path}/concordance/concordance.tsv
