# run concordance
# CHECK: {vcf_path}/snps_concordance/union_merged.bcftools.vcf.gz
# CHECK: {vcf_path}/snps_concordance/union_merged.bcftools.vcf.gz.csi
bcftools merge -m all --apply-filters PASS -O v `ls {vcf_path}/snps_concordance/*.union.filtered.bcftools.vcf.gz` | \
bcftools view --exclude 'FORMAT/SP > 0' --set-GTs . | \
bcftools view -X -O z > {vcf_path}/snps_concordance/union_merged.bcftools.vcf.gz
bcftools index {vcf_path}/snps_concordance/union_merged.bcftools.vcf.gz
