# run concordance
# CHECK: {vcf_path}/snps_concordance/union_merged.bcftools.vcf.gz
# CHECK: {vcf_path}/snps_concordance/union_merged.bcftools.vcf.gz.csi
bcftools merge -O v `ls {vcf_path}/snps_concordance/*.union.bcftools.vcf.gz` | bcftools view -O z > {vcf_path}/snps_concordance/union_merged.bcftools.vcf.gz
bcftools index {vcf_path}/snps_concordance/union_merged.bcftools.vcf.gz
