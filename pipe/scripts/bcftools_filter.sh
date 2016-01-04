# bcftools_filter
# CHECK: {vcf_path}/snps/{SM}.union.filtered.bcftools.vcf.gz
# CHECK: {vcf_path}/snps/{SM}.union.filtered.bcftools.vcf.gz.csi

min_depth=8
qual=30
mq=40
dv_dp=0.5

function filter_variants {{
    bcftools view $1 | \
    bcftools filter --mode + --soft-filter quality --include "QUAL >= ${{qual}}" |  \
    bcftools filter --mode + --soft-filter min_depth --include "INFO/DP > ${{min_depth}}" | \
    bcftools filter --mode + --soft-filter mapping_quality --include "INFO/MQ > ${{mq}}" | \
    vcffixup - | \
    bcftools filter --mode + --soft-filter dv_dp --include "(FORMAT/DV)/(FORMAT/DP) >= ${{dv_dp}} || FORMAT/GT == '0/0'" | \
    bcftools filter --mode + --soft-filter het --exclude 'AC==1' | \
    tb geno transfer-filter -
}}

filter_variants {vcf_path}/snps/{SM}.union.bcftools.vcf.gz | bcftools view -O z > {vcf_path}/snps/{SM}.union.filtered.bcftools.vcf.gz
bcftools index {vcf_path}/snps/{SM}.union.filtered.bcftools.vcf.gz
