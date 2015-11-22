# bcftools_filter
# CHECK: {vcf_path}/snps_concordance/{ID}.union.filtered.bcftools.vcf.gz
# CHECK: {vcf_path}/snps_concordance/{ID}.union.filtered.bcftools.vcf.gz.csi

min_depth=3
mq0f=5
qual=30

function max_depth {{
    bcftools query -f '%DP\n' $1 | awk 'BEGIN {{ x = 0; }} {{ x = x + $0}} END {{ print (x/NR)*1.9 }}'
}}

function filter_variants {{
    max=`max_depth $1`
    max_depth_filter="INFO/DP<${{max}}"
    bcftools view $1 | \
    bcftools filter --mode + --soft-filter quality --include "QUAL >= ${{qual}}" |  \
    bcftools filter --mode + --soft-filter min_depth --include "INFO/DP > ${{min_depth}}" | \
    bcftools filter --mode + --soft-filter max_depth --include "${{max_depth_filter}}" | \
    bcftools filter --mode + --soft-filter mq0f_lt_5 --include "MQ0F <= ${{mq0f}}" | \
    tb geno transfer-filter -
}}

filter_variants {vcf_path}/snps_concordance/{ID}.union.bcftools.vcf.gz | bcftools view -O z > {vcf_path}/snps_concordance/{ID}.union.filtered.bcftools.vcf.gz
bcftools index {vcf_path}/snps_concordance/{ID}.union.filtered.bcftools.vcf.gz
