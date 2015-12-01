# bcftools merge union
# CHECK: {vcf_path}/union_merged.vcf.gz
# CHECK: {vcf_path}/union_merged.vcf.gz.csi
# CHECK: {vcf_path}/union_merged.snpeff.vcf.gz
# CHECK: {vcf_path}/union_merged.snpeff.vcf.gz.csi
# CHECK: {vcf_path}/union_merged.clean.vcf.gz


source ~/.env
# Merge SNPs
parallel_bcftools_merge -m all `ls {vcf_path}/snps/*.union.filtered.bcftools.vcf.gz` | \
bcftools view -O v -  | \
tb filter REF --min=1 - | \
tb filter ALT --min=1 - | \
tb filter ALT --max=0.99 - | \
bcftools view -O z - > {vcf_path}/union_merged.vcf.gz
bcftools index -f {vcf_path}/union_merged.vcf.gz

# Variant Prediction
bcftools view -O v {vcf_path}/union_merged.vcf.gz | snpEff eff WS241 | bcftools view -O z > {vcf_path}/union_merged.snpeff.vcf.gz
bcftools index -f {vcf_path}/union_merged.snpeff.vcf.gz

# Generate imputation vcf
bcftools view -m 2 -M 2 --types snps {vcf_path}/union_merged.vcf.gz | \
bcftools filter --set-GTs . --exclude 'FORMAT/GF != "PASS"' | \
tb filter MISSING --max=0.10 - | \
tb filter HET --max=0.10 - | \
bcftools view -O z > {vcf_path}/union_merged.clean.vcf.gz
bcftools index {vcf_path}/union_merged.clean.vcf.gz

# Impute
# java -Xmx32768m -jar ~/.linuxbrew/bin/beagle.jar nthreads={cores} gt={vcf_path}/union_merged.clean.vcf.gz out={vcf_path}/union_merged.impute
# bcftools index -f {vcf_path}/union_merged.impute.vcf.gz
# tabix {vcf_path}/union_merged.impute.vcf.gz
# 
# # Kinship and mapping
# Rscript -e 'library(cegwas); kinship <- generate_kinship("{vcf_path}/union_merged.impute.vcf.gz"); save(kinship, file = "{vcf_path}/kinship.Rda"); '
# Rscript -e 'library(cegwas); snps <- generate_mapping("{vcf_path}/union_merged.impute.vcf.gz"); save(snps, file = "{vcf_path}/snps.Rda"); '# 