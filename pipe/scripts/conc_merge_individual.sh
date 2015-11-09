# union_snp_sites
# CHECK: {vcf_path}/snps_concordance/merged_sites.tsv
cat `ls {vcf_path}/snps_concordance/*.individual.sites.tsv` | sort -k1,1 -k2,2n | uniq > {vcf_path}/snps_concordance/merged_sites.tsv
rm {vcf_path}/snps_concordance/*.individual.sites.tsv
