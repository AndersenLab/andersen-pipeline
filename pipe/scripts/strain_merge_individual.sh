# bcftools merge individual
# CHECK: {vcf_path}/merged_sites.tsv
cat `ls {vcf_path}/strain_snps/*.individual.sites.tsv` | sort -k1,1 -k2,2n | uniq > {vcf_path}/strain_merged_sites.tsv
#rm {vcf_path}/strain_snps/*.individual.sites.tsv
