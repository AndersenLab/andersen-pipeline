# run concordance
# CHECK: {vcf_path}/snps_concordance/union_merged.bcftools.vcf.gz
# CHECK: {vcf_path}/snps_concordance/union_merged.bcftools.vcf.gz.csi
# CHECK: {results_path}/concordance/concordance.tsv

## Filter SP
parallel --gnu -j {cores} 'bcftools view {vcf_path}/snps_concordance/union_merged.bcftools.vcf.gz {{}} | \
bcftools filter --soft-filter "SP_filter" --mode + --exclude "FORMAT/SP > 0" --set-GTs . | \
bcftools gtcheck -G 1 -H - | \
awk -v chrom={{}} "/^CN/ {{ print \$0 \"\t\" chrom }} \$0 !~ /.*CN.*/ {{ print }} \$0 ~ /^# \[1\]CN/ {{ print \$0 \"\tchrom\"}}" - > {results_path}/concordance/{{}}.concordance.tsv' ::: I II III IV V X MtDNA

cat `find {results_path}/concordance/ | egrep '(I|II|III|IV|V|X|MtDNA)\.concordance\.tsv'` | \
grep 'CN' | \
awk 'NR == 1 && /Discordance/ {{ print }} NR > 1 && $0 !~ /Discordance/ {{ print }}' | \
awk '{{ gsub("(# |\\[[0-9]+\\])","", $0); gsub(" ","_", $0); print }}' | \
cut -f 2-7 > {results_path}/concordance/concordance_raw.tsv

datamash --headers --sort --group=Sample_i,Sample_j sum Discordance  sum Number_of_sites  < {results_path}/concordance/concordance_raw.tsv > {results_path}/concordance/concordance.tsv
