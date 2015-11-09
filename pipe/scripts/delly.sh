# delly
# CHECK: {vcf_path}/sv/delly.{svtype}.vcf.gz
# CHECK: {vcf_path}/sv/delly.{svtype}.vcf.gz.csi
delly -g {reference} --type {svtype} --outfile {vcf_path}/sv/delly.{svtype}.vcf.gz {merged_bam_list}
bcftools view -O z {results_path}/sv/delly.{svtype}.vcf  > {vcf_path}/sv/delly.{svtype}.vcf.gz
bcftools index {vcf_path}/sv/delly.{svtype}.vcf.gz


file={results_path}/sv/delly.{svtype}.vcf.gz.csi
minimumsize=10
actualsize=$(wc -c <"$file")
if [ $actualsize -ge $minimumsize ]; then
    rm {results_path}/sv/delly.{svtype}.vcf
fi
