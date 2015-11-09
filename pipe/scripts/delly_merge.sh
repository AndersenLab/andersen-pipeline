# delly merge
# CHECK: {vcf_path}/sv/delly.SV.vcf.gz
# CHECK: {vcf_path}/sv/delly.SV.vcf.gz.csi
bcftools merge -O z `ls {results_path}/sv/delly.*.vcf.gz` > {vcf_path}/sv/delly.SV.vcf.gz
bcftools index {vcf_path}/sv/delly.SV.vcf.gz

file={results_path}/sv/delly.SV.vcf.gz.csi
minimumsize=10
actualsize=$(wc -c <"$file")
if [ $actualsize -ge $minimumsize ]; then
    rm {results_path}/sv/delly.*.vcf.gz
fi