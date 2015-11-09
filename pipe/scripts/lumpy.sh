# lumpy
lumpyexpress -B {bam_merged} -T {vcf_path}/lumpy.{SM}.sv.tmp -o {vcf_path}/sv/lumpy.{SM}.sv.vcf
bcftools view -O z {vcf_path}/sv/lumpy.{SM}.sv.vcf > {vcf_path}/sv/lumpy.{SM}.sv.vcf.gz
bcftools index {vcf_path}/sv/lumpy.{SM}.sv.vcf.gz

file={vcf_path}/sv/lumpy.{SM}.sv.vcf.gz.csi
minimumsize=10
actualsize=$(wc -c <"$file")
if [ $actualsize -ge $minimumsize ]; then
    rm {vcf_path}/sv/lumpy.{SM}.sv.vcf
fi
