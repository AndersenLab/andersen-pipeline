# bcftools
# CHECK: {vcf_path}/snps/{SM}.union.bcftools.vcf.gz
# CHECK: {vcf_path}/snps/{SM}.union.bcftools.vcf.gz.csi
# CHECK: {vcf_path}/snps/merged_sites.tsv
contigs=`samtools view -H {bam_merged} | grep -Po 'SN:(.*)\t' | cut -c 4-1000`
parallel --gnu -j {cores} 'samtools mpileup -r {{}} {command_samtools.params} --fasta-ref {reference} {bam_merged} | bcftools call --insert-missed -T {vcf_path}/snps/merged_sites.tsv --skip-variants indels -O v {command_bcftools.params} - | tb geno het-polarization - | bcftools view -O z > {vcf_path}/snps/{SM}.{{}}.union.bcftools.vcf.gz' ::: ${{contigs}}
order=`echo $contigs | tr ' ' '\n' | awk -v vcf_path={vcf_path}/snps_concordance/{SM} '{{ print vcf_path "." $1".union.bcftools.vcf.gz" }}'`
bcftools concat ${{order}} -O z > {vcf_path}/snps/{SM}.union.bcftools.vcf.gz
bcftools index {vcf_path}/snps/{SM}.union.bcftools.vcf.gz
rm ${{order}}

