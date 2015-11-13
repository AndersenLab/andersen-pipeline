# conc snps union
# CHECK: {vcf_path}/snps_concordance/{ID}.union.bcftools.vcf.gz
# CHECK: {vcf_path}/snps_concordance/{ID}.union.bcftools.vcf.gz.csi
# CHECK: {vcf_path}/snps_concordance/merged_sites.tsv
contigs=`samtools view -H {bam_merged} | grep -Po 'SN:(.*)\t' | cut -c 4-1000`
parallel --gnu -j {cores} 'samtools mpileup -r {{}} {command_samtools.params} -G <(bam readgroups {bam_merged} --include={ID}) --fasta-ref {reference} {bam_merged} | bcftools call --skip-variants indels -m -T {vcf_path}/snps_concordance/merged_sites.tsv -O v - | tb rename --subst={SM}:{ID} - | tb geno het-polarization - | bcftools view -O z > {vcf_path}/snps_concordance/{ID}.{{}}.union.bcftools.vcf.gz' ::: ${{contigs}}
order=`echo $contigs | tr ' ' '\n' | awk -v vcf_path={vcf_path}/snps_concordance/{ID} '{{ print vcf_path "." $1".union.bcftools.vcf.gz" }}'`
bcftools concat ${{order}} -O z > {vcf_path}/snps_concordance/{ID}.union.bcftools.vcf.gz
bcftools index {vcf_path}/snps_concordance/{ID}.union.bcftools.vcf.gz
rm ${{order}}

