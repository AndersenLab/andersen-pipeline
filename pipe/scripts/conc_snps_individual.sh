# snps individual
# CHECK: {vcf_path}/snps_concordance/merged_sites.tsv
contigs=`samtools view -H {bam_merged} | grep -Po 'SN:(.*)\t' | cut -c 4-1000`
parallel --gnu -j {cores} 'samtools mpileup -r {{}} {command_samtools.params} -G <(bam readgroups {bam_merged} --include={ID}) --fasta-ref {reference} {bam_merged} | bcftools call --skip-variants indels -O z {command_bcftools.params} - > {vcf_path}/snps_concordance/{ID}.{{}}.individual.bcftools.vcf.gz' ::: ${{contigs}}
order=`echo $contigs | tr ' ' '\n' | awk -v vcf_path=/lscr2/andersenlab/dec211/WI_debug/vcf/snps_concordance/{ID} '{{ print vcf_path "." $1".individual.bcftools.vcf.gz" }}'`
bcftools concat ${{order}} -O v | bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' > {vcf_path}/snps_concordance/{ID}.individual.sites.tsv
rm ${{order}}

