# bcftools
contigs=`samtools view -H {bam_merged} | grep -Po 'SN:(.*)\t' | cut -c 4-1000`
parallel --gnu -j {cores} 'samtools mpileup -r {{}} {command_samtools.params} --fasta-ref {reference} {bam_merged} | tb geno het-polarization - | bcftools call -O z {command_bcftools.params} - > {vcf_path}/snps/{SM}.{{}}.individual.bcftools.vcf.gz' ::: ${{contigs}}
order=`echo $contigs | tr ' ' '\n' | awk -v vcf_path={vcf_path}/snps/{SM} '{{ print vcf_path "." $1".individual.bcftools.vcf.gz" }}'`

bcftools concat ${{order}} -O z > {vcf_path}/snps/{SM}.individual.bcftools.vcf.gz
bcftools index {vcf_path}/snps/{SM}.individual.bcftools.vcf.gz

bcftools concat ${{order}} -O v | \
bcftools view -M 2 -m 2 -O v - | \
tb geno het-polarization - | \
bcftools filter --include 'DP > 3' | \
egrep '(^#|1/1)' | \
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' > {vcf_path}/snps/{SM}.individual.sites.tsv