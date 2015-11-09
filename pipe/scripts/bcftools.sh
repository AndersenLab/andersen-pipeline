# bcftools
contigs=`samtools view -H {bam_merged} | grep -Po 'SN:(.*)\t' | cut -c 4-1000`
parallel --gnu -j {cores} 'samtools mpileup -r {{}} {command_samtools.params} --fasta-ref /lscr2/andersenlab/dec211/pyPipeline/genomes/WS245/c_elegans.PRJNA13758.WS245.genomic.fa.gz {bam_merged} | bcftools call -O z {command_bcftools.params} - > {vcf_path}/snps/{SM}.{{}}.bcftools.vcf.gz' ::: ${{contigs}}
order=`echo $contigs | tr ' ' '\n' | awk -v vcf_path=/lscr2/andersenlab/dec211/WI_debug/vcf/snps/{SM} '{{ print vcf_path "." $1".bcftools.vcf.gz" }}'`
bcftools concat ${{order}} -O z > {vcf_path}/snps/{SM}.bcftools.vcf.gz
bcftools index {vcf_path}/snps/{SM}.bcftools.vcf.gz
rm ${{order}}

