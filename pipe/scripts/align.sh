# Perform alignment and sort
bwa {command_bwa.type} {command_bwa.params} -R '{readgroup}' {reference} {FQ1} {FQ2} | \
sambamba view --nthreads={command_bwa.cores} --sam-input --format=bam --with-header /dev/stdin | \
sambamba sort --nthreads={command_bwa.cores} --show-progress --tmpdir={tmp_path} --nthreads={command_bwa.cores} --out={bam_individual}.sorted.bam /dev/stdin

# Mark duplicates
picard MarkDuplicates I={bam_individual}.sorted.bam O={bam_individual}.bam M={bam_individual}.dups.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
samtools index {bam_individual}.bam
mv {bam_individual}.dups.txt {results_path}/dups/{ID}.dups.txt
rm {bam_individual}.sorted.bam*

