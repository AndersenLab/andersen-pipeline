# Perform alignment and sort
# CHECK: {bam_merged}
# CHECK: {bam_merged}.bai
sambamba merge --nthreads={command_bwa.cores} --show-progress {bam_merged} {bamlist}
sambamba index --nthreads={command_bwa.cores} {bam_merged}


file={bam_merged}.bai
minimumsize=10
actualsize=$(wc -c <"$file")
if [ $actualsize -ge $minimumsize ]; then
    for i in {bamlist}; do
    	rm $i*
    done;
fi
