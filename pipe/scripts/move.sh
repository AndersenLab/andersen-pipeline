# Move and index
mv {bamlist} {bam_merged}
sambamba index --nthreads={command_bwa.cores} {bam_merged}


file={bam_merged}.bai
minimumsize=10
actualsize=$(wc -c <"$file")
if [ $actualsize -ge $minimumsize ]; then
    for i in {bamlist}; do
    	rm $i*
    done;
fi
