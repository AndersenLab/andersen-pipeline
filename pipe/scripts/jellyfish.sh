# Jellyfish
# CHECK: {results_path}/jellyfish/{ID}.kmer.tsv
# CLEANUP: {results_path}/jellyfish/{ID}.tmp.fa
if [ "{debug}" = "True" ]; then
	cat {FQ1} {FQ2} > {results_path}/jellyfish/{ID}.tmp.fa
	cat={results_path}/jellyfish/{ID}.tmp.fa
else
	cat <(zcat {FQ1}) <(zcat {FQ2}) > {results_path}/jellyfish/{ID}.tmp.fa
fi

jellyfish count {command_jellyfish.params} --output={results_path}/jellyfish/{ID}.jf --threads={cores} {results_path}/jellyfish/{ID}.tmp.fa
jellyfish dump {results_path}/jellyfish/{ID}.jf | tr '\n>' '\t\n' | awk -v ID={ID} 'NR > 1{{print ID "\t" $2 "\t" $1}}' > {results_path}/jellyfish/{ID}.kmer.tsv

