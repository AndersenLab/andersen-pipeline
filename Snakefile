from subprocess import check_output
from utils import sample

localrules: download_picard
configfile: "config.yaml"
chrom = ["I","II","III","IV","V","X"]


#===============#
# Setup Samples #
#===============#


rule all:
    input:
        "log/setup_genome.done",
        "log/picard.done",
        "log/dir_setup.done"

rule setup_dirs:
    output:
        config["result_dir"],
        config["result_dir"] + "/" + config["bam_dir"],
        config["result_dir"] + "/" + config["vcf_dir"],
        "log/dir_setup.done"
    shell:
        """
            mkdir -p {config[result_dir]}
            mkdir -p {config[result_dir]}/{config[bam_dir]}
            mkdir -p {config[result_dir]}/{config[vcf_dir]}
            touch log/dir_setup.done
        """

rule download_picard:
    output:
        "tools/picard.jar",
        "log/picard.done"
    shell:
        """
            mkdir -p tools
            wget -O - https://github.com/broadinstitute/picard/releases/download/1.138/picard-tools-1.138.zip > tools/picard-tools-1.138.zip
            unzip -d tools tools/picard-tools-1.138.zip
            mv tools/picard-tools-1.138/* tools/
            rm -d tools/picard-tools-1.138.zip tools/picard-tools-1.138/
            touch log/picard.done
        """


include: "rules/setup_genome.snake"
#include: "rules/concordance_analysis.snake"
#include: "rules/align.snake"
#include: "rules/call.snake"


