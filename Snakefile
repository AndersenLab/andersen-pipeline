from bunch import *
from pprint import pprint as print
import glob
from snakemake.utils import makedirs
import shutil
import os
import sys
print(sys.path)
from parsers import *



fastq_dir = "/Users/dancook/coding/git/pipeline/fq"

fq_set = ["BGI2-RET2-test1-paired", "EA-D06_CGCGCAG_L001"]

fq_list = expand("{sample}" + "-1", sample=fq_set) + \
          expand("{sample}" + "-2", sample=fq_set)


eav_file = "eav.txt"

dirs = Bunch()
dirs.logs = "logs"
dirs.results = "results"
dirs.fastqc = "results/fastqc"

fastqc_image_reports = ["adapter_content",
                        "duplication_levels",
                        "per_base_n_content",
                        "per_base_quality",
                        "per_base_sequence_content",
                        "per_sequence_gc_content",
                        "per_sequence_quality",
                        "per_tile_quality",
                        "kmer_profiles",
                        "sequence_length_distribution"]

rule all:
    input: "done.txt"

rule setup_directories:
    output:
        dirs.values(),
        expand("{qc_dir}/{img_dir}", qc_dir = dirs.fastqc, img_dir = fastqc_image_reports )
    message: "Creating Directories {output}"
    run:
        makedirs(output)

rule run_fastqc:
    input: 
        rules.setup_directories.output,
        fq_list = expand("{fastq_dir}/{fq_list}.fq.gz", fastq_dir = fastq_dir, fq_list = fq_list)
    log: "fastqc.log"
    output: 
        files = expand("{fastqc_dir}/{fq_list}.fq_fastqc/{files}.txt", fastqc_dir = dirs.fastqc,
                                                                       fq_list=fq_list,
                                                                       files = ["fastqc_data"]),
        zip_html = temp(expand("{fastqc_dir}/{fq_list}.fq_fastqc.{ext}", fastqc_dir = dirs.fastqc,
                                                                         fq_list=fq_list,
                                                                         ext = ["zip", "html"])),
    threads: 2
    version: shell("fastqc --version")
    shell: 
        """
            fastqc --outdir {dirs.fastqc} --extract --threads {threads} {input.fq_list}
            for f in {output.files} {output.zip_html}; do
                touch `ls $f`;
            done;
        """

rule process_fastqc:
    input: rules.run_fastqc.output.files
    message: "Moving Files"
    output:
        "done.txt"
    run:
        for fq in input:
            print(fq)
            fastqc(fq).to_eav()
        fastqc_dir = dirs.fastqc
        for fq in fq_list:
            # Save data to eav.
            for image in glob.glob("{fastqc_dir}/{fq}.fq_fastqc/Images/*.png".format(**locals())):
                img_report = os.path.basename(image).replace(".png","")
                mv_loc = "{fastqc_dir}/{img_report}/{fq}.png".format(**locals())
                shutil.move(image, mv_loc)
            #shutil.rmtree("{fastqc_dir}/{fq}.fq_fastqc/".format(**locals()))
        shell("touch done.txt")




