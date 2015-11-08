"""Andersen Pipeline

Usage:
  pipe genome [--search=<search> --download=<setup>]
  pipe <config.yaml>
  pipe concordance <config.yaml>

Options:
  -h --help     Show this screen.
  --version     Show version.

"""
from docopt import docopt
import yaml
from clint.textui import puts_err, colored, indent
import os
from utils import *
from utils.reference import chunk_genome
from commands import *
import pipe
import sys
from pprint import pprint as pp
import datetime
from subprocess import call
from slurmpy import node_cycle

__version__ = "0.0.1"

def getScriptPath():
    return os.path.dirname(pipe.__file__)

def main():
    args = docopt(__doc__,
                  version=__version__)
    if args["<config.yaml>"]:
        with open("long_script.sh", 'w+') as f:
            f.write("#" + str(datetime.datetime.now()) + "\n")

        # Load configuration
        config = check_file(args["<config.yaml>"])
        config = yaml.load(open(config, 'r'))
        fq_set = sample_file(config["sample_file"], config)
        nodes = node_cycle(config["nodelist"])

        if args["concordance"]:
            """
                For concordance analysis,
                the SM column becomes the ID column.
            """
            for i in fq_set.records:
                i["SM"] = i["ID"]
            config["vcf_path"] = config["vcf_path"] + "_concordance"
            # create temporary files so as to only output a single readgroup... (merge into bam utility...)

        #===========#
        # Alignment #
        #===========#
        merged_bam_list = []
        for SM in fq_set.iterate_SM_sets():
            deplist = []
            merged_deplist = []
            for fqs in SM["fq_set"]:
                # Jellyfish
                #if "command_jellyfish" in config:
                #    run_script("jellyfish", fqs["ID"], fqs, config)

                # Alignment
                dep = run_script("align", fqs["ID"], fqs, config)
                deplist.append(dep)
            if len(SM["fq_set"]) > 1:
                dep = run_script("merge", SM["SM"], SM, config, dependencies = deplist)
            else:
                dep = run_script("move", SM["SM"], SM, config, dependencies = deplist)
            merged_deplist.append(dep)

            # Run Telseq
            if "command_telseq" in config:
                run_script("telseq", SM["SM"], SM, config, dependencies = [dep])

            # Run Lumpy
            if "command_lumpy" in config:
                run_script("lumpy", SM["SM"], SM, config, dependencies = [dep])

            # Run bcftools
            if "command_bcftools" in config:
                run_script("bcftools", SM["SM"], SM, config, dependencies = [dep])

            # Run Coverage
            run_script("coverage", SM["SM"], SM, config, dependencies = [dep])

        for bam in fq_set.iterate_bams():
            merged_bam_list.append(bam)

        # Run delly
        if "command_delly" in config:
            delly_dep = []
            for svtype in ["DEL","DUP","INV","TRA","INS"]:
                delly_dep.append(run_script("delly",
                                 SM["SM"],
                                 {"svtype": svtype, "merged_bam_list": ' '.join(merged_bam_list)},
                                 config,
                                 dependencies = merged_deplist))
            run_script("delly_merge", SM["SM"], {}, config, dependencies = delly_dep)

        # Run 

    elif args['genome']:
        comm = ['python', getScriptPath() + '/' + "genome.py"] + sys.argv[1:]
        print comm
        exit(call(comm))
