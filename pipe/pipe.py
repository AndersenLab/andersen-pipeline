"""Andersen Pipeline

Usage:
  pipe genome [--search=<search> --download=<setup>]
  pipe <config.yaml>

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

        #===========#
        # Alignment #
        #===========#
        merged_bam_list = []
        concordance_deplist = []
        eav_deplist = []
        telseq_deplist = []
        for SM in fq_set.iterate_SM_sets():
            all_id_bams_ok = [] 
            all_sm_bams_ok = []
            for args in SM["fq_set"]:
                fq_script_args = { "name" : args["ID"],
                                   "args" : args,
                                   "config" : config}
                # Kmer programs
                run_script("jellyfish", **fq_script_args)
                run_script("kmers", **fq_script_args)

                # fastq analysis
                eav_deplist.append(run_script("fastq", **fq_script_args))
                
                # Alignment
                all_id_bams_ok.append(run_script("align", **fq_script_args))

            if len(SM["fq_set"]) > 1:
                SM_merged = run_script("merge", SM["SM"], SM, config, dependencies = all_id_bams_ok)
            else:
                SM_merged = run_script("move", SM["SM"], SM, config, dependencies = all_id_bams_ok)
            all_sm_bams_ok.append(SM_merged)
            #==========================#
            # Concordance Analysis (1) #
            #==========================#
            for fqs in SM["fq_set"]:
                # Generate individual snpset
                concordance_deplist.append(run_script("conc_snps_individual", fqs["ID"], fqs, config, dependencies = [SM_merged]))

            # Run Telseq
            telseq_deplist.append(run_script("telseq", SM["SM"], SM, config, dependencies = [SM_merged]))

            # Run Lumpy
            run_script("lumpy", SM["SM"], SM, config, dependencies = [SM_merged])

            # Run bcftools
            run_script("bcftools", SM["SM"], SM, config, dependencies = [SM_merged])

            # Run Coverage
            eav_deplist.append(run_script("coverage", SM["SM"], SM, config, dependencies = [SM_merged]))

        #==========================#
        # Concordance Analysis (2) #
        #==========================#
        snps_merged_dep = run_script("conc_merge_individual", "merge", fqs, config, dependencies = concordance_deplist)
        conc_union_deps = []
        for SM in fq_set.iterate_SM_sets():
            for fqs in SM["fq_set"]:
                conc_union_deps.append(run_script("conc_snps_union", fqs["ID"], fqs, config, dependencies = all_sm_bams_ok))
        run_script("conc_calculate", fqs["ID"], fqs, config, dependencies = conc_union_deps)

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
                                 dependencies = all_sm_bams_ok))
            run_script("delly_merge", SM["SM"], {}, config, dependencies = delly_dep)

        # Run aggregation
        run_script("aggregate_telseq", "ALL", {}, config, dependencies = telseq_deplist)
        run_script("aggregate_eav", "ALL", {}, config, dependencies = eav_deplist)

    elif args['genome']:
        comm = ['python', getScriptPath() + '/' + "genome.py"] + sys.argv[1:]
        print comm
        exit(call(comm))
