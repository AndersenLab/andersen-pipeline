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
from utils import sample_file
from utils.reference import chunk_genome
import pipe
import sys
from pprint import pprint as pp
from subprocess import call

__version__ = "0.0.1"

def getScriptPath():
    return os.path.dirname(pipe.__file__)

def main():
    args = docopt(__doc__, 
                  version=__version__)
    print(args)

    if args["<config.yaml>"]:

        # Load configuration
        config = args["<config.yaml>"]
        if not os.path.exists(config):
            with indent(4):
                puts_err(colored.red("\n" + config + "not found\n"))
        config = yaml.load(open(args["<config.yaml>"], 'r'))
        fq_set = sample_file(config["sample_file"], config)

        if args["concordance"]:
            """
                For concordance analysis,
                the SM column becomes the ID column.
            """
            print fq_set
            for i in fq_set.records:
                i["SM"] = i["ID"]
            config["vcf_path"] = config["vcf_path"] + "_concordance"
            print pp(config)
            # create temporary files so as to only output a single readgroup... (merge into bam utility...)


        for i in fq_set.iterate_fq():
            print i
        for i in fq_set.iterate_SM_sets():
            print i
        for i in fq_set.iterate_bams():
            print i
        for i in fq_set.iterate_vcfs():
            print i
    elif args['genome']:
        comm = ['python', getScriptPath() + '/' + "genome.py"] + sys.argv[1:]
        print comm
        exit(call(comm))

