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
from subprocess import call

__version__ = "0.0.1"

def getScriptPath():
    return os.path.dirname(pipe.__file__)

def main():
    args = docopt(__doc__, 
                  version=__version__)
    if args["<config.yaml>"]:

        # Load configuration
        config = check_file(args["<config.yaml>"])
        config = yaml.load(open(config, 'r'))
        fq_set = sample_file(config["sample_file"], config)


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
        for SM in fq_set.iterate_SM_sets():
            for fqs in SM["fq_set"]:
                print get_script("align", fqs, config)
    elif args['genome']:
        comm = ['python', getScriptPath() + '/' + "genome.py"] + sys.argv[1:]
        print comm
        exit(call(comm))

