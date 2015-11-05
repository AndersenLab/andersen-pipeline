"""Andersen Pipeline

Usage:
  pipe run <config.yaml>
  pipe genome
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

__version__ = "0.0.1"

def main():
    args = docopt(__doc__, version=__version__)
    print(args)
    if args["run"]:
        # Load configuration
        config = args["<config.yaml>"]
        if not os.path.exists(config):
            with indent(4):
              puts_err(colored.red("\n" + config + "not found\n"))
        config = yaml.load(open(args["<config.yaml>"], 'r'))
        fq_set = sample_file(config["sample_file"], config)
        for i in fq_set.iterate_fq():
            print i
        for i in fq_set.iterate_SM_sets():
            print i 
        for i in fq_set.iterate_bams():
            print i 
        for i in fq_set.iterate_vcfs():
            print i 

