import pipe
import os
from itertools import chain
from collections import namedtuple
from easydict import EasyDict as edict
from slurmpy import Slurm
from pprint import pprint as pp

def format_params(params):
    args, flags = "", ""
    if 'args' in params.keys():
        args = ' '.join(list(chain(*[[k,str(v)] for k,v in params["args"].items()])))
    if 'flags' in params.keys():
        flags = ' '.join(params["flags"])
    return args + " " + flags


def run_script(script_name, name, args, config, dependencies = "", node = ""):
    args.update(config)
    # Format Parameters
    for k,v in args.items():
        if k.startswith("command_") and type(v) != bool:
            args[k]["params"] = format_params(args[k])
    script = os.path.dirname(pipe.__file__) + "/scripts/" + script_name + ".sh"
    script = open(script, 'r').read()
    # Convert to DotMap for attribute access 
    args = edict(args)
    # Configure options:
    slurm_kwargs = {}
    if node:
        slurm_kwargs["node"] = node
    if dependencies:
        dependencies = ':'.join([str(x.job_id) for x in dependencies])
        slurm_kwargs["dependency"] = "afterok:" + dependencies
    print "Submitted " + script_name + "-" + name + "; " + dependencies 
    with open("long_script.sh", 'a+') as f:
        f.write(script.format(**args))
    print script.format(**args)
    return Slurm(script_name + "-" + name, slurm_kwargs).run(script.format(**args))

    
