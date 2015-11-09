import pipe
import os
from itertools import chain
from collections import namedtuple
from easydict import EasyDict as edict
from slurmpy import Slurm
from pprint import pprint as pp
from subprocess import check_output

def file_exists(filename):
    if os.path.isfile(filename) and os.path.getsize(filename) > 0:
        return True
    else:
        return False

def format_params(params):
    args, flags = "", ""
    if 'args' in params.keys():
        args = ' '.join(list(chain(*[[k,str(v)] for k,v in params["args"].items()])))
    if 'flags' in params.keys():
        flags = ' '.join(params["flags"])
    return args + " " + flags


def run_script(script_name, name, args, config, dependencies = "", node = ""):
    # Check if program is being run:
    command_prog = "command_" + script_name
    if command_prog in config:
        if not config[command_prog]:
            return False
    else:
        return False

    args.update(config)
    # Format Parameters
    for k,v in args.items():
        if k.startswith("command_") and type(v) != bool:
            args[k]["params"] = format_params(args[k])
    script = os.path.dirname(pipe.__file__) + "/scripts/" + script_name + ".sh"
    cleanup = os.path.dirname(pipe.__file__) + "/scripts/cleanup.sh"
    script = open(script, 'r').read()

    # Convert to DotMap for attribute access 
    args = edict(args)
    # Configure options:
    slurm_kwargs = {}
    if node:
        slurm_kwargs["node"] = node
    if dependencies:
        print dependencies
        dependencies = [x for x in dependencies if x is not None]
        print dependencies
        dependencies = ':'.join([str(x.job_id) for x in dependencies])
        if dependencies:
            slurm_kwargs["dependency"] = "afterok:" + dependencies
    else:
        dependencies = ""

    # Get headers:
    scriptf = script.format(**args)
    CHECK = [file_exists(x.split(":")[1].strip()) for x
             in scriptf.splitlines() if x.startswith("# CHECK:")]
    CLEANUP = [x.split(":")[1].strip() for x
             in scriptf.splitlines() if x.startswith("# CLEANUP:")]
    if not all(CHECK) or len(CHECK) == 0:
        with open("long_script.sh", 'a+') as f:
            f.write(scriptf)
        #print script.format(**args)
        job = Slurm(script_name + "-" + name, slurm_kwargs).run(scriptf)
        print "{script_name} - {name} ({job.job_id}) {dependencies}".format(**locals())
        for rm_file in CLEANUP:
            check_output(["sbatch","--dependency=afterany:" + str(job.job_id), cleanup, rm_file])
        return job
    else:
        return None



    
