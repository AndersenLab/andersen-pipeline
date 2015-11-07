import pipe
import os
from itertools import chain
from collections import namedtuple
from easydict import EasyDict as edict

def format_params(params):
    print params
    args, flags = "", ""
    if 'args' in params["args"].keys():
        args = ' '.join(list(chain(*[[k,str(v)] for k,v in params["args"].items()])))
    if 'flags' in params.keys():
        flags = ' '.join(params["flags"])
    return args + " " + flags


def get_script(script_name, args, config):
    args.update(config)
    # Format Parameters
    for k,v in args.items():
        if k.startswith("command_"):
            args[k]["params"] = format_params(args[k])
    script = os.path.dirname(pipe.__file__) + "/scripts/" + script_name + ".sh"
    script = open(script, 'r').read()
    # Convert to DotMap for attribute access 
    args = edict(args)
    return script.format(**args)

    
