r"""
# send in job name and kwargs for slurm params:
>>> s = Slurm("job-name", {"account": "ucgd-kp", "partition": "ucgd-kp"})
>>> print(str(s))
#!/bin/sh
<BLANKLINE>
#SBATCH -e logs/job-name.%J.err
#SBATCH -o logs/job-name.%J.out
#SBATCH -J job-name
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=compute
#SBATCH --mem=32768
#SBATCH --mem-per-cpu=2730
<BLANKLINE>
set -eo pipefail -o nounset
<BLANKLINE>
{script}

>>> s.run("do stuff", _cmd="ls", name_addition="")

"""
from __future__ import print_function

import sys
import os
import subprocess
import tempfile
import atexit
import hashlib
import datetime
from itertools import cycle

TMPL = """\
#!/bin/bash

#SBATCH -e logs/{name}.%J.err
#SBATCH -o logs/{name}.%J.out
#SBATCH -J {name}

{header}

set -eo pipefail -o nounset

"""

def set_slurm_kwarg_default(slurm_kwargs, key, value):
    if key not in slurm_kwargs.keys():
        slurm_kwargs[key] = str(value)
    return slurm_kwargs

def tmp(suffix=".sh"):
    t = tempfile.mktemp(suffix=suffix)
    atexit.register(os.unlink, t)
    return t

class node_cycle(object):
    """
    Utility class for
    cycling through nodes
    and appropriately distributing jobs.
    """
    def __init__(self, nodelist):
        self.nodelist = cycle(nodelist)

    def get(self):
        return self.nodelist.next()

class job(object):
    def __init__(self, job_id, name, dependencies):
        self.job_id = job_id
        self.name = name
        self.datetime = datetime.datetime.now()
        self.dependencies = dependencies

    def __repr__(self):
        return "< " + str(self.job_id) + " : " + self.name + " >"

class Slurm(object):
    def __init__(self, name, slurm_kwargs={}, dependencies = [], tmpl=None, date_in_name=True, scripts_dir="scripts/"):
        if slurm_kwargs is None:
            slurm_kwargs = {}
        if tmpl is None:
            tmpl = TMPL

        header = []

        # Set slurm_kwargs defaults:
        slurm_kwargs = set_slurm_kwarg_default(slurm_kwargs, 'nodes', '1')
        slurm_kwargs = set_slurm_kwarg_default(slurm_kwargs, 'ntasks', '1')
        slurm_kwargs = set_slurm_kwarg_default(slurm_kwargs, 'cpus-per-task', '8')
        slurm_kwargs = set_slurm_kwarg_default(slurm_kwargs, 'mem', '32768')
        slurm_kwargs = set_slurm_kwarg_default(slurm_kwargs, 'mem-per-cpu', '4096')
        slurm_kwargs = set_slurm_kwarg_default(slurm_kwargs, 'nice', '10000')

        for k, v in slurm_kwargs.items():
            if len(k) > 1:
                k = "--" + k + "="
            else:
                k = "-" + k + " "
            header.append("#SBATCH %s%s" % (k, v))
        self.header = "\n".join(header)
        self.name = "".join(x for x in name.replace(" ", "-") if x.isalnum() or x == "-")
        self.tmpl = tmpl
        self.slurm_kwargs = slurm_kwargs
        if scripts_dir is not None:
            self.scripts_dir = os.path.abspath(scripts_dir)
        else:
            self.scripts_dir = None
        self.date_in_name = bool(date_in_name)


    def __str__(self):
        return self.tmpl.format(name=self.name, header=self.header)

    def _tmpfile(self):
        if self.scripts_dir is None:
            return tmp()
        else:
            if not os.path.exists(self.scripts_dir):
                os.makedirs(self.scripts_dir)
            return "%s/%s.sh" % (self.scripts_dir, self.name)

    def run(self, command, name_addition=None, cmd_kwargs=None, _cmd="sbatch", dependencies = []):
        """
        command: a bash command that you want to run
        name_addition: if not specified, the sha1 of the command to run
                       appended to job name. if it is "date", the yyyy-mm-dd
                       date will be added to the job name.
        cmd_kwargs: dict of extra arguments to fill in command
                   (so command itself can be a template).
        _cmd: submit command (change to "bash" for testing).
        """
        if name_addition is None:
            name_addition = hashlib.sha1(command.encode("utf-8")).hexdigest()[0:8]

        if self.date_in_name:
            name_addition += "-" + str(datetime.date.today())
        name_addition = name_addition.strip(" -")

        if cmd_kwargs is None:
            cmd_kwargs = {}

        n = self.name
        self.name = self.name.strip(" -")
        self.name += ("-" + name_addition.strip(" -"))

        tmpl = str(self)

        if "logs/" in tmpl and not os.path.exists("logs/"):
            os.makedirs("logs")

        with open(self._tmpfile(), "w") as sh:
            sh.write(tmpl.format(**cmd_kwargs) + command)

        res = subprocess.check_output([_cmd, sh.name]).strip()
        self.name = n
        if not res.startswith(b"Submitted batch"):
            return None
        job_id = int(res.split()[-1])
        return job(job_id, self.name , dependencies)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
