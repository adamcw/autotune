import os
import itertools
import subprocess

def run_command(command, filename, folder="./", cwd="./"):
    if not os.path.isdir(folder):
        os.makedirs(folder)

    fp = open(os.path.join(folder, filename), "w+")
    command = [x for x in command if x]
    return subprocess.Popen(["nohup"] + command, stdin=fp, stdout=fp, stderr=fp, cwd=cwd).pid

def flatten(l, ltypes=(list, tuple)):
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)

def permute_runs(schema):
    runs = list()
    for run in schema:
        opts = []
        for key, val in run['options'].iteritems():
            if not isinstance(val, list):
                val = [val]

            if not opts:
                opts = list(val)
            else:
                opts = list(itertools.product(opts, val))

        for a in opts:
            b = dict(zip(run['options'].keys(), flatten(a)))
            runs.append({
                'executable': run['executable'].format(**b),
                'ems': run['ems'].format(**b),
                'filename': run['filename'].format(**b),
                'folder': run['folder'].format(**b),
                'options': b,
            })
    return runs
