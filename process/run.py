#!/usr/bin/env python
"""Run

Usage:
  run.py [--input=INPUT] [--dry] [--pbs] [--quiet]

Options:
  --input=INPUT         A JSON file containing the run information to be executed
  --dry                 Performs a dry run. Shows what will be run, without actually running.
  --pbs                 Create and submit a PBS script instead
  -q --quiet            Does not warn when there are t_check issues 
"""

import os
import sys
import json
import shutil
import itertools
from processruns.runner import permute_runs, run_command, pbs_command
from processruns.parser import parse_headers
from docopt import docopt

args = docopt(__doc__, version='Run')

if not args['--input']:
    fp = sys.stdin
else:
    if not os.path.isfile(args['--input']):
        sys.stderr.write("Could not find input file: %s\n" % args['--input'])
        sys.exit()
    fp = open(args['--input'], 'r')

try:
    RUNS = json.load(fp)
except ValueError:
    sys.stderr.write("Input JSON could not be read.\n")
    sys.exit()

if fp != sys.stdin:
    fp.close()

runs = permute_runs(RUNS)
for r in runs:
    if not os.path.isfile(r['executable']):
        sys.stderr.write("Could not find executable: %s\n" % r['executable'])
        sys.exit()

    if not os.path.isdir(r['ems']):
        sys.stderr.write("Could not find error model directory: %s\n" % r['ems'])
        sys.exit()

    # Create the folder recursively, naming each run accordingly
    for i in range(1, 1000):
        folder = os.path.join(r['folder'], 'run%d' % i)
        if not os.path.isdir(folder):
            if not args['--dry']:
                os.makedirs(folder)
            break
  
    # It is currently unsupported to be able to deal with data where two
    # identical runs are performed with differing t_check values. We therefore
    # will detect these events and either error, or gracefully correct where
    # possible. 
    if i > 1:
        filename = os.path.join(r['folder'], 'run1', 'out_raw')
        if not os.path.isfile(filename):
            if not args['--quiet']:
                sys.stderr.write("Could not find output of previous run: %s\n" % filename)
            continue

        try:
            headers = parse_headers(filename)
        except FileParserError as e:
            sys.stderr.write(e)
            continue

        # Did the header have to boot?
        h_boot = True if 'boot' in headers[0]['opts'] else False

        # Does this run need to boot?
        if '-boot' in r['options'] and r['options']['-boot'] == '1':
            r_boot = True
        else:
            r_boot = False

        # If the first run booted, then do not boot subsequent runs
        if h_boot and r_boot:
            if not args['--quiet']:
                sys.stderr.write("Warning: Already have a run with t_check, disabling boot option and continuing.\n")
            r['options']['-boot'] = '0'

        # Check that there are two booted headers
        if h_boot:
            if len(headers) != 2 or 'new t_check' not in headers[-1]['opts']:
                if not args['--quiet']:
                    sys.stderr.write("Error: Original boot process has not completed. %s\n" % filename)
                continue

            t_check = headers[1]['opts']['new t_check']
        else:
            if 't_check' not in headers[0]['opts']:
                if not args['--quiet']:
                    sys.stderr.write("Error: Unable to parse t_check from original run: %s\n" % filename)
                continue

            t_check = headers[0]['opts']['t_check']

        # If the t_check doesn't exist, or doesn't match, set it
        if '-t_check' not in r['options'] or r['options']['-t_check'] != t_check:
            r['options']['-t_check'] = t_check

    # Copy the executable to the new folder
    exe_name = os.path.join(".", os.path.basename(r['executable']))
    exe_src = r['executable']
    exe_dest = os.path.join(folder, exe_name)
    if not args['--dry']:
        shutil.copy2(exe_src, exe_dest)
    
    # Copy the error models
    ems_name = 'ems'
    ems_src = r['ems'] 
    ems_dest = os.path.join(folder, ems_name)
    if not args['--dry']:
        shutil.copytree(ems_src, ems_dest) 
    
    # Output Destination
    out_dest = os.path.join(folder, r['filename'])

    # Flatten options from dict to list to pass to the command
    opts = list(itertools.chain(*r['options'].items()))

    # Append the executable to the command
    command = [exe_name] + opts
    if args['--dry']:
        pid = 0
    elif args['--pbs']:
        pbs_file = 'submit.pbs'
        job_name = "%s_%s" % (r['options']['-d'], r['options']['-p'])
        
        if 'queue' not in r or 'wallhours' not in r:
            sys.stderr.write("You must include `wallhours` and `queue` in your `--input` when using PBS submit.\n");
            continue

        pid = pbs_command(command, r['filename'], pbs_file, job_name=job_name,
                folder=folder, cwd=folder, queue=r['queue'],
                hours=r['wallhours'])
    else:
        pid = run_command(command, r['filename'], folder=folder, cwd=folder)

    print("Started: %s (%s): %d" % (exe_name, out_dest, pid))
