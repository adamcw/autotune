#!/usr/bin/env python
"""Gather Runs

Usage:
  gatherruns.py [--input=INPUT] [--output=OUTPUT] [--threshold=THRESHOLD] [--quiet]

Options:
  --input=INPUT            The name of an input file of parsed run JSON
  --output=OUTPUT          The name of the output file to save the gathered results to
  --threshold=THRESHOLD    The number of X and Z changes required for a point to be included [default: 100]
  -q --quiet               Does not warn when multiple t_deletes are present
"""

import sys
import json
from collections import defaultdict
from docopt import docopt

args = docopt(__doc__, version='Gather Runs v0.1')

if not args['--input']:
    rows = json.loads(sys.stdin.read())
else:
    fp = open(args['--input'], 'r')
    rows = json.loads(fp.read())
    fp.close()

output = defaultdict(lambda: defaultdict(tuple))
structure = rows['structure']
for d in structure.keys():
    for p in structure[d].keys():
        t_delete = str(max(map(int, structure[d][p].keys())))

        if len(structure[d][p].keys()) > 1 and not args['--quiet']:
            sys.stderr.write(
                "More than one value of t_delete detected for d = %s, p = %s. "
                "Defaulting to use the highest t_delete: %s\n" % (d, p, t_delete))

        t_check = 0
        for run in structure[d][p][t_delete]:
            row_id = "%s_%s_%s_%s" % (d, p, t_delete, run)
            o = rows['rows'][row_id]
            
            if not t_check:
                t_check = o['options']['t_check']

            if o['options']['t_check'] != t_check:
                continue

            row = (
                o['results']['x_checks'], 
                o['results']['z_checks'], 
                o['results']['x_changes'], 
                o['results']['z_changes'],
                o['results']['t'],
                o['results']['time'],
            )

            if not output[d][p]:
                output[d][p] = {
                    'd': d, 
                    'p': p, 
                    't_check': o['options']['t_check'], 
                    'res': (0,) * len(row)
                }

            output[d][p]['res'] = tuple(map(sum, zip(output[d][p]['res'], row)))

# Open file in writing mode to overwrite it if it exists
if args['--output']:
    fp = open(args['--output'], 'w')
    fp.close()

for d in map(str, sorted(map(int, output.keys()))):
    for p in sorted(output[d].keys(), key=lambda x: float(x)):
        row = output[d][p]

        # Only include points with more than THRESHOLD changes, otherwise the point
        # is likely to be erroneous and mess up any results
        threshold = int(args['--threshold'])
        if row['res'][2] < threshold or row['res'][3] < threshold:
            continue

        output_line = "%s %s %s %s %s %s %s %s %s\n" % ((d, p, row['t_check']) + row['res'])
        if args['--output']:
            fp = open(args['--output'], 'a')
            fp.write(output_line)
            fp.close()
        else:
            sys.stdout.write(output_line)
