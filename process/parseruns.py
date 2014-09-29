#!/usr/bin/env python
"""Parse Runs

Usage:
  parseruns.py <dir> [--output=OUTPUT] [--quiet] [--incomplete | --complete] [--unbooted | --booted]

Options:
  dir               The runs directory
  output            The name of an output file to save JSON
  -i --incomplete   Incomplete jobs only    
  -c --complete     Complete jobs only
  -b --booted       Booted jobs only
  -u --unbooted     Unbooted jobs only
  -q --quiet        Does not warn when files cannot be read
"""

import os
import sys
import json
from docopt import docopt
from processruns.parser import *

args = docopt(__doc__, version='Parse Runs v0.1')

# Parse the runs and write to stdout
runs_dir = args['<dir>']

if not os.path.isdir(runs_dir):
    sys.stderr.write("Unable to find the given directory: %s\n" % runs_dir)
    sys.exit()

rows = parse_runs(runs_dir, args)
if not args['--output']:
    sys.stdout.write(json.dumps(rows))
else:
    fp = open(args['--output'], 'w+')
    fp.write(json.dumps(rows))
    fp.close()
