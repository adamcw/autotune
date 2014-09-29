#!/usr/bin/env python
"""List Runs

Usage:
  listruns.py [--input=INPUT] 

Options:
  --input=INPUT            The name of an input file of parsed run JSON
"""

import sys
import json
from docopt import docopt

args = docopt(__doc__, version='List Runs v0.1')

# The output listing format
OUTPUT_FORMAT = "{d:<3} {p:<7} {t_delete:<5} {time:<12} {t:<12} {t_check:<7} {xchecks:<9} {zchecks:<9} {xchanges:<9} {zchanges:<9} {booted:<7} {done:<5}"

print(OUTPUT_FORMAT.format(
    d="d",
    p="p",
    t_delete="t_del",
    run="run",
    time="time",
    t="t",
    t_check="t_chk",
    xchecks='X Checks', 
    zchecks='Z Checks', 
    xchanges='X Changes', 
    zchanges='Z Changes',
    booted='Booted',
    done='Done',
))

if not args['--input']:
    rows = json.loads(sys.stdin.read())
else:
    fp = open(args['--input'], 'r')
    rows = json.loads(fp.read())
    fp.close()

for o in rows['rows'].values():
    print(OUTPUT_FORMAT.format(
        d=o['d'],
        p=o['p'],
        t_delete=o['t_delete'],
        run=o['run'],
        time=o['results']['time'],
        t=o['results']['t'],
        t_check=o['options']['t_check'],
        xchecks=o['results']['x_checks'], 
        zchecks=o['results']['z_checks'], 
        xchanges=o['results']['x_changes'], 
        zchanges=o['results']['z_changes'],
        booted=o['status']['booted'],
        done=o['status']['complete'],
    ))
