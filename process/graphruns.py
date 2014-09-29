#!/usr/bin/env python
"""Graph Runs

Usage:
  graphruns.py [--input=INPUT] [--template=TEMPLATE] [--fileX=FILEX] [--fileZ=FILEZ] [--dir=DIR] [--clean] [--poly=<POLY>...] [--hide=<HIDE>...]

Options:
  --input=INPUT                 The name of an input file of converted run JSON
  --template=TEMPLATE           The name of the .gnuplot template file [default: template.gnuplot]
  --fileX=FILEX                 The output filename for the X logical graph [default: logicalX]
  --fileZ=FILEZ                 The output filename for the Z logical graph [default: logicalZ]
  --dir=DIR                     The name of the output directory [default: .]
  --poly=POLY                   A distance at which to calculate and draw a polynomial line
  --hide=HIDE                   A distance at which to hide from the graph
  --clean                       Deletes generated .dat and .gnuplot files.
"""

import os
import sys
import json
import math
import subprocess
from docopt import docopt
from string import Template

args = docopt(__doc__, version='Graph Runs v0.1')

if not args['--input']:
    rows = json.loads(sys.stdin.read())
else:
    if not os.path.isfile(args['--input']):
        sys.stdout.write("Count not read converted run JSON for graphing: %d\n", args['--input'])

    fp = open(args['--input'], 'r')
    rows = json.loads(fp.read())
    fp.close()

if args['--dir'] and not os.path.exists(args['--dir']): 
    os.makedirs(args['--dir'])

data_x = []
data_z = []
files = []

i = 0
for d in sorted(rows.keys(), key=lambda x: int(x)):
    stats = rows[d]
    d = int(d)
    first = True
    for stat in stats:
        #print stat
        (p, x, z) = map(float, (stat['p'], stat['pX'], stat['pZ']))
        
        #if x == 0 or z == 0:
        #    continue
        
        if first:    
            first = False
    
            filename = os.path.join(args['--dir'], "d%s.dat" % d)
            files.append(filename)
            out = open(filename, "w+")

            poly_min = 0.000001
            poly_max = 0.1

            dd = math.floor((d+1)/2)

            poly = pow(p, dd)
            A_x = x / poly
            A_z = z / poly

            poly = pow(poly_min, dd)
            x_1 = poly * A_x
            z_1 = poly * A_z
            
            poly = pow(poly_max, dd)
            x_2 = poly * A_x
            z_2 = poly * A_z

            filename = os.path.join(args['--dir'], "d%s_poly.dat" % d)
            files.append(filename)
            poly_out = open(filename, "w+")
            poly_out.write("%e %e %e\n" % (poly_min, x_1, z_1))
            poly_out.write("%e %e %e\n" % (poly_max, x_2, z_2))
            poly_out.close()
        
        out.write("%e %e %e\n" % (p, x, z))

    if not first:
        out.close()

    # Ensure that there is at least one data point
    if not first and str(d) not in args['--hide']:
        if str(d) in args['--poly']:
            data_x.append( "\"d%d_poly.dat\" using 1:2 notitle w l ls %d" % (d, i+19) )
            data_z.append( "\"d%d_poly.dat\" using 1:3 notitle w l ls %d" % (d, i+19) )
        
        data_x.append( "\"d%d.dat\" using 1:2 title \"d%d\" w lp ls %d" % (d, d, i+1) )
        data_z.append( "\"d%d.dat\" using 1:3 title \"d%d\" w lp ls %d" % (d, d, i+1) )
    
    i += 1

templ = open(args['--template'], "r")
t = Template(templ.read())
templ.close()

string_x = ", \\\n".join(data_x) + "\n"
filename = os.path.join(args['--dir'], "%s.gnuplot" % args['--fileX'])
files.append(filename)
l_x = open(filename, "w+")
l_x.write(t.substitute(filename=args['--fileX'], error_type='X', datastrings=string_x))
l_x.close()

string_z = ", \\\n".join(data_z) + "\n"
filename = os.path.join(args['--dir'], "%s.gnuplot" % args['--fileZ'])
files.append(filename)
l_z = open(filename, "w+")
l_z.write(t.substitute(filename=args['--fileZ'], error_type='Z', datastrings=string_z))
l_z.close()

os.chdir(args['--dir'])
subprocess.call(["gnuplot", "%s.gnuplot" % args['--fileX']])
subprocess.call(["gnuplot", "%s.gnuplot" % args['--fileZ']])

if args['--clean']:
    for f in files:
        os.remove(f)
