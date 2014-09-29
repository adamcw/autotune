import os
import sys
import parse
import subprocess
from collections import defaultdict

# Parser for the data lines of the output file
DEFAULT_PARSER = parse.compile(
    ">>> {time:g} {} {t:d} {} Last: {lastx:d} {lastz:d} {} Checks: {xchecks:d} {zchecks:d} {} Changes: {xchanges:d} {zchanges:d}")

# Parser for the option flags of the output file
DEFAULT_OPT_PARSER = parse.compile(
    "{opt} = {val}")

# The maximum number of header blocks (eg. before boot header and after boot header)
MAX_NUM_HEADERS = 2

# The name of the output file from the Autotune code
AUTOTUNE_OUTPUT = 'out_raw'

class FileParserError(Exception):
    pass

def parse_headers(filename, parser=None, opt_parser=None):
    '''Parses the header lines from the out_raw output of our example files.
    Puts lines of the form "X = Y" into a dict called 'opts', all other lines
    are placed in raw in 'misc'.'''

    if not parser:
        parser = DEFAULT_PARSER
    if not isinstance(parser, parse.Parser):
        raise ValueError("The provided parser is not of type parse.Parser\n")

    if not opt_parser:
        opt_parser = DEFAULT_OPT_PARSER
    if not isinstance(opt_parser, parse.Parser):
        raise ValueError("The provided option parser is not of type parse.Parser\n")

    fp = open(filename, 'r+')

    headers = []
    in_header = False
    for line in fp:
        line = line.strip()
        if parser.parse(line):
            # If we already have two headers, stop searching
            if len(headers) == MAX_NUM_HEADERS:
                break
            if headers[-1]['opts'].get('new t_check', False):
                break
            in_header = False
        else:
            # We have entered a new header block
            if not in_header:
                in_header = True
                headers.append({ 
                    'opts': {},
                    'misc': [],
                })

            # See if we have an option, add to the header appropriately
            opt = opt_parser.parse(line)
            if opt:
                headers[-1]['opts'][opt['opt']] = opt['val']
            else:
                headers[-1]['misc'].append(line)

    fp.close()
     
    return headers

def parse_tail(filename, parser=None):
    '''Reads the last line of a file and parses it with the given parser.
    Returns the parsed line.'''

    if not parser:
        parser = DEFAULT_PARSER
    if not isinstance(parser, parse.Parser):
        raise ValueError("The provided parser is not of type parse.Parser\n")

    out = subprocess.check_output(['tail', '-n', '1', filename]).strip()
    pa = parser.parse(out)
    
    if not pa:
        out = subprocess.check_output(['tail', '-n', '2', filename]).strip().split("\n")[0]
        pa = parser.parse(out)
        if not pa:
            raise FileParserError("An error occured reading the file: %s\n" % filename)
        
    return pa

def get_dirs(*args):
    return [x for x in os.listdir(os.path.join(*args)) 
        if os.path.isdir(os.path.join(*args + (x,)))]

def parse_runs(runs_dir, args):
    '''Parses a runs directory of the heirarchy: d -> p -> t_delete -> run.
    Descends the directory structure and parses the output of each run,
    returning a dictionary of rows indexed by d_p_{t_delete}_run.'''

    structure = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    rows = dict()

    try:
        distances = get_dirs(runs_dir)
    except OSError:
        sys.stderr.write("Could not find runs directory: %s\n" % path)
        sys.exit()

    for d in distances:
        probabilities = get_dirs(runs_dir, d)
        for p in probabilities:
            t_deletes = get_dirs(runs_dir, d, p)
            for t_delete in t_deletes:
                runs = get_dirs(runs_dir, d, p, t_delete)
                for run in runs:
                    #if not (int(run[3:]) > 10 and int(run[3:]) < 15):
                    #    continue

                    filename = os.path.join(runs_dir, d, p, t_delete, run, AUTOTUNE_OUTPUT)
                    if not os.path.isfile(filename):
                        if not args['--quiet']:
                            sys.stderr.write("Unable to read file: %s\n" % filename)
                        continue

                    try:
                        headers = parse_headers(filename)
                        pa = parse_tail(filename)
                    except FileParserError as e:
                        if not args['--quiet']:
                            sys.stderr.write(e.message)
                        continue
                        
                    opts = headers[-1]['opts']

                    # Grab the value of t_check
                    t_check = opts.get('new t_check', opts.get('t_check', False))

                    # Determine if the run has booted or not
                    boot = headers[0]['opts'].get('boot', False)
                    if not boot:
                        booted = "NA"
                    elif len(headers) == 2:
                        booted = "Yes"
                    else:
                        booted = "No"

                    # Determine if the run is completed
                    if int(pa['xchanges']) >= int(opts['max_num_X']) and int(pa['zchanges']) >= int(opts['max_num_Z']):
                        done = "Yes"
                    else:
                        done = "No"

                    # Filter results
                    filters = [
                        args['--complete'] and done == "No",
                        args['--incomplete'] and done == "Yes",
                        args['--booted'] and booted == "No",
                        args['--unbooted'] and booted == "Yes",
                    ]
                    if any(filters):
                        continue

                    structure[d][p][t_delete].append(run)
                    
                    row_id = '%s_%s_%s_%s' % (d, p, t_delete, run)
                    rows[row_id] = {
                        'd': int(d),
                        'p': float(p),
                        't_delete': int(t_delete),
                        'run': run,
                        'results': {
                            'x_checks': pa['xchecks'], 
                            'z_checks': pa['zchecks'], 
                            'x_changes': pa['xchanges'], 
                            'z_changes': pa['zchanges'],
                            't': pa['t'],
                            'time': pa['time'],
                        },
                        'status': {
                            'booted': booted,
                            'complete': done,
                        },
                        'options': {
                            't_check': int(t_check),
                            'big_t_max': long(opts['big_t_max']),
                            'max_num_X': int(opts['max_num_X']),
                            'max_num_Z': int(opts['max_num_Z']),
                            's0': int(opts['s0'][:-1]),
                            's1': int(opts['s1'][:-1]),
                        }
                    }
    
    return {'structure': structure, 'rows': rows}
