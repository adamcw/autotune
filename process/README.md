# Autotune Process Documentation #

## Overview ##

This is not an exhaustive guide to the tools in the pipeline. These scripts are
provided more as examples that can be used immediately with the Autotune
example file `ex1` to generate a graph of Logical Error Rate vs Physical Error
Rate. 

All tools in the pipeline provide usage and
help instructions when run with the flag `-h` or `--help`. These usage commands are highly
recommended for understanding options available outside of what is mentioned
here. All programs with `--input` or `--output` flags will read from `stdin` and 
`stdout`, respectively, by default.

The pipeline includes:

- run.py
- parseruns.py
- listruns.py 
- gatherruns.py
- convertruns
- graphruns.py

## Installation ##

The Autotune processing pipeline comprises of both code written in Python and
in C. You will need Python v2.7, and any version of `gcc` to compile the
C.

The processing pipeline should be automatically built during the `Autotune`
installation process. If this failed or you can install using:

    cd ./Autotune;
    make -C ./process install;

The pipeline depends on the Python packages `parse` and `docopt`, with
specifications available in `requirements.txt`. The installation of these
packages will depend on the environment in which they are being installed.

If you have `python` and `pip` installed, with the privileges to install new
packages, then the `make install` process will automatically install the
requirements for you. If you do not have permissions to install packages to the
global `python` installation, as in many cluster environments, you will need to
follow the virtualenv (virtual environment) process. `make install` will
indicate this if it fails.

If you succeeded at this step, you can skip down to "Example of Usage".

### Creating a Virtual Environment ###

A virtualenv is essentially just a bash script that will modify your PATH so
that commands like `python` and `pip` refer to locally installed versions
rather than globally installed ones. This allows the installation of packages
to the local `python` when permissions to install to the global `python` are
not available. 

If `pip` is not found, `make install` will inform you to use `make virtualenv` to
try and create a virtualenv (virtual environment) in which to install pip. This
can be done at any permission level provided `virtualenv` is provided to you.
If this works it will create a virtual environment containing `pip`, install
the required Python packages and then give you instructions on using the
virtualenv. 

If `virtualenv` is not found, `make virtualenv` will inform you to use `make
virtualenv_install` which will download `virtualenv` to the local directory,
install it, create a virtual environment containing `pip`, install the required
Python packages and then give you instructions on using the virtualenv. 

### Using a Virtual Environment ###

After successfully installing/creating a virtualenv, the build process should
show you the following instructions:

    *********************************
    PLEASE READ: USAGE INSTRUCTIONS 
    *********************************

    Type to update your PATH to use virtualenv (virtual environment):

        cd ./process;
        source AutotuneVE/bin/activate;

    To stop using this virtualenv.

        deactivate;

    Please see 'process/README.md' for more information.

In order to use a virtualenv you must "activate" it by running:

    cd ./process;
    source AutotuneVE/bin/activate;

This will update your path and from then on scripts looking for `python`, such
as those in the processing pipeline, will use the local version of `python` in
the virtualenv. 

To return to using the global version of `python` for any reason, the
virtualenv can be disabled using:

    deactivate;


## Example of Usage ##

    ./parseruns.py runs | ./gatherruns.py | ./convertruns | ./graphruns.py

## File Descriptions ##

### ex1 ###

`ex1` writes to a plain text file during its run time. The format specified in
this example file, `ex1`, had a number of lines known as headers, and then a
number of results lines, prefixed by ">>>".



### run.py ###

#### Usage ####

    Usage:
    run.py [--input=INPUT] [--dry] [--quiet]

    Options:
    --input=INPUT         A JSON file containing the run information to be executed
    --dry                 Performs a dry run. Shows what will be run, without actually running.
    -q --quiet            Does not warn when there are t_check issues 

#### Description ####

When using `run.py` to run `ex1` the output is placed into a known directory
structure which is then used by the following tool, `parseruns.py`.

These files are written into a directory structure as follows:

    {distance}/{probability}/{t_delete}/{run}

Distance is the distance of the code (-d), Probability is the overall physical
error rate (-p), t\_delete is the number of time steps that is kept in memory
during simulation (-t\_delete), and run is an increment counter of how many runs
of the given d/p/t\_delete combination have been performed.

This file takes input in JSON format which describes the options and paths to
execute. Multiple options can be provided in a list to create permutations of
the values provided allowing easy submission of parameter sweeps. An example
input is provided as `example_runs.json` in the directory.

#### Input Example ####

[
    { 
        "executable": "../ex/bin/ex1",
        "ems": "../ex/ems",
        "folder": "../../runs/{-d}/{-p}/{-t_delete}",
        "filename": "out_raw",
        "options": {
            "-d": [
                "3", "4", "5", "6"
            ],
            "-p": [
                "0.0001", "0.0002", "0.0005", 
                "0.001", "0.002", "0.005",
                "0.01", "0.02", "0.05"
            ],
            "-t_delete": "100",
            "-boot": "1",
            "-ems": "./ems/"
        }
    }
]

### parseruns.py ###

#### Usage ####

    Usage:
      parseruns.py <dir> [--output#OUTPUT] [--quiet] [--incomplete | --complete] [--unbooted | --booted]

    Options:
      dir               The runs directory
      output            The name of an output file to save JSON
      -i --incomplete   Incomplete jobs only    
      -c --complete     Complete jobs only
      -b --booted       Booted jobs only
      -u --unbooted     Unbooted jobs only
      -q --quiet        Does not warn when files cannot be read

#### Description ####

This script will parse the directory structure as generated through the use of
`run.py`. A different directory structure could be used by a modified `run.py`,
however this script will also need to be updated to reflect the change in
structure. This script parses all of the output files from the runs and outputs
them in a JSON format. This format contains two high level elements, structure
and rows. 

structure: A nested dictionary (object) in the order distance -> probability ->
t\_delete -> run. This is used for traversing the "directory" structure without
having to actually traverse any directories in anything that uses this JSON as
input.

rows: A dictionary (object) where each key is a row id, and the value of each
is a run. The id is {d}\_{p}\_{t\_delete}\_{run} which means a row can be looked up
easily if traversing the above mentioned structure by concatenating with
underscores.

Example of the row dictionary:

    'd',
    'p',
    't_delete',
    'run',
    'results': {
        'x_checks', 
        'z_checks', 
        'x_changes',
        'z_changes',
        't',
        'time'
    },
    'status': {
        'booted',
        'complete'
    },
    'options': {
        't_check',
        'big_t_max',
        'max_num_X',
        'max_num_Z',
        's0',
        's1'
    }



### listruns.py ###

#### Usage ####

    Usage:
      listruns.py [--input#INPUT] 

    Options:
      --input#INPUT            The name of an input file of parsed run JSON

#### Description ####

This is not part of the graphing pipeline but can be used to list all the runs
in a directory and see an overview of progress and completion. This program
takes the output of parseruns.py to operate. This script is useful while
performing a lot of runs to see at a glance the status of all the runs.  This
is an example of a piece of software that can be built on top of the output of
parseruns.py



### gatherruns.py ###

#### Usage ####

    Usage:
      gatherruns.py [--input#INPUT] [--output#OUTPUT] [--threshold#THRESHOLD] [--quiet]

    Options:
      --input#INPUT            The name of an input file of parsed run JSON
      --output#OUTPUT          The name of the output file to save the gathered results to
      --threshold#THRESHOLD    The number of X and Z changes required for a point to be included [default: 100]
      -q --quiet               Does not warn when multiple t_deletes are present

#### Description ####

The output from `parseruns.py` is raw, including all runs separately and all
versions of t\_delete. A graph can only use one value of t\_delete for each
value, so if multiple t\_delete values were run, only the highest will be
sampled by `gatherruns.py`. This file will also sum the results from multiple
runs to provide an overall number of changes, checks and time for each
combination. 

The output of this file is in plain text and is designed for use with `compile`. 



### convertruns ###

#### Usage ####

    Usage:
      convertruns [--input#INPUT] [--output#OUTPUT]

    Options:
      --input#INPUT         A plain-text input file
      --output#OUTPUT       A JSON-formatted output file

#### Description ####

This is compiled from `convertruns.c`, which is done simply `gcc convertruns.c
-o convertruns`. A Makefile has been provided as explained in the Installation.
This code will convert from the probability of discovering a logical error
every t\_check rounds of error correction, to the probability of discovering a
logical error per round. 

The output of this program is in JSON. A dictionary with distances as keys, and
lists of converted logical probabilities as values. 



### graphruns.py ###

#### Usage ####

    Usage:
      graphruns.py [--input#INPUT] [--template#TEMPLATE] [--fileX#FILEX] [--fileZ#FILEZ] [--dir#DIR] [--clean] [--poly#<POLY>...] [--hide#<HIDE>...]

    Options:
      --input#INPUT                 The name of an input file of converted run JSON
      --template#TEMPLATE           The name of the .gnuplot template file [default: template.gnuplot]
      --fileX#FILEX                 The output filename for the X logical graph [default: logicalX]
      --fileZ#FILEZ                 The output filename for the Z logical graph [default: logicalZ]
      --dir#DIR                     The name of the output directory [default: .]
      --poly#POLY                   A distance at which to calculate and draw a polynomial line
      --hide#HIDE                   A distance at which to hide from the graph
      --clean                       Deletes generated .dat and .gnuplot files.

#### Description ####

This is an example of a file which reads the JSON from `compile` and then
generates a graph of Logical Error Rate vs. Physical Error Rate for both X and
Z errors. This file makes use of `gnuplot` for plotting, however it could be
modified to use a plotting library or program of choice. 

In addition it provides the ability to generate data for polynomial lines that
can be specified as command line arguments.



## Manual Installation Process ##

### WARNING ###

The `make install` process should suffice for the large number of users, please
read the above instructions first. These instructions mostly detail what the
`make install` process already does automatically based on what it detects
available on your system. 

### With installation privileges ### 

If you can install new Python packages, installation is simple.

    pip install -r requirements.txt

This command will use pip to install all the required versions of the
libraries and you'll be good to go.

It should be noted, that to avoid package versioning conflicts between
different software, you may wish to use `virtualenv` detailed in the "Without
installation privileges" section anyway, instead of this direct method.

### Without installation privileges ###

Without installation privileges, there is still hope. A great piece of software
called `virtualenv` allows for the creation of a virtual user level Python
environment where packages can be installed. On environments where package
installation is generally not provided, you may have already come across this
in the past.

If `virtualenv` is already available in your environment, or can be installed
by a sysadmin, then you can simply type:

    virtualenv AutotuneVE;
    AutotuneVE/bin/pip install -r requirements.txt;

If not, you will need to download it:

    wget https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.10.tar.gz;
    tar xvfz virtualenv-1.10.tar.gz;
    python virtualenv-1.10/virtualenv.py AutotuneVE;
    AutotuneVE/bin/pip install -r requirements.txt;

This process with download and unpack `virtualenv` into the current directory.
It will then initialise a new virtual Python environment called `AutotuneVE`
and install the requirements into the environment. Do not worry if `pip` is not
available in your environment, `virtualenv` includes `pip`.

In order to use the pipline, you can "activate" the virtual environment, which
will essentially change the PATH to make use of the virtual environment. 

    source AutotuneVE/bin/activate;

Now any scripts calling `python`, such as these, will use the Python in the
virtual environment, and hence have access to the packages required.

To stop using the virtual environment you can simply type:

    deactivate;

If this process fails you can download the two libraries required manually,
however this is recommended only as a last resort:

    wget https://raw.github.com/r1chardj0n3s/parse/master/parse.py
    wget https://raw.github.com/docopt/docopt/master/docopt.py
