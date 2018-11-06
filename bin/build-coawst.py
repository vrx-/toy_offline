#!/usr/bin/env python
# encoding: utf-8
"""
build_roms.py

Created by Rob Hetland on 2008-03-31.
Copyright (c) 2008 Texas A&M Univsersity. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import sys
import os
import getopt

help_message = '''
build-roms.py is a script used to build the Regional Ocean Modeling System

Usage:

    build-roms.py [options] ROMS_APPLICATION

The ROMS_APPLICATION, the sole argument to the script, must be either specified
on the command line, or present as an environmental variable.  If present, the
argument will override the environmental variable.

Options:

    -h, --help        Print this message

    --fort=[ifort]    Specify fortran compiler
    --netcdf=[4, 3]   Version of NetCDF to use [default is NetCDF4]
    -j=n              Parallel make [n = number of procs, default is 4]

    --mpi             Compile using mpi
    --clean           Delete Build scratch directory before build
    --debug           Compile in debug mode
    --root            Specify ROMS_ROOT_DIR

Default values can be taken from environmental variables, if present.  Command
line options always override defaults.  Environmental variables are:

    FORT            for fortran compiler, overriden by --fort
    ROMS_ROOT_DIR   for ROMS root directory, overriden by --root
'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def main(argv=None):
    if argv is None:
        argv = sys.argv

    opts, args = getopt.getopt(argv[1:], "hj", ["help", "mpi", "fort=", "netcdf=", "debug", "clean", "root"])

    # Get default values from environment, if present

    if 'ROMS_ROOT_DIR' in os.environ:
        ROMS_ROOT_DIR = os.environ['ROMS_ROOT_DIR']
    else:
        ROMS_ROOT_DIR = os.path.join(os.environ['HOME'], 'COSMIO/COAWST')

    USE_MPI = ''
    USE_MPIF90 = ''

    USE_DEBUG = ''
    USE_LARGE = 'on'

    if 'FORT' in os.environ:
        FORT = os.environ['FORT']
    else:
        FORT = 'gfortran'

    USE_NETCDF4 = ''

    NP = 4

    # option processing
    for option, value in opts:
        if option in ("-h", "--help"):
            raise Usage(help_message)
        if option == "--mpi":
            USE_MPI = 'on'
            USE_MPIF90 = 'on'
        if option == "--fort":
            FORT = value
        if option == "--netcdf":
            if value != '4':
                USE_NETCDF4 = ''
        if option == "--root":
            ROMS_ROOT_DIR = value
        if option == "--debug":
            USE_DEBUG = 'on'
        if option == "-j":
            NP = value
        if option == "--clean":
            print('Removing Build directory...')
            os.system('rm -rf Build')

    if 'ROMS_APPLICATION' in os.environ:
        ROMS_APPLICATION = os.environ['ROMS_APPLICATION'].upper()
    else:
        assert len(args) == 1, 'ERROR: Must specify ROMS_APPLICATION'

    if len(args) == 1:
        # if 'upwelling' in args[0]:
        #     ROMS_APPLICATION = 'UPWELLING'
        # else:
        ROMS_APPLICATION = args[0].upper()

    headerfile = args[0] + '.h'
    assert os.path.exists(headerfile), \
        'ERROR: header file %s does not exist' % headerfile

    MY_PROJECT_DIR = os.getcwd()
    MY_HEADER_DIR = MY_PROJECT_DIR
    MY_ANALYTICAL_DIR = os.path.join(MY_PROJECT_DIR, 'Functionals')
    BINDIR = MY_PROJECT_DIR
    SCRATCH_DIR = os.path.join(MY_PROJECT_DIR, 'Build')

    ENV_STR = 'ROMS_APPLICATION=%s FORT=%s USE_MPI=%s USE_MPIF90=%s USE_NETCDF4=%s USE_DEBUG=%s USE_LARGE=%s MY_HEADER_DIR=%s MY_ANALYTICAL_DIR=%s BINDIR=%s SCRATCH_DIR=%s' % (ROMS_APPLICATION, FORT, USE_MPI, USE_MPIF90, USE_NETCDF4, USE_DEBUG, USE_LARGE, MY_HEADER_DIR, MY_ANALYTICAL_DIR, BINDIR, SCRATCH_DIR)

    print(ENV_STR.replace(' ', '\n'))

    os.system("cd %s ; %s make -j %s" % (ROMS_ROOT_DIR, ENV_STR, NP))


if __name__ == "__main__":
    sys.exit(main())

