#!/usr/bin/env python
"""
Utility function to call clipAndDecompose() from lt_gee_processing.py.

Usage:
    clip_and_decompose.py <chunkDir> <outDir> <clipFile> [--njobs=<int>] [--tileIdField=<str>] [--proj=<str>] [--returnOutDirs=<bool>]
    clip_and_decompose.py -h | --help


Options:
    -h --help                Show this screen.
    --njobs=<int>            Integer number of cores to use[default: 1]
    --tileIdField=<str>        Field in clipFile containing tile IDs. [default: name]
    --proj=<str>             OGR-interpretable projection string [default: EPSG:5070]
    --returnOutDirs=<bool>   Boolean indicating whether to return lowest level output subdirectories [default:'']
"""


import sys
import re
import lt_gee_processing as lt
import docopt


if __name__ == '__main__':
    try:
        cl_args = docopt.docopt(__doc__)
    except docopt.DocoptExit as e:
        import pdb; pdb.set_trace()
        print e.message
    
    # get rid of extra characters from doc string and 'help' entry
    args = {re.sub('[<>-]*', '', k): v for k, v in cl_args.iteritems()
            if k != '--help' and k != '-h'}  
    
    try:
        args['njobs'] = int(args['njobs'])
    except:
        raise ValueError('njobs not understood: ' + args['njobs'])

    if args['returnOutDirs']:        
        if args['returnOutDirs'].lower() == 'false':
            args['returnOutDirs'] = False
    
    sys.exit(lt.clipAndDecompose(**args))