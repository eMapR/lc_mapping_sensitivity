#!/usr/bin/env python
"""
Download files from a specified Google Drive folder (calls getGdriveFiles() from lt_gee_processing.py)

Usage:
    clip_and_decompose.py <gDirName> <outDirPath> [--njobs=<int>]
    clip_and_decompose.py -h | --help


Options:
    -h --help      Show this screen.
    --njobs=<int>  Integer number of cores to use[default: 1]

"""


import sys
import re
import lt_gee_processing as ltee
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
    
    sys.exit(ltee.getGdriveFiles(**args))