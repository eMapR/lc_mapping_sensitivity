#!/usr/bin/env python
"""

Run LandTrendr for the specified indices/dates (calls getGdriveFiles() from lt_gee_processing.py)

Usage:
    run_lt.py <param_file>
    run_lt.py -h | --help
    
    param_file - file ending in .py with all params specified

Options:
    -h --help      Show this screen.


"""


import os
import sys
import re
import docopt
import warnings
import pandas as pd
import ee
ee.Initialize()

import lt_gee_processing as ltee


def main(sleep=.5, njobs=10, silent=True):
    
    gDrive = ltee.authenticateGDrive()
    
    '''tasks = ltee.run_lt(featureCol,
                          featureKey,
                          featureValue,
                          aoiBuffer,
                          indexDict,
                          startYear,
                          endYear,
                          run_params,
                          startDay,
                          endDay,
                          mosaicType,
                          targetDay,
                          gDriveFolder,
                          outProj,
                          affine)#'''
    
    tasks = []
    if 'outDirPath' in globals():
        #ltee.listenAndDownladTasks(tasks, outDirPath, gDriveFolder, gDrive, outDir=outDirPath, clipFile=clipFile, silent=False, timestamp_offset=timestamp_offset)#'''
        ltee.listenAndDownladTasks(tasks, outDirPath, gDriveFolder, gDrive, silent=False, timestamp_offset=timestamp_offset, logFile=logFile)#'''



if __name__ == '__main__':
    try:
        cl_args = docopt.docopt(__doc__)
    except docopt.DocoptExit as e:
        import pdb; pdb.set_trace()
        print e.message
    
    # get rid of extra characters and 'help' entry from doc string 
    args = {re.sub('[<>-]*', '', k): v for k, v in cl_args.iteritems()
            if k not in ['--help','-h']} 

    sys.path.append(os.path.dirname(sys.argv[1]))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore') # raises SyntaxWarning because of *
        from run_lt_params import * # brings all variables into namespace
    
    sys.exit(main())
    
    