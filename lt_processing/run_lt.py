#!/usr/bin/env python
"""

Run LandTrendr for the specified indices/dates (calls getGdriveFiles() from lt_gee_processing.py)

Usage:
    run_LT.py <param_file>
    run_LT.py -h | --help

Options:
    -h --help      Show this screen.


"""


import os
import sys
import warnings
import pandas as pd
import ee
ee.Initialize()

import lt_gee_processing as ltgee


if __name__ == '__main__':
    
    sys.path.append(os.path.dirname(sys.argv[1]))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore') # raises SyntaxWarning because of *
        from run_lt_params import * # brings all variables into namespace
    sys.exit(ltgee.run_lt(featureCol,
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
                          affine))#'''
    