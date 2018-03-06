# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 17:36:22 2018

@author: shooper
"""
import os
import sys
sys.path.append('/vol/v1/general_files/script_library/earth_engine/')
from make_history import main as make_history
import lt_gee_processing as LT




def main():
    
    
    # find all dirs that have a temp dir but no perm dir to match
    # if no perm dir, call decompose
    # call make history
    clip_file = '/vol/v1/proj/stem_improv_paper/vector/regions/study_region_tiles.shp'
    njobs = 10
    search_dir = '/vol/v1/proj/stem_improv_paper/raster'
    out_dir = '/vol/v3/LTV3'
    tempdirs = [os.path.join(search_dir, d) for d in os.listdir(search_dir) if 'temp' in d]
    permdirs = [os.path.join(search_dir, d) for d in os.listdir(search_dir) if not 'temp' in d]
    all_dirs = tempdirs + permdirs
    complete = []
    
    for i, directory in enumerate(tempdirs):
        print 'Processing files in %s of %s directories' % (i, len(all_dirs))
        permdir = directory.replace('_temp', '')
        if permdir in complete:
            continue
        if not permdir in permdirs:
            print 'Clipping and decomposing stacks...\n'
            output_dirs = LT.clipAndDecompose(directory, out_dir, clip_file, njobs=njobs, returnOutDirs=True)
        
        print ''
        for directory in output_dirs:
            make_history(directory, search_str='*ftv*.tif', njobs=njobs)
        #complete.append(permdir)
    
    for d in complete:
        print d
    
    
    
    

if __name__ == '__main__':
    
    '''sys.path.append(os.path.dirname(sys.argv[1]))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore') # raises SyntaxWarning because of *
        from ee_test_params import * # brings all variables into namespace'''
    sys.exit(main())#'''