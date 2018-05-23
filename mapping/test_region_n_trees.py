# -*- coding: utf-8 -*-
"""
Created on Wed May 23 08:41:53 2018

@author: shooper
"""
import os, sys, time
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob
from multiprocessing import Pool


sys.path.append('/vol/v2/stem/stem-git/scripts')
from randomforest import test_n_trees


def write_params(train_txt):
        
    param_txt = os.path.join(os.path.dirname(train_txt), 'test_n_trees_params.txt')
    target_col = '_'.join(os.path.basename(train_txt).split('_')[:2])
    
    with open(param_txt, 'w') as txt:
        txt.write('Paramater File for test_n_trees.py\n\n')
        txt.write('target_col; %s\n' % target_col)
        txt.write('sample_txt; %s\n' % train_txt)
        txt.write('var_txt; /vol/v2/stem/stem-git/param_files/var_info_conus.txt\n')
        txt.write('max_trees; 500\n')
        txt.write('step; 25')
    
    return param_txt


def test(args):
    
    i, n_tests, t0, pattern, param_txt = args
    
    
    results = test_n_trees.main(param_txt, silent=True)
    
    # Configure the terminal message
    max_bar_size = 50
    lapsed_time = (time.time() - t0)/60
    percent_done = float(i)/n_tests
    n_chars = int(percent_done * max_bar_size)
    n_blank = max_bar_size - n_chars
    progress_bar = u'\u2588' * n_chars + ' '  * n_blank + '|'
    msg = pattern % (i + 1, n_tests, percent_done * 100, progress_bar, lapsed_time)
    n_returns = len(re.findall('\n', msg))
    prefix = '\033[A' * n_returns
    
    sys.stdout.write(prefix + msg)
    sys.stdout.flush()
    
    return results


def main(out_png):
    
    train_samples = glob('/vol/v1/proj/stem_improv_paper/region_models/region_*/samples/nlcd*/*predictors.txt') + \
                    glob('/vol/v1/proj/stem_improv_paper/region_models/region_*/samples/canopy*/*predictors.txt') + \
                    glob('/vol/v1/proj/stem_improv_paper/region_models/region_*/samples/imperv*/*predictors.txt')
    
    n_tests = len(train_samples)
    t0 = time.time()
    pattern = '\nFinished testing %s of %s (%.1f%%)\n%s\n%.1f minutes\n'
    print '\nTesting number of trees for %s samples...%s\n' % (n_tests, '\n' * len(re.findall('\n', pattern)))
    
    args = [[i, n_tests, t0, pattern, write_params(train_txt)] for i, train_txt in enumerate(train_samples)]
    pool = Pool(min(n_tests, 40))
    results = pool.map(test, args, 1)
    pool.close()
    pool.join()
    try:
        
        for n_trees, errors in results:
            plt.plot(n_trees, errors, '-', alpha=.2)
        
        plt.savefig(out_png, dpi=300)
    except:
        import pdb; pdb.set_trace()


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))