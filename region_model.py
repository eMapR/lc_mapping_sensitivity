# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 10:47:05 2017

@author: shooper
"""

import os, sys, time, shutil
import pandas as pd

import wall_to_wall_confusion_matrix as confusion
import plot_importance

sys.path.append('/vol/v2/stem/stem-git/scripts')
import stem
import get_stratified_random_pixels as gsrp
import extract_xy_by_mosaic as extract
from randomforest import train_rf
from randomforest import predict_rf



def write_params(out_txt, params, param_names, opt_params=[]):
    
    with open(out_txt, 'w') as txt:
        for p in param_names:
            txt.write(p + '; ' + params[p] + '\n')
        
        if len(opt_params) > 0:
            txt.write('\nOptional parameters\n')
            for p in opt_params:
                txt.write(p + '; ' + params[p] + '\n')
    

def write_extract_params(out_txt, params):
    
    # Read in, then write out, var_txt
    var_txt = params['extract_var_txt']
    pd.read_csv(var_txt, sep='\t', index_col='var_name').to_csv(out_txt, sep='\t')
    
    extract_params = ['mosaic_path', 'xy_txt', 'years', 'target_col', 'tile_id_field']
    with open(out_txt, 'a') as txt:
        txt.write('\n')
        for p in extract_params:
            txt.write(p + '; ' + params[p] + '\n')


def main(param_txt):
    
    t0 = time.time()
    inputs = gsrp.read_params(param_txt)
    for var in inputs:
        exec ("{0} = str({1})").format(var, inputs[var])
    
    params = {k:v.replace('"','') for k, v in inputs.iteritems()}
    
    # If xy_txt isn't in params, there's no sample yet for extract_xy_by_mosaic
    if 'xy_txt' not in params and 'sample_txt' not in params:
        xy_param_txt = os.path.join(scratch_dir, 'get_stratified_random_pixels_params.txt')
        params['n_tiles'] = n_tiles_sample
        sample_params = ['raster_path', 'col_name', 'data_band', 'nodata', 'n_sample', 'bins', 'out_dir', 'pct_train', 'data_type', 'n_tiles', 'sampling_scheme', 'min_sample', 'max_sample']

        if 'bin_scale' in inputs:
            sample_params.append('bin_scale')
        write_params(xy_param_txt, params, sample_params)
    
        # Make random sample
        xy_txt = gsrp.main(xy_param_txt)
        params['xy_txt'] = xy_txt
        
    # If sample_txt isn't in params, get predictor values
    if 'sample_txt' not in params:
        extract_param_txt = os.path.join(scratch_dir, 'extract_xy_by_mosaic_params.txt')
        write_extract_params(extract_param_txt, params)
        sample_txt = extract.main(extract_param_txt)
        params['sample_txt'] = sample_txt
        
    
    # train
    train_param_txt = os.path.join(scratch_dir, 'train_rf_params.txt')
    train_params = ['target_col', 'sample_txt', 'var_txt', 'out_dir']
    optional_params = ['n_trees', 'n_jobs', 'max_features', 'model_type']
    write_params(train_param_txt, params, train_params, optional_params)
    rf_path = train_rf.main(train_param_txt)
    model_dir = os.path.dirname(rf_path)
    var_info_txt = os.path.join(model_dir, os.path.basename(params['var_txt']))
    plot_importance.main(var_info_txt)
    
    # predict
    train_params = os.path.join(model_dir, 'train_rf_params.txt')
    params['train_params'] = train_params
    params['rf_path'] = rf_path
    params['out_dir'] = model_dir
    predict_param_txt = os.path.join(scratch_dir, 'predict_rf_params.txt')
    
    shutil.copy2(param_txt, model_dir)
    predict_params = ['train_params', 'var_txt', 'mask_path', 'nodata', 'rf_path', 'out_dir']
    optional_params = ['n_tiles']
    params['n_tiles'] = n_tiles_predict
    write_params(predict_param_txt, params, predict_params, optional_params)
    pred_path = predict_rf.main(predict_param_txt)
    
    #evalutate
    eval_dir = os.path.join(model_dir, 'evaluation')
    if not os.path.isdir(eval_dir):
        os.mkdir(eval_dir)
    out_cm_txt = os.path.join(eval_dir, 'confusion.txt')
    confusion.main(raster_path, pred_path, out_cm_txt, target_col, int(nodata), int(nodata))
    
    
    print '\nTotal time: %.1f minutes' % ((time.time() - t0)/60)
    

if __name__ == '__main__':
    
    sys.exit(main(*sys.argv[1:]))
    