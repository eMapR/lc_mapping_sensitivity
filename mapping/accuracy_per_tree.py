# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 17:27:20 2017

@author: shooper
"""

import os, sys, time, random
import cPickle as pickle
import pandas as pd
import numpy as np
from osgeo import gdal

sys.path.append('/vol/v2/stem/stem-git/scripts')
import stem
from evaluation.evaluation import confusion_matrix_by_area
from get_stratified_random_pixels import parse_bins, read_params
from extract_xy_by_mosaic import extract_var


def get_predictors(tile_ar, var_info, tx, nodata, coords, temp_nodata=-9999):
    
    tile_mask = tile_ar == 0
    tile_ar[tile_mask] = nodata
    # Get the ids of tiles this kernel covers
    tile_ids = np.unique(tile_ar)
    #tile_strs = ['0' + str(tile) for tile in tile_ids if tile!=nodata]
    tile_strs = [str(tile) for tile in tile_ids if tile!=nodata]

    # Get an array of predictors where each column is a flattened 2D array of a
    #   single predictor variable
    predictors = stem.get_predictors(var_info, tx, tile_strs, tile_ar, coords, tile_mask, temp_nodata, 1)
    predictor_mask = ~ np.any(predictors==temp_nodata, axis=1)
    
    return predictors, predictor_mask


def replace_val(array, replace_dict):
    
    max_val = np.max(array)
        
    map_array = np.arange(0, max_val + 1)
    keys = [k for k in sorted(replace_dict.keys()) if k in map_array]
    map_array[keys] = [replace_dict[k] for k in keys]
    ar = map_array[array]

    return ar


'''def tree_accuracy(dt, predictors, bins, classes, ar_ref, ar_tile, predictor_mask, mask, r_nodata, target_col, temp_nodata=-9999, dtype=np.int16):
    

    predictions = dt.predict(predictors[~predictor_mask]).astype(dtype)
    ar_pred = np.full(predictors.shape[0], r_nodata, dtype=dtype)
    ar_pred[~predictor_mask] = predictions
    class_dict = dict(zip(dt.classes_, classes))
    class_dict[r_nodata] = r_nodata
    ar_pred = replace_val(ar_pred.reshape(ar_ref.shape), class_dict)
    
    sample = pd.DataFrame({target_col: ar_ref[~mask], 'prediction': ar_pred[~mask]})
    cm = confusion_matrix_by_area(ar_pred, ar_ref, sample, r_nodata, r_nodata, mask, bins=bins, target_col=target_col)
    
    labels = ['%s_%s' % (l, u) for l, u in bins]
    p_acc = cm.ix['producer', labels]
    u_acc = cm.ix[labels, 'user']
    p_acc.set_axis(0, ['p_%s' % i.split('_')[1] for i in p_acc.index])
    u_acc.set_axis(0, ['u_%s' % i.split('_')[1] for i in u_acc.index])
    s = pd.concat([p_acc, u_acc])
    s['accuracy'] =  cm.ix['producer','user']
    s['kappa'] = cm.ix['producer', 'kappa']

    return s'''
    

def get_test_sample(sample, ar, bins, n_per_bin, nodata, target_col, tx, ar_tile, var_info, test_sample=None):
    
    # If test_sample not given, make one
    if not type(test_sample) == pd.core.frame.DataFrame:
        # Set all training pixels to nodata
        ar[sample.row, sample.col] = nodata
        nodata_mask = ar != nodata
        rows, cols = np.indices(ar.shape)
        test_rows = []
        test_cols = []
        # Sample for each bin
        for this_min, this_max in bins:
            mask = (ar > this_min) & (ar <= this_max) & nodata_mask
            masked_rows = rows[mask]
            masked_cols = cols[mask]
            try:
                these_rows = random.sample(masked_rows, n_per_bin)
                these_cols = random.sample(masked_cols, n_per_bin)
            except:
                print (('Not enough pixels between {0} and {1} to generate {2}' +\
                    ' random samples. Returning all {3} pixels for this bin.'))\
                    .format(this_min, this_max, n_per_bin, masked_rows.size)
                these_rows = masked_rows
                these_cols = masked_cols
            test_rows.extend(these_rows)
            test_cols.extend(these_cols)
        
        ul_x, x_res, x_rot, ul_y, y_rot, y_res = tx
        test_sample = pd.DataFrame(columns=sample.columns)
        test_sample[target_col] = ar[test_rows, test_cols]
        test_sample['row'] = test_rows
        test_sample['col'] = test_cols
        test_sample['y'] = [int(ul_y + r * y_res) for r in test_rows]  
        test_sample['x'] = [int(ul_y + c * y_res) for c in test_cols]
        test_sample['obs_id'] = test_sample.index
        test_sample.set_index('obs_id')
        
        # Sample predictors
        test_sample['tile_id'] = [('000' + str(t))[-4:] for t in ar_tile[test_rows, test_cols]]
        test_sample['tile_str'] = test_sample.tile_id
    
    df_tile = test_sample.drop_duplicates('tile_id')[['tile_id', 'tile_str']]
    n_years = 1
    year = 2001
    n_files = (len(df_tile) * len(var_info[var_info.by_tile == 1])) * n_years +\
    n_years * (len(var_info[var_info.data_band < 0]) + len(var_info[var_info.data_band > 0]))
    
    point_dict = {p:g.index.tolist() for p, g in test_sample.groupby('tile_id')}
    t1 = time.time()
    df_yr = pd.DataFrame(index=test_sample.index)
    last_file = 1
    for var_name, var_row in var_info.iterrows(): #var_name is index col
        if var_name in test_sample.columns: # Test sample already exists
            continue
        # Get variables from row
        search_str  = var_row.search_str
        basepath    = var_row.basepath
        path_filter = var_row.path_filter
        data_type   = var_row.data_type
        by_tile      = var_row.by_tile
        path_filter = ''
        #data_band   = var_row.data_band 
        if np.isnan(var_row.nodata):
            nodata = None
        else:
            nodata = int(var_row.nodata)
        if int(var_row.data_band) < 0: 
            data_band = year - 1984 + 1
        else: 
            data_band = int(var_row.data_band)

        df_var, last_file = extract_var(year, var_name, by_tile, data_band,
                                         data_type, df_tile, test_sample, point_dict, basepath,
                                         search_str, path_filter, tx,
                                         last_file, n_files, nodata, False)
        df_yr[df_var.columns] = df_var  
    
    for c in df_yr.columns:
        test_sample[c] = df_yr[c]
    
    return test_sample[sample.columns.tolist() + var_info.index.tolist()]


def tree_accuracy(dt, test_sample, predictors, bins, classes, ar_ref, mask, r_nodata, target_col, bin_counts=None, temp_nodata=-9999, dtype=np.int16):
    
    
    class_dict = dict(zip(dt.classes_, classes))
    class_dict[r_nodata] = r_nodata
    predictions = dt.predict(predictors).astype(dtype)
    predictions = replace_val(predictions, class_dict)
    
    sample = pd.DataFrame({target_col: test_sample[target_col], 'prediction': predictions})
    cm = confusion_matrix_by_area(ar_ref, ar_ref, sample, r_nodata, r_nodata, mask, bins=bins, target_col=target_col, get_totals=False, total_counts=bin_counts, silent=True) #Using ar_ref for ar_p is only a problem because total counts are derived from ar_ref instead of a prediction map
    
    labels = ['%s_%s' % (l, u) for l, u in bins]
    p_acc = cm.ix['producer', labels]
    u_acc = cm.ix[labels, 'user']
    p_acc.set_axis(0, ['p_%s' % i.split('_')[1] for i in p_acc.index])
    u_acc.set_axis(0, ['u_%s' % i.split('_')[1] for i in u_acc.index])
    s = pd.concat([p_acc, u_acc])
    s['accuracy'] =  cm.ix['producer','user']
    s['kappa'] = cm.ix['producer', 'kappa']

    return s
    
    
def main(params_txt, pred_path=None, bins=None, r_nodata=255, temp_nodata=-9999, n_per_bin=200):
    
    inputs = read_params(params_txt)
    for var in inputs:
        exec ("{0} = str({1})").format(var, inputs[var])
    r_nodata = int(r_nodata)
    temp_nodata = int(temp_nodata)
    n_per_bin = int(n_per_bin)
    
    if type(forest) == str:
        try:
            with open(forest) as f:
                forest = pickle.load(f)
        except BaseException as e:
            print e.message
    
    if type(tile_mosaic) == str:
        ds_t = gdal.Open(tile_mosaic)
        ar_tile = ds_t.ReadAsArray()
        ds_t = None
    else:
        ar_tile = tile_mosaic
    
    ds_r = gdal.Open(ref_path)
    ar_ref = ds_r.ReadAsArray()
    tx = ds_r.GetGeoTransform()
    ds_r = None
    
    #bins = parse_bins(bins)
    if not bins:
        unique_vals = np.unique(ar_ref)
        bins = zip(unique_vals - 1, unique_vals)
    
    bin_counts = None
    if pred_path:
        # Get the number of pixels in each bin in the prediction map
        bin_counts = {}
        ds_p = gdal.Open(pred_path)
        ar_p = ds_p.ReadAsArray()
        for l, u in bins:
            bin_mask = (ar_p > l) & (ar_p <= u)
            bin_counts['%s_%s' % (l, u)] = ar_p[bin_mask].size
        ds_p = None
        ar_p = None
            
    var_info = pd.read_csv(var_txt, sep='\t', index_col='var_name')
    predict_cols = sorted(var_info.index)
    var_info.reindex(predict_cols)
    
    out_dir = os.path.dirname(out_txt)
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    if 'test_txt' in inputs:
        test_sample = pd.read_csv(test_txt)#, sep='\t', index+col='obs_id')
    elif 'train_txt' in inputs:
        train_sample = pd.read_csv(train_txt, sep='\t', index_col='obs_id')
        test_bn = os.path.basename(train_txt).replace('.txt', '_test.txt')
        test_txt = os.path.join(out_dir, test_bn)
        test_sample = None
        if os.path.exists(test_txt): # append to the existing
            test_sample = pd.read_csv(test_txt)
        test_sample = get_test_sample(train_sample, ar_ref, bins, n_per_bin, r_nodata, target_col, tx, ar_tile, var_info, test_sample)
        test_bn = os.path.basename(train_txt).replace('.txt', '_test.txt')
        test_txt = os.path.join(out_dir, test_bn)
        if os.path.exists(test_txt): # append to the existing 
            test_sample.to_csv(train_txt.replace('.txt', '_test.txt'), sep='\t')
    else:
        raise RuntimeError('Either test_txt or train_txt must be specified. Input parameters: \n\t' + '\n\t'.join(inputs.keys()))
    

    # Get predictors
    '''t1 = time.time()
    coords = [ tx[0], tx[3], tx[0] + (ar_ref.shape[0] * tx[1]), tx[3] + (ar_ref.shape[1] * tx[5])]
    tile_mask = ar_tile == 0
    ar_tile[tile_mask] = r_nodata
    tile_ids = np.unique(ar_tile)
    tile_strs = [str(tile) for tile in tile_ids if tile!=r_nodata]
    predictors = stem.get_predictors(var_info, tx, tile_strs, ar_tile, coords, tile_mask, temp_nodata, 1)
    predictor_mask = np.any(predictors==temp_nodata, axis=1)
    mask = tile_mask | (ar_ref == r_nodata) | predictor_mask.reshape(ar_ref.shape)
    print '\nGetting predictors: %.1f minutes\n' % ((time.time() - t1)/60)'''
    
    mask = (ar_ref == r_nodata) | (ar_tile == 255)
    predictors = test_sample[predict_cols]    
    
    # Loop through each tree and get accuracies
    per_tree = []
    n_trees = forest.n_estimators
    classes = forest.classes_
    t1 = time.time()
    for i, dt in enumerate(forest.estimators_):
        #print i
        #accuracy = tree_accuracy(dt, predictors, bins, classes, ar_ref, ar_tile, predictor_mask, mask, r_nodata, target_col)
        accuracy = tree_accuracy(dt, test_sample, predictors, bins, classes, ar_ref, mask, r_nodata, target_col, bin_counts)
        per_tree.append(accuracy)
        msg_inputs = i + 1, n_trees, (i + 1)/float(n_trees) * 100, (time.time() - t1)/60
        sys.stdout.write('\rFinished %s of %s trees (%.1f%%). Cumulative time: %.1f minutes' % msg_inputs)
        sys.stdout.flush()
    
    df = pd.DataFrame(per_tree)
    df.to_csv(out_txt, sep='\t', index=False)
    
    print '\n\nText file written to', out_txt
    

if __name__ == '__main__':
    
    sys.exit(main(*sys.argv[1:]))
    