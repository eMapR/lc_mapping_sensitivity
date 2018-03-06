# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 10:02:43 2017

@author: shooper

Wall-to-wall confusion matrix
"""

import os
import sys
from osgeo import gdal
import numpy as np
import pandas as pd

sys.path.append('/vol/v2/stem/stem-git/scripts')
from evaluation.evaluation import confusion_matrix_by_area
from lthacks import createMetadata


def main(ref_path, pred_path, out_txt, target_col, ref_nodata, pred_nodata):
    
    ref_nodata = int(ref_nodata)
    pred_nodata = int(pred_nodata)
    
    ds_t = gdal.Open(ref_path)
    ar_t = ds_t.ReadAsArray()
    
    ds_p = gdal.Open(pred_path)
    ar_p = ds_p.ReadAsArray()
    
    sample = pd.DataFrame({target_col: ar_t.ravel()})
    sample['prediction'] = ar_p.ravel()
    rows, cols = np.indices(ar_t.shape)
    sample['row'] = rows.ravel()
    sample['col'] = cols.ravel()
    sample.drop(sample.index[(sample[target_col]==ref_nodata) | (sample.prediction == pred_nodata)])
    
    unique_vals = np.unique(ar_t[ar_t != ref_nodata])
    bins = zip(unique_vals - 1, unique_vals)
    
    cm, cm_smp = confusion_matrix_by_area(ar_p, ar_t, sample, pred_nodata, ref_nodata, bins=bins, target_col=target_col, out_txt=out_txt, match='best')
    #cm.to_csv(out_txt, sep='\t')
    
    accuracy = cm.ix['producer','user']
    kappa = cm.ix['producer', 'kappa']
    
    print 'overall accuracy .............. ', accuracy
    print 'kappa ......................... ', kappa
    
    print 'Text file written to', out_txt


if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:]))