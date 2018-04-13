# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 13:48:29 2017

@author: shooper

Make masks/reference rasters clipped to the shape of each study region

Usage:
    make_region_masks.py <region_path> <tile_path> <reference_path> <out_dir> [--id_field=<str>] [--ref_basename=<str>]

    region_path - path to vector file where each feature is the boundary of a region
    tile_path - path to vector file where each feature is a tile
    reference_path - path to reference raster (e.g. NLCD)
    out_dir - top level output directory

Options:
    -h, --help              show this screen
    --id_field=<str>        field in region_path vector containing unique IDs[default: region_id]
    --ref_basename=<str>    name of reference dataset (e.g., nlcd, canopy) [default: nlcd]

example:
make_region_masks.py /vol/v1/proj/stem_improv_paper/vector/regions/study_regions.shp /vol/v1/proj/stem_improv_paper/vector/regions/study_region_tiles.shp /vol/v1/general_files/datasets/spatial_data/nlcd/nlcd_2001_v2/nlcd_2001_landcover_clipped_to_conus_tiles.tif /vol/v1/proj/stem_improv_paper/region_models

"""

import os
import sys
import docopt
from osgeo import gdal, ogr
import re
import pandas as pd
import numpy as np

sys.path.append('/vol/v2/stem/stem-git/scripts')
from stem import coords_to_shp
from mosaic_by_tsa import kernel_from_shp, get_offset_array_indices
from lthacks import attributes_to_df, array_to_raster


def main(region_path, tile_path, reference_path, out_dir, id_field='region_id', ref_basename='nlcd'):
    
    df = attributes_to_df(region_path) 
    tile_info = attributes_to_df(tile_path)
    tile_info['ul_x'] = tile_info.xmin
    tile_info['lr_x'] = tile_info.xmax
    tile_info['ul_y'] = tile_info.ymax
    tile_info['lr_y'] = tile_info.ymin
    
    _, vector_ext = os.path.splitext(region_path)
    region_ids = df[id_field].unique()
    n_regions = len(region_ids)
    
    region_ds = ogr.Open(region_path)
    region_lyr = region_ds.GetLayer()
    
    
    for i, r_id in enumerate(region_ids):
        print 'Making region dir for %s (%s of %s)' % (r_id, i, n_regions)
        df_r = df[df.region_id == r_id]
        id_str = ('0' + str(r_id))[-2:]
        
        fid = df_r.index[0]
        region_feature = region_lyr.GetFeature(fid)
        xmin, xmax, ymin, ymax = region_feature.GetGeometryRef().GetEnvelope()
        region_feature.Destroy()
        df_r['ul_x'] = xmin
        df_r['lr_x'] = xmax
        df_r['ul_y'] = ymax
        df_r['lr_y'] = ymin
        clip_coords = df_r.loc[fid, ['ul_x', 'lr_x', 'ul_y', 'lr_y']]
        
        region_dir = os.path.join(out_dir, 'region_%s' % id_str)
        if not os.path.exists(region_dir):
            os.mkdir(region_dir)
        
        # Make a shapefile of the tiles
        out_vector = os.path.join(region_dir, 'tile_{0}{1}'.format(id_str, vector_ext))
        if not os.path.exists(out_vector):
            ''' switch to selection by min/max of coords '''
            region_tiles = tile_info[tile_info[id_field] == r_id]
            coords_to_shp(region_tiles, region_path, out_vector)
        
        # Make a map of reference NLCD
        ds = gdal.Open(out_vector.replace(vector_ext, '.tif'))
        mask = ds.ReadAsArray() == 255
        ds = None
        nlcd_year = re.search('\d\d\d\d', reference_path).group() # finds the first one (potentially buggy)
        out_ref_map = os.path.join(region_dir, '%s_%s_%s.tif' % (ref_basename, nlcd_year, id_str))
        if not False:#os.path.exists(out_ref_map):
            ref_ds = gdal.Open(reference_path)
            ref_tx = ref_ds.GetGeoTransform()
            ref_shape = ref_ds.RasterYSize, ref_ds.RasterXSize
            
            col_off = (ref_tx[0] - clip_coords.ul_x)/ref_tx[1]
            row_off = (ref_tx[3] - clip_coords.ul_y)/ref_tx[5]
            n_cols = abs((clip_coords.ul_x - clip_coords.lr_x)/ref_tx[1])
            n_rows = abs((clip_coords.ul_y - clip_coords.lr_y)/ref_tx[1])
            
            ar_inds, ref_inds = get_offset_array_indices((n_rows, n_cols), ref_shape, (row_off, col_off))
            ref_n_cols = ref_inds[1] - ref_inds[0]
            ref_n_rows = ref_inds[3] - ref_inds[2]
            
            ar_ref = ref_ds.ReadAsArray(ref_inds[2], ref_inds[0], ref_n_cols, ref_n_rows)
            ar = np.full((n_rows, n_cols), 255)
            ar[ar_inds[0]:ar_inds[1], ar_inds[2]:ar_inds[3]] = ar_ref
            ar[mask] = 255
            
            tx = clip_coords.ul_x, 30, 0, clip_coords.ul_y, 0, -30
            prj = ref_ds.GetProjection()
            driver = gdal.GetDriverByName('gtiff')
            array_to_raster(ar, tx, prj, driver, out_ref_map, nodata=255)
            
        # Make a clipped raster of the tiles
        out_raster = out_vector.replace(vector_ext, '.tif')
        if not os.path.exists(out_raster):
            tiles = ogr.Open(tile_path)
            tile_lyr = tiles.GetLayer()
            tx = clip_coords.ul_x, 30, 0, clip_coords.ul_y, 0, -30
            tile_array, _ = kernel_from_shp(tile_lyr, clip_coords, tx, 255, val_field='name')
            tile_array[ar == 255] = 255
            driver = gdal.GetDriverByName('gtiff')
            prj = tile_lyr.GetSpatialRef().ExportToWkt()
            array_to_raster(tile_array, tx, prj, driver, out_raster, nodata=255)
            tiles.Destroy()


if __name__ == '__main__':
    
    try:
        cl_args = docopt.docopt(__doc__)
    except docopt.DocoptExit as e:
        import pdb; pdb.set_trace()
        print e.message
    
    # get rid of extra characters and 'help' entry from doc string 
    args = {re.sub('[<>-]*', '', k): v for k, v in cl_args.iteritems()
            if k not in ['--help','-h']} 
    sys.exit(main(**args))
    #sys.exit(main(*sys.argv[1:]))
        
    
    
    
    


