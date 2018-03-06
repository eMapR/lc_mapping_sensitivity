# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 10:29:40 2017

@author: shooper

Finds all contiguous features and assigns a unique ID preserved in attribute
table of out_path
"""

import sys
import os
import ogr
import numpy as np

ogr.UseExceptions()

def main(in_path, out_dir=None):
    
    ds = ogr.Open(in_path)
    lyr_in = ds.GetLayer()
    
    # For each feature, find all features touching it
    touching = {}
    fids = [f.GetFID() for f in lyr_in]
    print 'Finding adjacent features...\n'
    for i in fids:
        #fids.remove(i) # Don't search for this feature anymore.
        search_feature = lyr_in.GetFeature(i)
        try:
            search_geom = search_feature.GetGeometryRef()
        except:
            import pdb; pdb.set_trace()
        
        # Search all other features to see if they touch search_feature
        this_touching = [i]
        for j in fids:
            if j == i:
                continue
            this_feature = lyr_in.GetFeature(j)
            this_geom = this_feature.GetGeometryRef()
            if this_geom.Touches(search_geom):
                this_touching.append(j)
            this_feature = None
        touching[i] = this_touching
    
    # Go through each feature's list and find any features contiguous with
    #   features in their lists and so on
    region_id = 0
    regions = {}
    closed = False
    searched = []
    print '\nFinding regions...'
    for i, adjacent in touching.iteritems():
        if i in searched:
            continue

        # Add all features to a working list from any of the features in adjacent
        these_adjacent = []#adjacent #Add these just in case they're not touching other feature than i
        for j in adjacent:
            these_adjacent.extend(touching[j])
            searched.append(j) # Also mark each feature as searched
        these_adjacent.extend(adjacent)
        
        # Get a unique list of features from the working list. If any features
        #   in the unique list have not yet been searched, add their contigous
        #   features to the working list and check again. Repeat until all 
        #   features touching any other feature in the working list have been
        #   searched.
        while not closed:
            unique = np.unique(these_adjacent)
            not_searched = [f for f in unique if f not in searched]
            if len(not_searched) == 0:
                regions[region_id] = unique
                closed = True
                region_id += 1
            else:
                for k in not_searched:
                    these_adjacent.extend(touching[k])
                    searched.append(k)   
        closed = False # Reset for next region
    
    # Make output dataset and add field for region id
    print '%s contiguous regions found. Saving regions...' % len(regions)
    if out_dir:
        driver = ds.GetDriver()
        srs = lyr_in.GetSpatialRef()
        lyr_in_def = lyr_in.GetLayerDefn()
        basename = os.path.basename(in_path)
        _, ext = os.path.splitext(basename)
        out_path = os.path.join(out_dir, basename.replace(ext, '_contiguous%s' % ext))
        ds_out = driver.CreateDataSource(out_path)
        #except: print 'Could not create shapefile with out_shp: \n', out_path
        lyr_out = ds_out.CreateLayer(os.path.basename(out_path)[:-4], srs, geom_type=lyr_in.GetGeomType())
        
        for i in range(lyr_in_def.GetFieldCount()):
            field_def = lyr_in_def.GetFieldDefn(i)
            lyr_out.CreateField(field_def)
        lyr_out.CreateField(ogr.FieldDefn('region_id', ogr.OFTInteger))
        
        lyr_out_def = lyr_out.GetLayerDefn()
        for r_id, f_ids in regions.iteritems():
            for fid in f_ids:
                feat_in = lyr_in.GetFeature(fid)
                feat_out = ogr.Feature(lyr_out_def)
                feat_out.SetFID(fid)
                feat_out.SetField('region_id', r_id)
                for i in range(lyr_in_def.GetFieldCount()):
                    feat_out.SetField(lyr_in_def.GetFieldDefn(i).GetName(), feat_in.GetField(i))
                geom = feat_in.GetGeometryRef()
                feat_out.SetGeometry(geom.Clone())
                lyr_out.CreateFeature(feat_out)
                feat_out.Destroy()
                feat_in.Destroy()
            
        ds_out.Destroy()
        print '\nNew %s created at: %s' % (driver.GetName(), out_path)
    
    else:
        lyr_in.CreateField(ogr.FieldDefn('region_id', ogr.OFTInteger))
        for r_id, f_ids in regions.iteritems():
            for fid in f_ids:
                feat_in = lyr_in.GetFeature(fid)
                feat_in.SetField('redion_id', r_id)
                feat_in.Destroy()
        ds.Destroy()
        print '\nIDs of contiguous regions added to input dataset: ', in_path


if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:]))
        
                    