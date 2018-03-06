import os, sys, time, fnmatch
from glob import glob
from multiprocessing import Pool
from osgeo import gdal, ogr
import numpy as np


from lthacks import createMetadata#, array_to_raster


def write_raster(data, filename, prj, origin, driver, dtype=gdal.GDT_Int16, pixel_size=30):
    # origin = [X, Y] offset
    # srs = template.GetProjectionRef()

    bands, rows, cols = data.shape

    #driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(filename, cols, rows, bands, dtype)

    outRaster.SetGeoTransform((origin[0], pixel_size, 0, origin[1], 0, -pixel_size))
    outRaster.SetProjection(prj)

    for i in range(bands):
        band = outRaster.GetRasterBand(i+1)
        band.WriteArray(data[i,:,:])
    print 'Raster written to ', filename
    

def make_history_vars(vert_year_raster, ftv_raster, out_dir=None, make_ts=False, nodata=-9999):
    
    if not make_ts:
      print('    SKIPPING')
      return
    
    ds_years = gdal.Open(vert_year_raster)
    vert_years = ds_years.ReadAsArray()
    tx = ds_years.GetGeoTransform()
    prj = ds_years.GetProjection()
    driver = ds_years.GetDriver()
    ds_years = None
    
    #ds_ftv = gdal.Open(ftv_raster)
    #ftv_vals = ds_ftv.ReadAsArray()
    #ds_ftv = None
    
    n_verts, ysize, xsize = vert_years.shape
    
    start_year = 1984
    end_year = 2016    
    n_years = end_year-start_year+1

    #Years from vertex
    if make_ts:
        time_since = np.full((n_years, xsize, ysize), nodata, np.int16)
    #Change from vertex
    #history = np.full((n_years, xsize, ysize), nodata, np.int16)

    for x in xrange(xsize):
        for y in xrange(ysize):
            
            if (vert_years[:, x, y] == 0).all():
                if make_ts:
                    time_since[:, x, y] = 0
                #history[:, x, y] = 0
                continue

            s0 = vert_years[0,x,y]-start_year
            #history[:s0 + 1, x, y] = 0 #In case 1st year isn't 1984
            if make_ts:
                time_since[:s0 + 1, x, y] = range(s0 + 1)
            
            max_year = 0
            for n in range(1,7):
                s1 = vert_years[n, x, y]-start_year
                if (s1 == -start_year) or (s1 > n_years):
                    # Fill to end in case last year isn't 2016
                    if max_year < n_years - 1:
                        if make_ts:
                            time_since[s0 + 1 : n_years, x, y] = range(1, n_years - s0)
                        #history[s0 + 1 : n_years, x, y] = 0
                    break
                if make_ts:
                    time_since[s0 + 1 : s1 + 1, x, y] = range(1, s1 - s0 + 1)
                #history[s0 + 1 : s1 + 1, x, y] = ftv_vals[s0 + 1 : s1 + 1, x, y] - ftv_vals[s0, x, y]
                
                max_year = max([max_year, s1])
                s0 = s1
            #import pdb; pdb.set_trace()

    if not out_dir:
        out_dir = os.path.dirname(ftv_raster)
    basename = os.path.basename(ftv_raster)
    var_name = basename.split('_')[-1]
    out_basename = basename.replace('ftv_' + var_name, 'delta_%s' % var_name)
    out_path = os.path.join(out_dir, out_basename)
    origin = tx[0], tx[3]
    #write_raster(history, out_path, prj, origin, driver)
    #desc = ('Change from vertex calculated per pixel.\n' +\
    #        '\tVert year raster: {0}\n' + \
    #        '\tFitted val raster: {1}').format(vert_year_raster, ftv_raster)
    #createMetadata(sys.argv, out_path, description=desc)
    
    if make_ts:
        ext = var_name.split('.')[-1]
        out_basename = basename.replace('ftv_' + var_name, 'time_since.%s' % ext)
        out_path = os.path.join(out_dir, out_basename)
        write_raster(time_since, out_path, prj, origin, driver, dtype=gdal.GDT_Byte)
        desc = 'Time since last vertex calculated per pixel.\n' +\
                '\tVert year raster: %s\n' % vert_year_raster
        createMetadata(sys.argv, out_path, description=desc)


def par_make_history_vars(args):
    
    n_files, (vert_year_raster, ftv_raster, out_dir, make_ts, nodata, n) = args
    
    t0 = time.time()
    print 'Making history for %s of %s files...' % (n, n_files)
    try:
        make_history_vars(vert_year_raster, ftv_raster, out_dir, make_ts, nodata)
    except:
        print vert_year_raster
        sys.exit()
    print 'Time for file %s of %s: %.1f minutes' % (n, n_files, (time.time() - t0)/60)


def main(search_dir, out_dir=None, make_time_since=True, nodata=-9999):
    ######################################################    
    #search_dir='/vol/v1/proj/stem_improv_paper/test_regions/tiles/tiles_seg' 
    #out_dir=None 
    #make_time_since=True
    #nodata=-9999    
    ######################################################
    
    
    t0 = time.time()

    # Set up args per pixel to pass to pool.map()
    args = []
    j = 0
    #tile_dirs = sorted(os.listdir(search_dir)) # jdb comment out 11/5/17 - added next 5 lines
    rasters = []
    for root, dirnames, filenames in os.walk(search_dir):
        for filename in fnmatch.filter(filenames, '*ftv_nbr.bsq'):
            rasters.append(os.path.join(root, filename))

    tile_dirs = list(set([os.path.dirname(raster) for raster in rasters]))

    #n_tiles = len(tile_dirs)
    this_out_dir = out_dir
    for i, tile_dir in enumerate(tile_dirs):
        ######################################################        
        #i=0
        #tile_dir = tile_dirs[0]
        ######################################################
    
    
    
    
    
       #tile_dir = os.path.join(search_dir, t) # jdb comment out 11/5/17 - not needed, it is now in the for line
        ftv_rasters = glob(os.path.join(tile_dir, '*ftv*.bsq'))
        # vert_year_raster = glob(os.path.join(tile_dir, '*yrs.bsq'))[0] # jdb comment out 6/29/17 added in line just below to replace
        timeWindows = list(set([os.path.basename(ftv_raster)[38:46] for ftv_raster in ftv_rasters]))       
        timeWindwosDone = [];        
        make_ts = make_time_since
        if not out_dir:
            this_out_dir = tile_dir
        for ftv_raster in ftv_rasters:          
            
            make_ts = False
            thisWindow = os.path.basename(ftv_raster)[38:46]
            if thisWindow not in timeWindwosDone:
              timeWindwosDone.append(thisWindow)
              make_ts = True
            
            ftv_bname = os.path.basename(ftv_raster)           
            ftv_dname = os.path.dirname(ftv_raster)            
            vert_year_raster = os.path.join(ftv_dname,ftv_bname[0:ftv_bname.index('ftv_')]+'vert_yrs.bsq') # jdb added 6/29/17 to match different ftv run in same dir
            args.append([vert_year_raster,
                         ftv_raster,
                         this_out_dir,
                         make_ts,
                         nodata,
                         j + 1])
            j += 1
            #make_ts = False # Set it so we only make it once with the first ftv raster found in this tile
    n_files = len(args)
    args = [(n_files, a) for a in args]
    #import pdb; pdb.set_trace()
    #args = args[:35]
    pool = Pool(10)
    pool.map(par_make_history_vars, args, 1)#'''
    
    print '\nTotal time for fitting: %.1f minutes\n' % ((time.time() - t0)/60)
    

if __name__ == '__main__':
    
    sys.exit(main(*sys.argv[1:]))#'''
    


#out_dir = '/home/server/student/homes/shooper/delete'
#search_dir = '/vol/v2/conus_tiles/tiles'

#main(search_dir) #out_dir=out_dir, make_time_since=True)'''

