# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 10:15:46 2017

@author: braatenj
"""

import json
import subprocess
import multiprocessing
import os
from glob import glob


def run_cmd(cmd):
  print(cmd)  
  return subprocess.call(cmd, shell=True)

def get_tile_id_and_coords(feature):          
  xmax = feature['properties']['xmax']
  xmin = feature['properties']['xmin']
  ymax = feature['properties']['ymax']
  ymin = feature['properties']['ymin']
  coords = ' '.join([str(coord) for coord in [xmin,ymax,xmax,ymin]])
   
  # prepend 0's to tileID so they all have 4 digits
  tileID = str(feature['properties']['id'])
  zeros = '0'*(4-len(tileID))
  tileID = zeros+tileID        
  
  return (coords, tileID) 


#vrtFile =  '/vol/v1/proj/stem_improv_paper/test_regions/tiles/staging/vrt/ltee_pecora_multiwindow_20171102_04150715_17_all_stack.vrt' 
#tileFile = '/vol/v1/proj/stem_improv_paper/test_regions/tiles/staging/regions/region_17.geojson'

vrtDir = '/vol/v1/proj/stem_improv_paper/test_regions/tiles/staging/vrt/'
tileDir = '/vol/v1/proj/stem_improv_paper/test_regions/tiles/staging/regions/'

#################################################################################################################33333

vrtFiles = glob(vrtDir+'*17_all_stack.vrt')

tileFiles = []
for vrtFile in vrtFiles:
  region = os.path.basename(vrtFile)[42:44]
  tileFile = glob(tileDir+'*'+region+'*.geojson')[0]
  tileFiles.append(tileFile)

outDir = '/vol/v1/proj/stem_improv_paper/test_regions/tiles/tiles_seg'
indexID = 'nbr'
nVert = 7
startYear = 1984
endYear = 2016
proj = 'EPSG:5070'


for vrtFile, tileFile in zip(vrtFiles, tileFiles):
  
  name = os.path.basename(vrtFile)[0:41]
  region = os.path.basename(vrtFile)[42:44]
  
  
  # make sure path parts are right
  if outDir[-1] != '/':
    outDir += '/'
    
  # define the list of replacement types in the new out images    
  outTypes = ['vert_yrs.bsq',
              'vert_src.bsq',
              'vert_fit.bsq',
              'ftv_'+indexID+'.bsq',            
              'ftv_tcb.bsq',
              'ftv_tcg.bsq',
              'ftv_tcw.bsq']
  
  # make a list of band ranges for each out type
  vertStops = []
  for vertType in range(4):  
    vertStops.append(vertType*nVert+1)
  
  nYears = (endYear - startYear) + 1
  ftvStops = []  
  for ftvType in range(1,5):  
    ftvStops.append(ftvType*nYears+vertStops[-1])
    
  bandStops = vertStops+ftvStops
  bandRanges = [range(bandStops[i],bandStops[i+1]) for i in range(len(bandStops)-1)]
  
  
  # load the tile features
  with open(tileFile) as f:
    features = json.load(f)['features']
   
   
  # make a list of all the gdal_translate commands needed for the ee conus chunk
  cmdList = []
  for feature in features:   
    coords, tileID = get_tile_id_and_coords(feature)
    
    regionOutDir = outDir+'region_'+region
    if not os.path.isdir(regionOutDir):
      os.mkdir(regionOutDir)    
    tileOutDir = regionOutDir+'/'+tileID
    if not os.path.isdir(tileOutDir):
      os.mkdir(tileOutDir)
    
    for i in range(len(outTypes)):
      outFile = tileOutDir+'/'+tileID+'_'+name+'_'+outTypes[i]
      bands = ' -b '+' -b '.join([str(band) for band in bandRanges[i]])
      cmd = 'gdal_translate -q --config GDAL_DATA "/usr/lib/anaconda/share/gdal" -of ENVI -a_srs ' + proj + bands + ' -projwin '+coords+' '+ vrtFile + ' ' + outFile    
      cmdList.append(cmd)  
  
  
  # run the commands in parallel
  processes=len(features)
  if processes > 10: processes=10
  pool = multiprocessing.Pool(processes=processes)
  pool.map(run_cmd, cmdList)  
  pool.close()