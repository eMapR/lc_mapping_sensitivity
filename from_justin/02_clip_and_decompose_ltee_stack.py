
import os
import time
import fnmatch
import subprocess
import multiprocessing
import sys
from osgeo import ogr


#def run_cmd(cmd):
#  print(cmd)  
#  return subprocess.call(cmd, shell=True)

def get_info(name):
  pieces = name.split('-')[0:6]
  return {'key': pieces[0],
          'value': pieces[1],
          'indexID': pieces[2],
          'nVert': int(pieces[3]),
          'startYear':int(pieces[4][0:4]),
          'endYear':int(pieces[4][4:8])}


########### added 2/1/18 by smh ###########
def call(command):
    outFile, command = command
    print('Processing ' + outFile)
    subprocess.call(command, shell=True)# pretty sure you want shell=False, the default
##########################################

t0 = time.time()
# get the arguments
arg = sys.argv
chunkDir = arg[1]
outDir = arg[2]  
clipFile = arg[3]
#key = arg[4]
#value = arg[5] 
#indexID = arg[6]
#nVert = int(arg[7])
#startYear = int(arg[8])
#endYear = int(arg[9])
#proj = arg[10]
#"""


#chunkDir = '/vol/v1/proj/stem_improv_paper/raster/prep'
#outDir = '/vol/v1/proj/stem_improv_paper/raster/r0701' 
#clipFile = '/vol/v1/proj/stem_improv_paper/vector/regions/study_regions.shp'
#key = 'epa_code'
#value = 'r0701'
#indexID = 'nbr'
#nVert = 7
#startYear = 1984
#endYear = 2016


######################################################################

# define the projection
proj = 'EPSG:5070'

njobs = 9

# make sure path parts are right
if chunkDir[-1] != '/':
  chunkDir += '/'
if outDir[-1] != '/':
  outDir += '/'


# find the tif chunks
tifs = []
for root, dirnames, filenames in os.walk(chunkDir):
  for filename in fnmatch.filter(filenames, '*.tif'):
    tifs.append(os.path.join(root, filename))

# set the unique names 
names = list(set(['-'.join(fn.split('-')[0:6]) for fn in tifs])) 


# loop through each unique names, find the matching set, merge them as vrt, and then decompose them
commands = [] #record all commands to run
for name in names:
  #name = names[0]  
  runName = os.path.basename(name)
  # find the files that belong to this set  
  matches = []
  for tif in tifs:   
    if name in tif:
      matches.append(tif)

  # get info about the GEE run
  info = get_info(runName)

  # make a list of tile tifs
  vrtFile = chunkDir+runName+'.vrt'
  tileListFile = vrtFile.replace('.vrt', '_filelist.txt')
  tileList = open(tileListFile, 'w')
  for match in matches:
    tileList.write(match+'\n')
  tileList.close()


  # create vrt
  cmd = 'gdalbuildvrt -input_file_list '+tileListFile+' '+vrtFile
  subprocess.call(cmd, shell=True)


  # define the list of replacement types in the new out images    
  outTypes = ['vert_yrs.bsq',
              'vert_src.bsq',
              'vert_fit.bsq',
              'ftv_idx.bsq',            
              'ftv_tcb.bsq',
              'ftv_tcg.bsq',
              'ftv_tcw.bsq',
              'ftv_ndvi.bsq',
              'ftv_ndsi.bsq']


  # make a list of band ranges for each out type
  vertStops = []
  for vertType in range(4):  
    vertStops.append(vertType*info['nVert']+1)
  
  nYears = (info['endYear'] - info['startYear']) + 1
  ftvStops = []  
  for ftvType in range(1,len(outTypes)-2):  
    ftvStops.append(ftvType*nYears+vertStops[-1])
    
  bandStops = vertStops+ftvStops
  bandRanges = [range(bandStops[i],bandStops[i+1]) for i in range(len(bandStops)-1)]


  # get the driver from the inShape file ext
  ext = str.lower(os.path.splitext(clipFile)[-1])
  drivers = {'.shp'    :'ESRI Shapefile', 
             '.geojson': 'GeoJSON'}           
  driver = ogr.GetDriverByName(drivers[ext])
  
  
  # read in the inShape file and get the extent of the feature defined by key and value
  dataSource = driver.Open(clipFile, 0)
  layer = dataSource.GetLayer()
  layer.SetAttributeFilter(info['key']+" = '"+info['value']+"'")
  feature = layer.GetNextFeature()  
  extent = feature.GetGeometryRef().GetEnvelope()
  
  # format the exent as -projwin arguments for gdal translate
  projwin = '{} {} {} {}'.format(extent[0], extent[3], extent[1], extent[2])    
  #projwin = '-2166434 2699450 -2164169 2697218'  
  
  
  # make a list of all the gdal_translate commands needed for the ee conus chunk
  if not os.path.isdir(outDir):
    os.mkdir(outDir)
  
  
  # loop through the datasets and pull them out of the mega stack and write them to the define outDir
  ##### added 2/1/18 by smh to run concurrent processes ######################
  for i in range(len(outTypes)):   
    outBname = runName+'-'+outTypes[i]   
    outFile = outDir+outBname
    if os.path.exists(outFile):
        continue
    #print(outBname)
    bands = ' -b '+' -b '.join([str(band) for band in bandRanges[i]])
    cmd = 'gdal_translate -q --config GDAL_DATA "/usr/lib/anaconda/share/gdal" -of ENVI -a_srs ' + proj + bands + ' -projwin '+projwin+' '+ vrtFile + ' ' + outFile
    commands.append((outFile, cmd))

# safer to put inside if __name__ == '__main__' block:
#   https://docs.python.org/2/library/multiprocessing.html#windows
if __name__ == '__main__':
  if njobs > 1:
    pool = multiprocessing.Pool(njobs)
    pool.map(call, sorted(commands))#sorting groups cmds by run type
    pool.close()
    pool.join()
    print 'Processing time: %.1f days' % ((time.time() - t0)/86400)
  else:
    for cmd in commands: subprocess.call(cmd, shell=True)


##########################################################################
  
  '''for i in range(len(outTypes)):   
    outBname = runName+'-'+outTypes[i]   
    print(outBname)
    outFile = outDir+outBname
    bands = ' -b '+' -b '.join([str(band) for band in bandRanges[i]])
    cmd = 'gdal_translate -q --config GDAL_DATA "/usr/lib/anaconda/share/gdal" -of ENVI -a_srs ' + proj + bands + ' -projwin '+projwin+' '+ vrtFile + ' ' + outFile       
    subprocess.call(cmd, shell=True)'''
"""
for thisFile in os.listdir(chunkDir):
  thisFilePath = os.path.join(chunkDir, thisFile)
  os.unlink(thisFilePath)
"""

