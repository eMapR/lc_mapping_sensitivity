
import os
import time
import shutil
import fnmatch
import warnings
import subprocess
import multiprocessing
import sys
from datetime import datetime
from osgeo import ogr

from lthacks import createMetadata


def getInfo(name):
    pieces = name.split('-')[0:6]
    info = {'key': pieces[0],
            'value': pieces[1],
            'indexID': pieces[2].lower(),
            'nVert': int(pieces[3]),
            'startYear':int(pieces[4][0:4]),
            'endYear':int(pieces[4][4:8]),
            'dateRange':pieces[5],
            'processingDate': str(datetime.now().date()).replace('-','')
            }
    
    return info
    

def call(callInfo):
    outFile, command, count, nCommands, startTime = callInfo
    cumTime = time.time() - startTime
    formatObj = count, nCommands, float(count)/nCommands * 100, cumTime/60
    sys.stdout.write('\rProcessing %s of %s files (%.1f%%) | cum. time: %.1f mins.' % formatObj)
    sys.stdout.flush()
    subprocess.call(command, shell=True)
    desc = 'LandTrendr data created with Google Earth Engine. Filename is in the form {version}_{index}_{dateRange}_{tileID}_{nVertices}_{processingDate}_{LTdataType}_{LTdataContent}.tif'
    createMetadata(sys.argv, outFile, description=desc)


def main(chunkDir, outDir, clipFile,  njobs=0, tileIdField='name', proj='EPSG:5070'):
    
    t0 = time.time()
    
    njobs = int(njobs)
    
    if not os.path.isdir(outDir):
        warnings.warn('outDir does not exist... Creating outDir %s' % outDir)
        os.mkdir(outDir)
    
    # define the list of replacement types in the new out images    
    outTypes = ['vert_yrs.tif',
                'vert_src.tif',
                'vert_fit.tif',
                'ftv_idx.tif',
                'ftv_tcb.tif',
                'ftv_tcg.tif',
                'ftv_tcw.tif',
                'ftv_ndvi.tif',
                'ftv_ndsi.tif']
    
    # find the tif chunks
    tifs = []
    for root, dirnames, filenames in os.walk(chunkDir):
        for filename in fnmatch.filter(filenames, '*.tif'):
            tifs.append(os.path.join(root, filename))    
    
    # set the unique names 
    names = list(set(['-'.join(fn.split('-')[0:6]) for fn in tifs]))
    nTiles = 9
    nCommands = len(names) * len(outTypes) * nTiles
    # loop through each unique names, find the matching set, merge them as vrt, and then decompose them
    c = 1 #counter for total number of commands submitted
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
        info = getInfo(runName)

        # make a list of tile tifs
        vrtFile = os.path.join(chunkDir, runName + '.vrt')
        tileListFile = vrtFile.replace('.vrt', '_filelist.txt')
        with open(tileListFile, 'w') as tileList:
            for match in matches:
                tileList.write(match + '\n')
    
        # create vrt
        cmd = 'gdalbuildvrt -input_file_list %s %s' % (tileListFile, vrtFile)
        subprocess.call(cmd, shell=True)
    
        # make a list of band ranges for each out type
        vertStops = []
        for vertType in range(4):  
            vertStops.append(vertType*info['nVert']+1)

        nYears = (info['endYear'] - info['startYear']) + 1
        ftvStops = []  
        for ftvType in range(1, len(outTypes)-2):  
            ftvStops.append(ftvType*nYears+vertStops[-1])

        bandStops = vertStops+ftvStops
        bandRanges = [range(bandStops[i],bandStops[i+1]) for i in range(len(bandStops)-1)]

        
        # Check that the directory structure exists. If not, build it.
        indexDir = os.path.join(outDir, '%s_%s' % (info['indexID'], info['dateRange']))
        if not os.path.isdir(indexDir):
            os.mkdir(indexDir)
        
        # Open the tile layer and loop through each tile within this region
        dataset = ogr.Open(clipFile)
        layer = dataset.GetLayer()
        filterString = "%s = '%s'" % (info['key'], info['value'])
        layer.SetAttributeFilter(filterString)
        feature = layer.GetNextFeature()
        i = 1
        while feature:
            tileID = feature.GetField(tileIdField)
            tileDir = os.path.join(indexDir, tileID)
            if not os.path.isdir(tileDir):
                os.mkdir(tileDir)

            # format the exent as -projwin arguments for gdal translate
            extent = feature.GetGeometryRef().GetEnvelope()
            projwin = '{} {} {} {}'.format(extent[0], extent[3], extent[1], extent[2])    
         
            # make a list of all the gdal_translate commands needed for the ee conus chunk
            if not os.path.isdir(outDir):
                os.mkdir(outDir)
            # loop through the datasets and pull them out of the mega stack and write them to the define outDir
            info['tile'] = tileID
            for i in range(len(outTypes)):   
                #outBname = runName+'-'+outTypes[i]
                info['outType'] = outTypes[i]
                outBname = 'LTV3_{indexID}_{dateRange}_{tile}_{nVert}_{processingDate}_{outType}'.format(**info)
                outFile = os.path.join(tileDir, outBname)
                if os.path.exists(outFile):
                    continue
                bands = ' -b '+' -b '.join([str(band) for band in bandRanges[i]])
                cmd = 'gdal_translate -q --config GDAL_DATA "/usr/lib/anaconda/share/gdal" -of GTiff -a_srs ' + proj + bands + ' -projwin '+projwin+' '+ vrtFile + ' ' + outFile
                commands.append((outFile, cmd, c, nCommands, t0))
                c += 1
            feature = layer.GetNextFeature()
            i += 1
        if i != nTiles:
            raise RuntimeError('%s is not a valid filter string for clipFile %s' % (filterString, clipFile))
    
    if njobs > 1:
        pool = multiprocessing.Pool(njobs)
        pool.map(call, commands, chunksize=1)
        pool.close()
        pool.join()
        
    else:
        for cmd in commands: subprocess.call(cmd, shell=True)
    
    # Delete source files
    shutil.rmtree(chunkDir)
    
    print '\n\nTotal processing time: %.1f minutes' % ((time.time() - t0)/60)


if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:]))