# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 11:27:18 2018

@author: shooper
"""


import os
import sys
import time
import multiprocessing
import shutil
import fnmatch
import warnings
import subprocess
from datetime import datetime
from osgeo import ogr
from functools import partial
from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive
import pandas as pd

from lthacks import createMetadata

# define function to download tif files found in gDrive folder - called by multiprocessing.Pool().map()
def download_file(fileInfo, outDirPath, nFiles, startTime):
    
    i, fileInfo = fileInfo
    
    fileName = fileInfo['title']
    getFile = drive.CreateFile({'id': fileInfo['id']})

    try:
        getFile.GetContentFile(os.path.join(outDirPath, fileName)) # Download file
    except Exception as e:
        #print('Problem with file %s' % fileName)
        logFile = os.path.join(outDirPath, 'error_log.txt')
        with open(logFile, 'a') as f:
            f.write('%s\n%s\nfile: %s\n\n' % (time.ctime(time.time()), e, fileName))
    
    cumTime = time.time() - startTime
    formatObj = i + 1, nFiles, float(i + 1)/nFiles * 100, cumTime/60
    sys.stdout.write('\rDownloaded file %s of %s (%.1f%%) | cum. time: %.1f mins.' % formatObj)
    sys.stdout.flush()
    


def getGdriveFiles(gDirName, outDirPath, njobs=4):
    # example of inputs
    #gDirName = "ltgee_nbr_ltr_tile28_06200910_prevOneYearOffTest_ltr_stack"
    #outDirPath = "/vol/v1/proj/ltee_vs_ltidl/raster/region_28/ltgee/seg_prevOneYearOff/prep"
    t0 = time.time()
    njobs = int(njobs)
    os.chdir('/vol/v1/general_files/script_library/earth_engine/') #GoogleAuth looks in here for an authorization file - could pass the file as an argument and the get the os.path.dirname

    # authenticate gDrive application and request access to gDrive account
    gauth = GoogleAuth()
    gauth.LocalWebserverAuth() # creates local webserver and auto handles authentication.
    global drive 
    drive = GoogleDrive(gauth)
    
    # find files in the specified gDrive folder
    gDir = drive.ListFile({'q': "mimeType='application/vnd.google-apps.folder' and title contains '"+gDirName+"'"}).GetList()
    if len(gDir) == 1:
        # create the output folder if it does not already exist
        if not os.path.isdir(outDirPath):
            os.mkdir(outDirPath)
        
        # List files in gDrive folder
        query = "'%s' in parents and title contains '.tif'" % gDir[0]['id']
        fileList = pd.DataFrame(drive.ListFile({'q': query}).GetList())
        
        def getIncompleteFiles(fileInfo):
            #fileInfo = pd.DataFrame(fileList)
            fileInfo['localPath'] = fileInfo.title.apply(
                            lambda z: os.path.join(outDirPath, z))
            existing = fileInfo[fileInfo.localPath.apply(
                            lambda z: os.path.exists(z))]
            unfinished = existing[existing.fileSize.astype(int) !=\
                          existing.localPath.apply(lambda z: os.stat(z).st_size)]
            # find files that either didn't finish or that don't exist locally
            files = pd.concat([unfinished, 
                               fileInfo[~fileInfo.title.isin(existing.title)]],
                               ignore_index=True)
            #return [dict(r) for i, r in fileList.iterrows()]#back to what download_files expects
            return files
        
        files = getIncompleteFiles(fileList)
        
        if len(files) == 0:
            sys.exit('No incomplete or un-downloaded files found in %s' % gDirName)
        # wait 10 seconds to start - if the folder is created in the line above
        # then the download won't start, rerunning the script will get it to start
        # could be that the folder is not fully registered before pool.map(func, fileList)
        time.sleep(10)
        
        # continue trying to download files untill all are complete
        nRemaining = len(files)
        print '\n%s files to download...' % nRemaining
        while nRemaining > 0:
            t1 = time.time()
            if njobs > 1:
                pool = multiprocessing.Pool(processes=njobs)
                func = partial(download_file, outDirPath=outDirPath, nFiles=nRemaining, startTime=t1)
                pool.map(func, files.iterrows(), chunksize=1)
                pool.close()
                pool.join()
            else:
                for f in files.iterrows(): download_file(f, outDirPath, nRemaining, t1)
            
            # Check if any files didn't complete
            files = getIncompleteFiles(files)
            nRemaining = len(files)
            njobs = min(njobs, nRemaining)
            
            # Add 2 lines in the error log as separators so it's clear what files
            #   have errors after an additional attempt
            if nRemaining > 0:
                print '\nProblems with %s files. Retrying download...' % nRemaining
                errorLog = os.path.join(outDirPath, 'error_log.txt')
                with open(errorLog, 'a') as f: f.write('\n\t\t##########\t\t' * 2)
            
    else:
        raise RuntimeError('Could not find %s' % gDir)
    
    print '\n\nProcessing time: %.1f minutes' % ((time.time() - t0)/60)


def getInfoFromName(name):
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


def clipAndDecompose(chunkDir, outDir, clipFile,  njobs=0, tileIdField='name', proj='EPSG:5070', returnOutDirs=False):
    
    
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
        info = getInfoFromName(runName)

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
    
    if returnOutDirs:
        return outputDirs
