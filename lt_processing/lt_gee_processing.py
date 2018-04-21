# -*- coding: utf-8 -*-
"""
Function library for LT processing on Google Earth Engine.

Created on Fri Feb 23 11:27:18 2018

@author: shooper
"""


import os
import sys
import time
import glob
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

import ee
ee.Initialize()

from lthacks import createMetadata, attributes_to_df

# #######################################################################################
# ###### ANNUAL SR TIME SERIES STACK BUILDING FUNCTIONS #################################
# #######################################################################################


OAUTH_DIR = '/vol/v1/general_files/script_library/earth_engine/'

def harmonizationRoy(oli):
    ''' Align L8 and l7 
    slope and intercept citation: Roy, D.P., Kovalskyy, V., Zhang, H.K., Vermote, E.F., Yan, L., Kumar, S.S, Egorov, A., 2016, Characterization of Landsat-7 to Landsat-8 reflective wavelength and normalized difference vegetation index continuity, Remote Sensing of Environment, 185, 57-70.(http:#dx.doi.org/10.1016/j.rse.2015.12.024); Table 2 - reduced major axis (RMA) regression coefficients
    '''
    slopes = ee.Image.constant([0.9785, 0.9542, 0.9825, 1.0073, 1.0171, 0.9949]); # create an image of slopes per band for L8 TO L7 regression line - David Roy
    itcp = ee.Image.constant([-0.0095, -0.0016, -0.0022, -0.0021, -0.0030, 0.0029]); # create an image of y-intercepts per band for L8 TO L7 regression line - David Roy
    # select OLI bands 2-7 and rename them to match L7 band names
    # ...resample the L8 bands using bicubic
    # ...multiply the y-intercept bands by 10000 to match the scale of the L7 bands then apply the line equation - subtract the intercept and divide by the slope
    # ...set the output system:time_start metadata to the input image time_start otherwise it is null
    y = oli.select(['B2','B3','B4','B5','B6','B7'],['B1', 'B2', 'B3', 'B4', 'B5', 'B7'])\
            .resample('bicubic')\
            .subtract(itcp.multiply(10000)).divide(slopes)\
            .set('system:time_start', oli.get('system:time_start')); 
    
    return y.toShort(); # set image to signed 16-bit integer 




# ------ DEFINE FUNCTION TO RETRIEVE A SENSOR SR COLLECTION -----------------------------
 
def getSRcollection(year, startDay, endDay, sensor, box):
    
    # get surface reflectance images
    # ...filter them by a bounding box
    # ...filter their dates from June 1st - Sept 30th
    srCollection = ee.ImageCollection('LANDSAT/'+ sensor + '/C01/T1_SR')\
                     .filterBounds(box)\
                     .filterDate(str(year)+'-'+startDay, str(year)+'-'+endDay); 
    
    def function(img):
        dat = ee.Image(
                ee.Algorithms.If(
                    sensor == 'LC08', # condition - if image is OLI
                    harmonizationRoy(img.unmask()), # true - then apply the L8 TO L7 alignment function after unmasking pixels that were previosuly masked (why/when are pixels masked)
                    # false - else select out the reflectance bands from the non-OLI image and unmask any previous pixels
                    # ...resample by bicubic 
                    # ...set the output system:time_start metadata to the input image time_start otherwise it is null
                    img.select(['B1', 'B2', 'B3', 'B4', 'B5', 'B7']).unmask()\
                    .resample('bicubic')\
                    .set('system:time_start', img.get('system:time_start'))
                    ));
        qa = img.select('pixel_qa')
        #Shadow
        # ..snow
        # ..cloud
        mask = qa.bitwiseAnd(8).eq(0)\
                .And(qa.bitwiseAnd(16).eq(0))\
                .And(qa.bitwiseAnd(32).eq(0))
        datMasked = dat.mask(mask)        
        cloudMask = img.select('pixel_qa').bitCount().lte(2); # select out the fmask layer and create a 0/1 mask - set 0,1 to 1, all else to 0; 0=clear; 1=water; 2=shadow; 3=snow; 4=cloud
        datMasked = dat.mask(cloudMask); #apply the mask - 0's in mask will be excluded from computation and set to opacity=0 in display
        return datMasked; # return     
         
    srCollection = srCollection.map(function);
    
    return srCollection; # return 



# ------ DEFINE FUNCTION TO COMBINE LT5, LE7, & LC8 COLLECTIONS -------------------------

def getCombinedSRcollection(year, startDay, endDay, box):
    lt5 = getSRcollection(year, startDay, endDay, 'LT05', box); # get TM collection for a given year and bounding area
    le7 = getSRcollection(year, startDay, endDay, 'LE07', box); # get ETM+ collection for a given year and bounding area
    lc8 = getSRcollection(year, startDay, endDay, 'LC08', box); # get OLI collection for a given year and bounding area
    
    return  ee.ImageCollection(lt5.merge(le7).merge(lc8)); # merge the individual sensor collections into one imageCollection object



# ------ DEFINE FUNCTION TO REDUCE COLLECTION TO SINGLE IMAGE BY MEDOID -----------------

# medoid composite with equal weight among indices
# Medoids are representative objects of a data set or a cluster with a data set whose average dissimilarity to all the objects in the cluster is minimal. Medoids are similar in concept to means or centroids, but medoids are always members of the data set.
def medoidMosaic(inCollection, dummyCollection):
    
    imageCount = inCollection.toList(1).length();
    finalCollection = ee.ImageCollection(ee.Algorithms.If(imageCount.gt(0), inCollection, dummyCollection));
    median = ee.ImageCollection(finalCollection).median(); # calculate the median of the annual image collection - returns a single 6 band image - the collection median per band
    
    def function(img):
        diff = ee.Image(img).subtract(median).pow(ee.Image.constant(2)); # get the difference between each image/band and the corresponding band median and take to power of 2 to make negatives positive and make greater differences weight more
        return diff.reduce('sum').addBands(img); #per image in collection, sum the powered difference across the bands - set this as the first band add the SR bands to it - now a 7 band image collection
    medoid = finalCollection.map(function);
    return ee.ImageCollection(medoid).reduce(ee.Reducer.min(7)).select([1,2,3,4,5,6], ['B1','B2','B3','B4','B5','B7']); # find the powered difference that is the least - what image object is the closest to the median of teh collection - and then subset the SR bands and name them - leave behind the powered difference band



# ------ DEFINE FUNCTION TO REDUCE COLLECTION TO SINGLE IMAGE BY DISTANCE FROM MEDIAN DAY ---
def targetDayMoasic(inCollection, targetDay):
    
    def function(image):
        day = ee.Date(image.get('system:time_start')).getRelative('day', 'year');
        delta = image.select(None).addBands(day.subtract(targetDay).abs().multiply(-1)).int16();
        return delta.select([0], ['delta']).addBands(image);

    inCollectionDelta = inCollection.map(function)
    
    return inCollectionDelta.qualityMosaic('delta')\
                          .select([1,2,3,4,5,6]);




# ------ DEFINE FUNCTION TO APPLY A MOSAIC FUNCTION TO A COLLECTION -------------------------------------------

def buildMosaic(year, startDay, endDay, box, mosaicType, targetDay, dummyCollection):
    
    tmp = [] # create a temp variable to hold the upcoming annual mosiac
    collection = getCombinedSRcollection(year, startDay, endDay, box); # get the SR collection
    if(mosaicType == "medoid"):
        tmp = medoidMosaic(collection, dummyCollection); # reduce the collection to single image per year by medoid 
    elif (mosaicType == "targetDay"):
        tmp = targetDayMoasic(collection, targetDay); # reduce the collection to single image per year by medoid
    img = tmp.set('system:time_start', time.mktime((year, 8, 1, 0,0,0,0,0,0))) #(new Date(year,8,1)).valueOf()); # add the year to each medoid image
    
    return ee.Image(img); # return as image object



# ------ DEFINE FUNCTION TO BUILD ANNUAL MOSAIC COLLECTION ------------------------------
def buildMosaicCollection(startYear, endYear, startDay, endDay, box, mosaicType, targetDay, dummyCollection):
    imgs = []; #create empty array to fill
    for i in range(startYear, endYear + 1): # for each year from hard defined start to end build medoid composite and then add to empty img array
        tmp = buildMosaic(i, startDay, endDay, box, mosaicType, targetDay, dummyCollection); # build the medoid mosaic for a given year
        imgs.append(tmp.set('system:time_start', time.mktime((i, 8, 1, 0,0,0,0,0,0)))); # concatenate the annual image medoid to the collection (img) and set the date of the image - hardwired to the year that is being worked on for Aug 1st
    
    return ee.ImageCollection(imgs); #return the array img array as an image collection

# #######################################################################################
# #######################################################################################
# #######################################################################################






# #######################################################################################
# ###### INDEX CALCULATION FUNCTIONS ####################################################
# #######################################################################################

# TASSELLED CAP
def tcTransform(img):
    b = ee.Image(img).select(["B1", "B2", "B3", "B4", "B5", "B7"]); # select the image bands
    brt_coeffs = ee.Image.constant([0.2043, 0.4158, 0.5524, 0.5741, 0.3124, 0.2303]); # set brt coeffs - make an image object from a list of values - each of list element represents a band
    grn_coeffs = ee.Image.constant([-0.1603, -0.2819, -0.4934, 0.7940, -0.0002, -0.1446]); # set grn coeffs - make an image object from a list of values - each of list element represents a band
    wet_coeffs = ee.Image.constant([0.0315, 0.2021, 0.3102, 0.1594, -0.6806, -0.6109]); # set wet coeffs - make an image object from a list of values - each of list element represents a band
  
    sum = ee.Reducer.sum(); # create a sum reducer to be applyed in the next steps of summing the TC-coef-weighted bands
    brightness = b.multiply(brt_coeffs).reduce(sum); # multiply the image bands by the brt coef and then sum the bands
    greenness = b.multiply(grn_coeffs).reduce(sum); # multiply the image bands by the grn coef and then sum the bands
    wetness = b.multiply(wet_coeffs).reduce(sum); # multiply the image bands by the wet coef and then sum the bands
    # stack TCG and TCW behind TCB with .addBands, use select() to name the bands
    tc = brightness.addBands(greenness)\
                     .addBands(wetness)\
                     .select([0,1,2], ['TCB','TCG','TCW'])\
                     .set('system:time_start', img.get('system:time_start'));
    return tc;

# NBR
def nbrTransform(img):
    # calculate normalized difference of B4 and B7. orig was flipped: ['B7', 'B4']
    #... scale results by 1000
    #... name the band
    nbr = img.normalizedDifference(['B4', 'B7'])\
                 .multiply(1000)\
                 .select([0], ['NBR'])\
                 .set('system:time_start', img.get('system:time_start'));
    return nbr;

# NDVI
def ndviTransform(img):
    # calculate normalized dif between band 4 and band 3 (B4-B3/B4_B3)
    # ... scale results by 1000
    # ... name the band
    ndvi = img.normalizedDifference(['B4', 'B3'])\
                .multiply(1000)\
                .select([0], ['NDVI'])\
                .set('system:time_start', img.get('system:time_start'));
    return ndvi;
                
# NDSI
def ndsiTransform(img):
    # calculate normalized dif between band 4 and band 3 (B4-B3/B4_B3)
    #... scale results by 1000
    #... name the band
    ndsi = img.normalizedDifference(['B2', 'B5'])\
                .multiply(1000)\
                .select([0], ['NDSI'])\
                .set('system:time_start', img.get('system:time_start'));
    return ndsi;


# MAKE AN LT STACK
def makeLtStack(img):
    
    
    tc = tcTransform(img);
    if index == 'NBR':    indexImg = nbrTransform(img).multiply(indexFlip)
    elif index == 'B5':   indexImg = img.select(['B5']).multiply(indexFlip)
    elif index == 'NDVI': indexImg = ndviTransform(img).multiply(indexFlip)
    elif index == 'TCB':  indexImg = tc.select(['TCB']).multiply(indexFlip)
    elif index == 'NDSI': indexImg = ndsiTransform(img).multiply(indexFlip)
    elif index == 'TCG':  indexImg = tc.select(['TCG']).multiply(indexFlip)
    elif index == 'B3':   indexImg = img.select(['B3']).multiply(indexFlip)
    else:
        raise RuntimeError('The index you provided is not supported: %s' % index)
    
    tc = tc.select(['TCB', 'TCG', 'TCW'],['FTV_TCB', 'FTV_TCG', 'FTV_TCW']);
    ndvi = ndviTransform(img).select(['NDVI'],['FTV_NDVI']);
    ndsi = ndsiTransform(img).select(['NDSI'],['FTV_NDSI']);

    indexImgFtv = indexImg.select([0], ['FTV_IDX']).multiply(indexFlip);
    allStack = indexImg.addBands(indexImgFtv)\
                        .addBands(tc)\
                        .addBands(ndvi)\
                        .addBands(ndsi)\
                        .set('system:time_start', img.get('system:time_start'));
    
    return allStack;
        

# #######################################################################################
# ###### LANDTRENDR #####################################################################
# #######################################################################################

# ------ DEFINE FUNCTION TO EXTRACT VERTICES FROM LT RESULTS AND STACK BANDS ------------

def getLTvertStack(LTresult, run_params):
    
    emptyArray = [];
    vertLabels = [];
    for i in range(1, run_params['maxSegments']+2):
        vertLabels.append("Vert" + str(i));
        emptyArray.append(0);
    zeros = ee.Image(ee.Array([emptyArray,
                               emptyArray,
                               emptyArray]));
    lbls = [['year','raw','ftv'], vertLabels] #['Vert1','Vert2','Vert3','Vert4','Vert5','Vert6','Vert7']]; # labels for 2 dimensions of the array that will be cast to eachother in the final step of creating the vertice output
    
    vmask = LTresult.arraySlice(0,3,4); # slices out the 4th row of a 4 row x N col (N = number of years in annual stack) matrix, which identifies vertices - contains only 0s and 1s, where 1 is a vertex (referring to spectral-temporal segmentation) year and 0 is not . 
    
    # uses the sliced out isVert row as a mask to only include vertice in this data - after this a pixel will only contain as many "bands" are there are vertices for that pixel - min of 2 to max of 7.
    #... from the vertOnly data subset slice out the vert year row, raw spectral row, and fitted spectral row
    #... adds the 3 row x 7 col 'zeros' matrix as a band to the vertOnly array - this is an intermediate step to the goal of filling in the vertOnly data so that there are 7 vertice slots represented in the data - right now there is a mix of lengths from 2 to 7
    #... concatenates the 3 row x 7 col 'zeros' matrix band to the vertOnly data so that there are at least 7 vertice slots represented - in most cases there are now > 7 slots filled but those will be truncated in the next step
    #... 7 # before this line runs the array has 3 rows and between 9 and 14 cols depending on how many vertices were found during segmentation for a given pixel. this step truncates the cols at 7 (the max verts allowed) so we are left with a 3 row X 7 col array
    #... this takes the 2-d array and makes it 1-d by stacking the unique sets of rows and cols into bands. there will be 7 bands (vertices) for vertYear, followed by 7 bands (vertices) for rawVert, followed by 7 bands (vertices) for fittedVert, according to the 'lbls' list
    #import pdb; pdb.set_trace()
    #ltVertStack = LTresult.arrayMask(vmask).arraySlice(0, 0, 3).addBands(zeros).toArray(1).arraySlice(1, 0, run_params['maxSegments']+1).arrayFlatten(lbls, '')#'''
                      #.toShort();'''    
    ltVertStack = LTresult.arrayMask(vmask)\
                          .arraySlice(0, 0, 3)\
                          .addBands(zeros)\
                          .toArray(1)\
                          .arraySlice(1, 0, run_params['maxSegments']+1)\
                          .arrayFlatten(lbls, '')
                          #.toShort();'''
    return ltVertStack;


def run_lt(featureCol, featureKey, featureValue, aoiBuffer, indexDict, startYear, endYear, run_params, startDay, endDay, mosaicType, targetDay, gDriveFolder, outProj, affine, trackTasks=True):
    t0 = time.time()
    # get geometry stuff
    # # load the aoi
    #... filter by park - left arg is the key name, right arg is the value to match
    aoi = ee.FeatureCollection('ft:' + featureCol)\
            .filter(ee.Filter.stringContains(featureKey, featureValue))\
            .geometry()\
            .buffer(aoiBuffer);
    
    box = aoi.bounds(); # get the bounds from the drawn polygon
    
    
    # make a dummy collection for filling missing years (if there are any)
    dummyCollection = ee.ImageCollection([ee.Image([0,0,0,0,0,0]).mask(ee.Image(0))]); # make a dummy collection to fill in potentially missing years
    #indexList = [['NBR', -1], ['B5', 1], ['NDVI', -1], ['TCB', 1], ['NDSI', -1], ['TCG', -1], ['B3', 1]];
    
    years = ['yr%s' % year for year in range(startYear, endYear + 1)]
    
    vertYearLabels = [];
    rawVertLabels = [];
    ftvVertLabels = [];
    for i in range(1, run_params['maxSegments'] + 2): # should go up to maxSeg. + 1
        vertYearLabels.append("yearVert" + str(i));
        rawVertLabels.append("rawVert" + str(i));
        ftvVertLabels.append("ftvVert" + str(i));
        
    # make sure makeLtStack can use these. Since imageCollection.map()
    #   can't accept additional args, just set them as globals.
    global index 
    global indexFlip
    
    tasks = []
    for i in [0]:#range(len(startDay)):
        # build the annual SR collection
        annualSRcollection = buildMosaicCollection(startYear, endYear, startDay[i], endDay[i], box, mosaicType, targetDay, dummyCollection); # put together the cloud-free medoid surface reflectance annual time series collection
        #for index, indexFlip in indexDict.iteritems():
        for index, indexFlip in [['NBR', 1]]:
            #get the index and the index flipper
            #index, indexFlip = indexList[j]; # pull out the index to segment on and value to flip index or not
            
            # make the collection for this index and add the collection to the run parameters
            ltCollection = annualSRcollection.map(makeLtStack); # make the LT collection for this run
            run_params['timeSeries'] = ltCollection; # add the single spectral index annual time series collection to the segmentation run parameter object
            # run LT
            lt = ee.Algorithms.Test.LandTrendr(**run_params); # run LandTrendr spectral temporal segmentation algorithm
            
            # pull out the LT vert layers
            ltVertStack = getLTvertStack(lt.select("LandTrendr"), run_params); # extract the year, raw spectral value, and fitted values for vertices as a stacked bands

            # rearrange the vert layers
            #('yearVert1','yearVert2','yearVert3','yearVert4','yearVert5','yearVert6','yearVert7');
            #... ('rawVert1','rawVert2','rawVert3','rawVert4','rawVert5','rawVert6','rawVert7')
            #...#('ftvVert1','ftvVert2','ftvVert3','ftvVert4','ftvVert5','ftvVert6','ftvVert7')
            vertYear = ltVertStack.select(vertYearLabels); 
            vertRaw = ltVertStack.select(rawVertLabels).multiply(indexFlip);
            vertFtv = ltVertStack.select(ftvVertLabels).multiply(indexFlip);
            # pull out the FTV layers
            ftvIDX  = lt.select(['FTV_IDX_fit' ]).arrayFlatten([years]);
            ftvTCB  = lt.select(['FTV_TCB_fit' ]).arrayFlatten([years]);
            ftvTCG  = lt.select(['FTV_TCG_fit' ]).arrayFlatten([years]);
            ftvTCW  = lt.select(['FTV_TCW_fit' ]).arrayFlatten([years]);
            ftvNDVI = lt.select(['FTV_NDVI_fit']).arrayFlatten([years]);
            ftvNDSI = lt.select(['FTV_NDSI_fit']).arrayFlatten([years]);
          
            # stack all the layers up
            allStack = vertYear.addBands(vertRaw)\
                               .addBands(vertFtv)\
                               .addBands(ftvIDX)\
                               .addBands(ftvTCB)\
                               .addBands(ftvTCG)\
                               .addBands(ftvTCW)\
                               .addBands(ftvNDVI)\
                               .addBands(ftvNDSI)\
                               .round()\
                               .toShort();
        
            # make a file name
            nVert = int(run_params['maxSegments']) + 1;
            fNameDict = {'featureKey': featureKey,
                         'featureValue': featureValue,
                         'index': index,
                         'nVert': nVert,
                         'startYear': startYear,
                         'endYear': endYear,
                         'startDay': startDay[i].replace('-',''),
                         'endDay': endDay[i].replace('-','')
                         }
            fileNamePrefix = ('{featureKey}-{featureValue}-{index}-{nVert}-'+\
                              '{startYear}{endYear}-{startDay}{endDay}')\
                              .format(**fNameDict)

            task = ee.batch.Export.image.toDrive(allStack.clip(aoi),
                                         description=fileNamePrefix, 
                                         folder=gDriveFolder,
                                         fileNamePrefix=fileNamePrefix,
                                         crs=outProj,
                                         crsTransform=affine,
                                         region=aoi.bounds().getInfo()['coordinates'],
                                         maxPixels=1e13
                                         )#'''
            task.start()
            print task.status()['id']
            tasks.append(task)
    
    def getStatus():
        allTasks = pd.DataFrame(ee.data.getTaskList())
        
        # Create a dataframe index based on the task creation time.
        taskStatus = allTasks.copy()
        taskStatus.index = pd.to_datetime(taskStatus['creation_timestamp_ms'], unit='ms')
        del(taskStatus['creation_timestamp_ms'])
        
        # Convert the start and update timestamps to Python data types.
        taskStatus.start_timestamp_ms = pd.to_datetime(taskStatus.start_timestamp_ms, unit='ms')
        taskStatus = taskStatus.loc[taskStatus.index > pd.to_datetime(t0, unit='s')]
        return taskStatus
    taskStatus = getStatus()
    task_ids = taskStatus['id']#'''
    
    if trackTasks:
        taskStatus = pd.DataFrame([t.status() for t in tasks])
        taskStatus['active'] = True
        nTasks = float(len(tasks))
        SLEEP = .5
        t1 = time.time()
        # Wait until all tasks have finished running 
        #while taskStatus['active'].any():
        while taskStatus.state.isin(['READY', 'RUNNING']).any():
            time.sleep(SLEEP)
            #inactive = taskStatus[~taskStatus.active]
            taskStatus = getStatus()
            inactive = taskStatus[~taskStatus.state.isin(['READY', 'RUNNING'])]
            failed = inactive[inactive.state.isin(['FAILED', 'CANCELED', 'CANCEL_REQUESTED'])]
            complete = inactive[inactive.state == 'COMPLETE']
            nFinished = len(complete) + len(failed)
            cumTime = time.time() - t1
            msg = '\r%s of %d (%.1f%%) tasks done | cum. time: %.1f mins' % (nFinished, nTasks, (nFinished/nTasks) * 100, cumTime/60)
            sys.stdout.write(msg)
            sys.stdout.flush()
            '''taskStatus = pd.DataFrame([t.status() for t in tasks])#retrieve new status
            taskStatus['active'] = [t.active() for t in tasks] #update active column'''
            
        if len(failed) > 0:
            print '\n\nIDs of failed tasks:\n\t' + '\n\t'.join(failed['id'].tolist())#'''

    return tasks, gDriveFolder


def authenticateGDrive():
    
    workingDir = os.getcwd()
    os.chdir(OAUTH_DIR) # oAuth looks for client secrets file here
    
    # authenticate gDrive application and request access to gDrive account
    gauth = GoogleAuth()
    gauth.LocalWebserverAuth() # creates local webserver and auto handles authentication.
    gDrive = GoogleDrive(gauth)
    
    os.chdir(workingDir) # Switch back to workingDir
    
    return gDrive

def getStatus(startTime, tasks=None, taskFilter='*', columns=None):
    
    if tasks is None:
        allTasks = pd.DataFrame(ee.data.getTaskList())
    else:
        allTasks = pd.DataFrame([t.status() for t in tasks])
    
    if '!' in taskFilter:
        pattern, notpattern = taskFilter.split('!')
        matches = fnmatch.filter(allTasks.description, pattern)
        notmatches = fnmatch.filter(allTasks.description, notpattern)
        matches = [m for m in matches if m not in notmatches]
    else:
        matches = fnmatch.filter(allTasks.description, taskFilter)
    
    # Create a dataframe index based on the task creation time.
    taskStatus = allTasks[allTasks.description.isin(matches)].copy()
    taskStatus.index = pd.to_datetime(taskStatus['creation_timestamp_ms'], unit='ms')
    del(taskStatus['creation_timestamp_ms'])
    
    # Convert the start and update timestamps to Python data types.
    taskStatus.start_timestamp_ms = pd.to_datetime(taskStatus.start_timestamp_ms, unit='ms')
    taskStatus = taskStatus.loc[taskStatus.index > pd.to_datetime(startTime, unit='s')]
    
    return taskStatus
    
    
def _downloadFromQueue(inQueue, downloadList):#, out_queue):

    while True:
        args = inQueue.get() #Blocks until something is ready
        if args is None:
            print '\n#######None########\n'
            break
        return_args = download_file(*args)
        if return_args is not None: #Download failed
            inQueue.put(args) # Put it back in the queue and try again
        
        # Let out_queue know the file has been processed
        else:
            outFile = args[0][1]['title']#first arg, 1th item is fileInfo series
            downloadList.append(outFile)#'''


def _callTranslateFromQueue(inQueue):#, out_queue):

    while True:
        cmd = inQueue.get() #Blocks until something is ready
        if cmd is None:
            print '\n#########None############\n'
            break
        callTranslate(cmd)
            
            
def listenAndDownladTasks(tasks, downloadDir, gDriveFolder, gDrive=None, outDir=None, clipFile=None, sleep=.5, njobs=10, silent=True, timestamp_offset=30, logFile=None):
    
    # check that clipFile exists if outDir was specified
    if outDir:
        try:
            clipFileCheck = ogr.Open(clipFile)
        except:
             raise RuntimeError('clipFile is not a valid OGR readable file: %s' % clipFile)
        if not clipFileCheck:
            raise RuntimeError('clipFile is not a valid OGR readable file: %s' % clipFile)
        clipFileCheck = None
        if not os.path.isdir(outDir):
            os.mkdir(outDir)
        
    # If gdrive hasn't been authenticated yet, authenticate
    if gDrive == None:
        gDrive = authenticateGDrive()

    if not os.path.isdir(downloadDir):
        os.mkdir(downloadDir)
    
    if logFile:
        logFile = os.path.join(downloadDir, 'download_log.csv')
    
    mgr = multiprocessing.Manager()
    downloadQueue = mgr.Queue()
    downloadedList = mgr.list()
    pool = multiprocessing.Pool(njobs)
    # prime the pool to process the queues
    pool.apply_async(_downloadFromQueue, (downloadQueue, downloadedList))
    if outDir:
        decomposeQueue = mgr.Queue()
        pool.apply_async(_callTranslateFromQueue, (decomposeQueue, ))
    
    # Get task progress info
    t0 = time.time() - timestamp_offset * 60
    
    taskStatus = getStatus(t0)#(tasks, taskFilter='*!*info')
    
    #taskStatus = pd.DataFrame([t.status() for t in tasks])
    taskStatus['active'] = True
    nTasks = float(len(taskStatus))
    #nTasks = 7
    t1 = time.time()
    # Wait until all tasks have finished running 
    #while taskStatus['active'].any():
    taskStatus['downloading'] = False
    taskStatus['downloadDone'] = False
    downloadDict = {}#pd.DataFrame(columns=['taskName', '])
    
    FAILED_STATES = ['FAILED', 'CANCELLED', 'CANCEL_REQUESTED']
    while not taskStatus.downloadDone.all():
        #inactive = taskStatus[~taskStatus.active]
        '''# Update status
        taskStatus.state = pd.DataFrame([t.status() for t in tasks])['state']#update state
        taskStatus['active'] = [t.active() for t in tasks] #update active column'''
        taskStatus['state'] = getStatus(t0)['state']
        
        # check if all tasks are finished
        #if not (taskStatus.state == 'COMPLETED').all():
        failed = taskStatus.loc[taskStatus.state.isin(FAILED_STATES)]
        complete = taskStatus.loc[taskStatus.state == 'COMPLETED']
        downloadReady = complete.loc[~complete.downloading]# EE complete, but not downloading yet
        nFinished = len(complete) + len(failed)
        taskStatus.loc[failed.index, 'downloadDone'] = True #Mark as done

        # Check if download hasn't yet started for any complete files 
        if len(downloadReady) > 0:
            searchStrs = [task.description + '*' for i, task in downloadReady.iterrows()]
            readyInfo = listDriveFiles(gDrive, gDriveFolder, searchPatterns=searchStrs)
            for fileInfo in readyInfo.iterrows():
                downloadQueue.put([fileInfo, gDrive, downloadDir, 0])
            # Set download status to True
            taskStatus.loc[downloadReady.index, 'downloading'] = True
            for s in searchStrs:
                downloadDict[s[:-1]] = fnmatch.filter(readyInfo.title, s)
        # If all tasks are done, break the while loop in _downloadFromQueue
        if nFinished == nTasks:
            downloadQueue.put(None)
            taskStatus.downloading = True
        
        # Check to see if all files from any task are done
        completeImgs = complete.loc[complete.description.apply(lambda x: not x.endswith('info')) & ~complete.downloadDone]
        if outDir is not None: # If it's None, don't decompose
            for i, taskName in completeImgs.description.iteritems():
                # Check that text files have been downloaded
                # First check that the tasks are complete
                csvState = getStatus(t0, taskFilter=taskName+'*info').state
                if (csvState == 'COMPLETED').all():
                    csvInfo = listDriveFiles(gDrive, gDriveFolder, searchPatterns=taskName+'*info*')
                    
                    # If the files aren't completely downloaded yet, try the next task
                    if len(getIncompleteFiles(csvInfo, outDir, sizeDifTolerance=500)) > 0:
                        continue
                    runCsv =  csvInfo[csvInfo.title.apply(lambda x: x.endswith('run_info.csv'))].iloc[0]
                    bandCsv = csvInfo[csvInfo.title.apply(lambda x: x.endswith('band_info.csv'))].iloc[0]
                    
                    # Check that the images have all completed downlaoding
                    downloaded = sorted(fnmatch.filter(downloadedList, taskName + '*.tif'))#assumes gee always exports tifs
                    fromGDrive = sorted(fnmatch.filter(downloadDict[taskName], '*.tif'))
                    if  downloaded == fromGDrive:#all([True else False for name in taskName if name in downloadedList]):
                        #fullPaths = [os.path.join(downloadDir, f) for f in downloaded]
                        runPath = os.path.join(downloadDir, runCsv.title)
                        bandPath = os.path.join(downloadDir, bandCsv.title)
                        commands = getDecomposeCommands(downloadDir, runPath, bandPath, outDir, clipFile)
                        for cmd in commands:
                            decomposeQueue.put(cmd)
                        #import pdb; pdb.set_trace()
                        taskStatus.loc[i, 'downloadDone'] = True
                # If either of the CSVs or the image failed, mark the download as done
                elif (csvState.isin(FAILED_STATES)).any() or taskStatus.loc[i, 'state'].isin(FAILED_STATES):
                    taskStatus.loc[i, 'downloadDone'] = True
                else:
                    continue
            if taskStatus.downloadDone.all(): #Send termination signal to queue
                decomposeQueue.put(None)
        # If outDir is None, 
        else:
            for i, taskName in completeImgs.description.iteritems():
                if len(listDriveFiles(gDrive, gDriveFolder, downloadDir, searchPatterns=['%s*.tif' % taskName], sizeDifTolerance=500)) == 0:
                    taskStatus.loc[i, 'downloadDone'] = True
                    
            # check if any downloads are done sure
            if len(listDriveFiles(gDrive, gDriveFolder, downloadDir, sizeDifTolerance=500)) == 0 and taskStatus.downloading.all():
                taskStatus.downloadDone = True #this will break out of for loop
        #print 'downloads done:', len(taskStatus[taskStatus.downloadDone])
        if logFile:
            taskStatus.to_csv(logFile)
        
        if not silent:
            cumTime = time.time() - t1
            msg = '\r%s of %d (%.1f%%) tasks done | cum. time: %.1f mins' % (nFinished, nTasks, (nFinished/nTasks) * 100, cumTime/60)
            sys.stdout.write(msg)
            sys.stdout.flush()
            
        '''if (taskStatus.downloadDone | taskStatus.index.isin(failed.index)).all():
            break'''
    
    # Close all processes
    for _ in range(njobs): 
        downloadQueue.put(None) # add kill switch to break out of while loop
    pool.close()
    pool.join()
        
    if len(failed) > 0:
        print '\n\nIDs of failed tasks:\n\t' + '\n\t'.join(failed['id'].tolist())


# define function to download tif files found in gDrive folder - called by multiprocessing.Pool().map()
def download_file(fileInfo, gDrive, outDirPath, startTime, nFiles=None):
    
    i, fileInfo = fileInfo
    
    fileName = fileInfo['title']
    getFile = gDrive.CreateFile({'id': fileInfo['id']})
    
    try:
        getFile.GetContentFile(os.path.join(outDirPath, fileName)) # Download file
    except Exception as e:
        #print('Problem with file %s' % fileName)
        logFile = os.path.join(outDirPath, 'error_log.txt')
        with open(logFile, 'a') as f:
            f.write('%s\n%s\nfile: %s\n\n' % (time.ctime(time.time()), e, fileName))
        return [fileInfo, gDrive, outDirPath, startTime, nFiles]
    
    if nFiles:
        cumTime = time.time() - startTime
        formatObj = i + 1, nFiles, float(i + 1)/nFiles * 100, cumTime/60
        sys.stdout.write('\rDownloaded file %s of %s (%.1f%%) | cum. time: %.1f mins.' % formatObj)
        sys.stdout.flush()

  
def getIncompleteFiles(fileInfo, outDirPath, sizeDifTolerance=0):
    '''Find any files that have not been completely downloaded'''
    
    fileInfo['localPath'] = fileInfo.title.apply(
                    lambda z: os.path.join(outDirPath, z))
    existing = fileInfo.loc[fileInfo.localPath.apply(
                    lambda z: os.path.exists(z))].copy()
    existing['sizeDif'] = existing.fileSize.astype(int) - existing.localPath.apply(lambda z: os.stat(z).st_size)
    unfinished = existing.loc[existing.sizeDif >= sizeDifTolerance]

    # find files that either didn't finish or that don't exist locally
    fileInfo = pd.concat([unfinished, 
                       fileInfo[~fileInfo.title.isin(existing.title)]],
                       ignore_index=True)

    return fileInfo

  
def listDriveFiles(gDrive, gDirName, outDirPath=None, searchPatterns=None, sizeDifTolerance=0):
    '''List files in Google Drive folder `gDriveFolder`. If `outDirPath` 
    is given, return only files that have not been completely downloaded. If
    `searchPatterns` is given, return only files that match any of patterns'''
    
    # List files in gDrive folder
    gDriveFolder = gDrive.ListFile({'q': "mimeType='application/vnd.google-apps.folder' and title contains '"+gDirName+"'"}).GetList()
    query = "'%s' in parents" % gDriveFolder[0]['id']
    fileInfo = pd.DataFrame(gDrive.ListFile({'q': query}).GetList())

    if outDirPath:
        fileInfo = getIncompleteFiles(fileInfo, outDirPath, sizeDifTolerance=sizeDifTolerance)
    
    # Filter out files that don't match one of the searchPatterns
    if searchPatterns:
        # If a single pattern is given, make it a list
        if not hasattr(searchPatterns, '__iter__'):
            searchPatterns = [searchPatterns]
        matches = []
        for pattern in searchPatterns:
            matches.extend(fnmatch.filter(fileInfo.title, pattern))
        fileInfo = fileInfo[fileInfo.title.isin(matches)]
    
    return fileInfo
        

def getGdriveFiles(gDirName, outDirPath, njobs=4, gDrive=None, searchPatterns=None):
    t0 = time.time()
    njobs = int(njobs)
    outDirPath = os.path.abspath(outDirPath)
    
    # If gdrive hasn't been authenticated yet, authenticate
    if gDrive == None:
        gDrive = authenticateGDrive()
    
    # find files in the specified gDrive folder
    gDir = gDrive.ListFile({'q': "mimeType='application/vnd.google-apps.folder' and title contains '"+gDirName+"'"}).GetList()
    if len(gDir) == 1:
        # create the output folder if it does not already exist
        if not os.path.isdir(outDirPath):
            os.mkdir(outDirPath)
        
        # List files in gDrive folder        
        files = listDriveFiles(gDrive, gDir, outDirPath, searchPatterns)
        
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
                func = partial(download_file, gDrive=gDrive, outDirPath=outDirPath, startTime=t1, nFiles=nRemaining)
                pool.map(func, files.iterrows(), chunksize=1)
                pool.close()
                pool.join()
            else:
                for f in files.iterrows(): download_file(f, outDirPath, t1, nRemaining)
            
            # Check if any files didn't complete
            files = getIncompleteFiles(files, outDirPath)
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
    

def callTranslate(callInfo, silent=True):
    outFile, command = callInfo[:2]
    if len(callInfo) > 2:
        outFile, command, count, nCommands, startTime = callInfo
        silent = False
        cumTime = time.time() - startTime
        formatObj = count, nCommands, float(count)/nCommands * 100, cumTime/60
    if not silent:
        sys.stdout.write('\rProcessing %s of %s files (%.1f%%) | cum. time: %.1f mins.' % formatObj)
        sys.stdout.flush()
    subprocess.call(command, shell=True)
    desc = 'LandTrendr data created with Google Earth Engine. Filename is in the form {version}_{index}_{dateRange}_{tileID}_{nVertices}_{processingDate}_{LTdataType}_{LTdataContent}.tif'
    createMetadata(sys.argv, outFile, description=desc)


'''def getDecomposeCommands(eeTaskName, chunkDir, pieces, outImgTypes, outDir, clipFile, tileIdField='name', proj='EPSG:5070', nTiles=9, count=None, nCommands=None, startTime=None):
    
    
    # get info about the GEE run
    info = getInfoFromName(eeTaskName)
    import pdb; pdb.set_trace()
    
    # make a list of tile tifs
    vrtFile = os.path.join(chunkDir, eeTaskName + '.vrt')
    tileListFile = vrtFile.replace('.vrt', '_filelist.txt')
    with open(tileListFile, 'w') as tileList:
        for piece in pieces:
            tileList.write(piece + '\n')

    # create vrt
    cmd = 'gdalbuildvrt -input_file_list %s %s' % (tileListFile, vrtFile)
    subprocess.call(cmd, shell=True)

    # make a list of band ranges for each out type
    vertStops = []
    for vertType in range(4):  
        vertStops.append(vertType*info['nVert']+1)

    nYears = (info['endYear'] - info['startYear']) + 1
    ftvStops = []  
    for ftvType in range(1, len(outImgTypes)-2):  
        ftvStops.append(ftvType*nYears+vertStops[-1])

    bandStops = vertStops+ftvStops
    bandRanges = [range(bandStops[i],bandStops[i+1]) for i in range(len(bandStops)-1)]

    
    # Check that the directory structure exists. If not, build it.
    indexDir = os.path.join(outDir, '%s_%s' % (info['indexID'], info['dateRange']))
    if not os.path.isdir(indexDir):
        os.mkdir(indexDir)
    
    # Open the tile layer and loop through each tile within this region
    #   Use attribute value instead of geographic bounds of vrt because GEE
    #   often includes nodata buffer of unpredctable size
    dataset = ogr.Open(clipFile)
    layer = dataset.GetLayer()
    filterString = "%s = '%s'" % (info['key'], info['value'])
    layer.SetAttributeFilter(filterString)
    feature = layer.GetNextFeature()
    i = 1
    
    commands = []
    c = count
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
        for i in range(len(outImgTypes)):   
            #outBname = eeTaskName+'-'+outImgTypes[i]
            info['outType'] = outImgTypes[i]
            outBname = 'LTV3_{indexID}_{dateRange}_{tile}_{nVert}_{processingDate}_{outType}'.format(**info)
            outFile = os.path.join(tileDir, outBname)
            if os.path.exists(outFile):
                continue
            bands = ' -b '+' -b '.join([str(band) for band in bandRanges[i]])
            cmd = 'gdal_translate -q --config GDAL_DATA "/usr/lib/anaconda/share/gdal" -of GTiff -co "TILED=YES" -co "INTERLEAVE=BAND" -co "BIG_TIFF=YES" -a_srs ' + proj + bands + ' -projwin '+projwin+' '+ vrtFile + ' ' + outFile
            if count is None:
                commands.append((outFile, cmd))
            else:
                commands.append((outFile, cmd, c, nCommands, startTime))
                c += 1
        feature = layer.GetNextFeature()
        i += 1
    if i != nTiles:
        raise RuntimeError('%s is not a valid filter string for clipFile %s' % (filterString, clipFile))
    
    return commands'''
def getDecomposeCommands(chunkDir, runInfoTxt, bandInfoTxt, outDir, clipFile, tileIdField='name', proj='EPSG:5070', nTiles=9, count=None, nTasks=1, nCommands=None, startTime=None, overwrite=False):
    
    t0 = time.time()

    if not os.path.isdir(outDir):
        warnings.warn('outDir does not exist... Creating outDir %s' % outDir)
        os.mkdir(outDir)
    
    # get info about the GEE run
    runInfo = pd.read_csv(runInfoTxt).loc[0]
    runInfo['dateRange'] = (runInfo.startDay + runInfo.endDay).replace('-', '')
    runInfo['segIndex'] = runInfo.segIndex.lower()
    runInfo['processingDate'] = datetime.now().strftime('%Y%m%d')
    bandInfo = pd.read_csv(bandInfoTxt)
    if 'run_name' in runInfo:
        eeTaskName = runInfo.run_name
    elif 'runName' in runInfo:
        eeTaskName = runInfo.runName
    else:# try to infer from 
        raise ValueError('no run name specified in runInfo: %s' % runInfoTxt)
    
    # find the tif chunks
    tifs = []
    for root, dirnames, filenames in os.walk(chunkDir):
        for filename in fnmatch.filter(filenames, eeTaskName + '*.tif'):
            tifs.append(os.path.join(root, filename))    
    
    # make a list of tile tifs
    vrtFile = os.path.join(chunkDir, eeTaskName + '.vrt')
    tileListFile = vrtFile.replace('.vrt', '_filelist.txt')
    with open(tileListFile, 'w') as tileList:
        for piece in tifs:
            tileList.write(piece + '\n')

    # create vrt
    tempFile = os.path.join(chunkDir, 'tempt.txt')
    # direct stdout to file so it buildvrt doesn't print to screen
    cmd = 'gdalbuildvrt -input_file_list %s %s > %s' % (tileListFile, vrtFile, tempFile)
    subprocess.call(cmd, shell=True)
    os.remove(tempFile)

    # make a list of band ranges for each out type
    outImgTypes = bandInfo.name.unique()
    tiles = attributes_to_df(clipFile)
    nTiles = len(tiles[tiles[runInfo.featureKey] == runInfo.featureValue])
    nCommands = len(outImgTypes) * nTiles * nTasks
    
    # Check that the directory structure exists. If not, build it.
    if not os.path.isdir(outDir):
        os.mkdir(outDir)
    indexDir = os.path.join(outDir, '{segIndex}_{dateRange}'.format(**runInfo))
    if not os.path.isdir(indexDir):
        os.mkdir(indexDir)
    
    # Open the tile layer and loop through each tile within this region
    #   Use attribute value instead of geographic bounds of vrt because GEE
    #   often includes nodata buffer of unpredctable size
    dataset = ogr.Open(clipFile)
    layer = dataset.GetLayer()
    filterString = "{featureKey} = {featureValue}".format(**runInfo)
    layer.SetAttributeFilter(filterString)
    
    #import pdb; pdb.set_trace()
    if layer.GetFeatureCount() != nTiles:
        raise RuntimeError('%s is not a valid filter string for clipFile %s' % (filterString, clipFile))
    
    feature = layer.GetNextFeature()
    commands = []
    c = count
    while feature: 
        tileID = feature.GetField(tileIdField)
        tileDir = os.path.join(indexDir, tileID)
        if not os.path.isdir(tileDir):
            os.mkdir(tileDir)

        # format the exent as -projwin arguments for gdal translate
        extent = feature.GetGeometryRef().GetEnvelope()
        projwin = '{} {} {} {}'.format(extent[0], extent[3], extent[1], extent[2])    
     
        # make a list of all the gdal_translate commands needed for the ee conus chunk
        # loop through the datasets and pull them out of the mega stack and write them to the define outDir
        info = runInfo.to_dict()
        info['tile'] = tileID
        for outType in outImgTypes:
            info['outType'] = outType
            outBname = 'LTV3_{segIndex}_{dateRange}_{tile}_{nVerts}_{processingDate}_{outType}.tif'.format(**info)
            outFile = os.path.join(tileDir, outBname)
            if os.path.exists(outFile) and not overwrite:
                continue
            bands = bandInfo.loc[bandInfo.name == outType, 'band']
            bandStr = ' -b ' + ' -b '.join([str(b) for b in bands])
            cmd = 'gdal_translate -q --config GDAL_DATA "/usr/lib/anaconda/share/gdal" -of GTiff -co "TILED=YES" -co "INTERLEAVE=BAND" -co "BIGTIFF=YES" -a_srs ' + proj + bandStr + ' -projwin '+projwin+' '+ vrtFile   + ' ' + outFile
            if count is None:
                commands.append((outFile, cmd))
            else:
                commands.append((outFile, cmd, c, nCommands, startTime))
                c += 1
        
        feature = layer.GetNextFeature()
    
    #os.remove(tempFile) #remove it now because subprocess doesn;t block
    return commands


def clipAndDecompose(chunkDir, outDir, clipFile, searchStr='', njobs=1, tileIdField='name', proj='EPSG:5070', returnOutDirs=False, overwrite=False):
    
    t0 = time.time()
    njobs = int(njobs)
    # find all the runInfo txt files
    runInfoFiles = glob.glob(os.path.join(chunkDir, searchStr + '*run_info.csv'))
    commands = []
    c = 1
    nTasks = len(runInfoFiles)
    #import pdb; pdb.set_trace()
    for runInfoTxt in runInfoFiles:
        bandInfoTxt = runInfoTxt.replace('run_info', 'band_info')
        if len(glob.glob(runInfoTxt.replace('-run_info.csv', '*.tif'))) == 0:
            nTasks -= 1
            continue
        cmds = getDecomposeCommands(chunkDir, runInfoTxt, bandInfoTxt, outDir, clipFile, count=c, nTasks=nTasks, startTime=t0, overwrite=overwrite)
        c =+ len(cmds)
        commands.extend(cmds)

    #import pdb; pdb.set_trace()
    
    if njobs > 1:
        pool = multiprocessing.Pool(njobs)
        pool.map(callTranslate, commands, chunksize=1)
        pool.close()
        pool.join()
        
    else:
        for cmd in commands: callTranslate(cmd)
    
    # Delete source files
    #shutil.rmtree(chunkDir)
    
    print '\n\nTotal processing time: %.1f minutes' % ((time.time() - t0)/60)

