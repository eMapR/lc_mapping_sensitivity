# -*- coding: utf-8 -*-
"""
Function library for LT processing on Google Earth Engine.

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

import ee
ee.Initialize()

from lthacks import createMetadata

# #######################################################################################
# ###### ANNUAL SR TIME SERIES STACK BUILDING FUNCTIONS #################################
# #######################################################################################


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
# #######################################################################################
# #######################################################################################


'''def launchTask(task, log_file):
    with open(log_file, 'a') as txt:'''
        


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


def run_lt(featureCol, featureKey, featureValue, aoiBuffer, indexDict, startYear, endYear, run_params, startDay, endDay, mosaicType, targetDay, gDriveFolder, outProj, affine):   
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
    for i in range(len(startDay)):
        # build the annual SR collection
        annualSRcollection = buildMosaicCollection(startYear, endYear, startDay[i], endDay[i], box, mosaicType, targetDay, dummyCollection); # put together the cloud-free medoid surface reflectance annual time series collection
        for index, indexFlip in indexDict.iteritems():#range(len(indexList)):
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
    
    
    taskStatus = pd.DataFrame([t.status() for t in tasks])
    taskStatus['active'] = True
    nTasks = float(len(tasks))
    sleep = .5
    t0 = time.time()
    # Wait until all tasks have finished running 
    while taskStatus['active'].any():
        time.sleep(sleep)
        inactive = taskStatus[~taskStatus.active]
        failed = inactive[inactive.state.isin(['FAILED', 'CANCELED', 'CANCEL_REQUESTED'])]
        complete = inactive[inactive.state == 'COMPLETE']
        nFinished = len(complete) + len(failed)
        cumTime = time.time() - t0
        msg = '\r%s of %d (%.1f%%) tasks done | cum. time: %.1f mins' % (nFinished, nTasks, (nFinished/nTasks) * 100, cumTime/60)
        sys.stdout.write(msg)
        sys.stdout.flush()
        taskStatus = pd.DataFrame([t.status() for t in tasks])#retrieve new status
        taskStatus['active'] = [t.active() for t in tasks] #update active column
    
    if len(failed) > 0:
        print '\n\nIDs of failed tasks:\n\t' + '\n\t'.join(failed['id'].tolist())


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
    t0 = time.time()
    njobs = int(njobs)
    outDirPath = os.path.abspath(outDirPath)
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
                cmd = 'gdal_translate -q --config GDAL_DATA "/usr/lib/anaconda/share/gdal" -of GTiff -co "TILED=YES" -co "INTERLEAVE=BAND" -co "BIG_TIFF=YES" -a_srs ' + proj + bands + ' -projwin '+projwin+' '+ vrtFile + ' ' + outFile
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
