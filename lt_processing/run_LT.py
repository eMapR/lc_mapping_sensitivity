# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 17:40:51 2018

@author: shooper
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 12:11:52 2018

@author: shooper
"""


import os
import sys
import time
import shutil
import warnings
import pandas as pd
import ee
ee.Initialize()


def read_params(txt, sep='='):
    '''
    Return a dictionary and a dataframe from parsed parameters in txt
    '''
    if not os.path.exists(txt):
        print 'Param file does not exist:\n%s' % txt
        return None

    # Read in the rest of the text file line by line
    #
    try:
        with open(txt) as f:
            input_vars = [line.split(sep) for line in f]
    except:
        print 'Problem reading parameter file:\n', txt
        return None

    # Set the dictionary key to whatever is left of the ";" and the value
    #   to whatever is to the right. Strip whitespace too.
    '''d = {}
    for var in input_vars:
        if len(var) == 2:
            d[var[0].replace(" ", "")] =\
                '"{0}"'.format(var[1].strip(" ").replace("\n", ""))'''
    
    for var in input_vars:
        if len(var) == 2:
            name = var[0].replace(' ', '')
            value = var[1].strip(' ').replace('\n', '').split('#')[0]
            exec ('%s = str("%s")' % (name, value))
            exec ('global %s' % name)
    import pdb; pdb.set_trace()
    print '\nParameters read from:\n', txt, '\n'
    
    #return d


# #######################################################################################
# ###### ANNUAL SR TIME SERIES STACK BUILDING FUNCTIONS #################################
# #######################################################################################

# ------ DEFINE L8 to L7 ALIGN FUNCTION ---------

# slope and intercept citation: Roy, D.P., Kovalskyy, V., Zhang, H.K., Vermote, E.F., Yan, L., Kumar, S.S, Egorov, A., 2016, Characterization of Landsat-7 to Landsat-8 reflective wavelength and normalized difference vegetation index continuity, Remote Sensing of Environment, 185, 57-70.(http:#dx.doi.org/10.1016/j.rse.2015.12.024); Table 2 - reduced major axis (RMA) regression coefficients
def harmonizationRoy(oli):
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
    ltVertStack = LTresult.arrayMask(vmask).arraySlice(0, 0, 3).addBands(zeros).toArray(1).arraySlice(1, 0, run_params['maxSegments']+1).arrayFlatten(lbls, '')#'''
                      #.toShort();'''    
    '''ltVertStack = LTresult.arrayMask(vmask)\
                          .arraySlice(0, 0, 3)\
                          .addBands(zeros)\
                          .toArray(1)\
                          .arraySlice(1, 0, run_params['maxSegments']+1)\
                          .arrayFlatten(lbls, '')
                          #.toShort();'''
    return ltVertStack;


# ------ RUN LANDTRENDR -----------------------------------------------------------------
def main():
    
    # get geometry stuff
    # # load the aoi
    #... filter by park - left arg is the key name, right arg is the value to match
    aoi = ee.FeatureCollection('ft:'+featureCol)\
            .filter(ee.Filter.stringContains(featureKey, featureValue))\
            .geometry()\
            .buffer(aoiBuffer);
    
    box = aoi.bounds(); # get the bounds from the drawn polygon
    
    
    # make a dummy collection for filling missing years (if there are any)
    dummyCollection = ee.ImageCollection([ee.Image([0,0,0,0,0,0]).mask(ee.Image(0))]); # make a dummy collection to fill in potentially missing years
    indexList = [['NBR', -1], ['B5', 1], ['NDVI', -1], ['TCB', 1], ['NDSI', -1], ['TCG', -1], ['B3', 1]];
    
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
    #index, indexFlip = '0', 0
    global index 
    global indexFlip
    
    tasks = []
    for i in range(len(startDay)):
        # build the annual SR collection
        annualSRcollection = buildMosaicCollection(startYear, endYear, startDay[i], endDay[i], box, mosaicType, targetDay, dummyCollection); # put together the cloud-free medoid surface reflectance annual time series collection
        for j in range(len(indexList)):
            #get the index and the index flipper
            index, indexFlip = indexList[j]; # pull out the index to segment on and value to flip index or not
            
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
            nVert = int(run_params['maxSegments'])+1;
            fileNamePrefix = featureKey+'-'+featureValue+'-'+index+'-'+str(nVert)+'-'+str(startYear)+str(endYear) + '-' + startDay[i].replace('-', '') + endDay[i].replace('-', '');
            export_params = {'region': aoi,
                             'folder': gDriveFolder,
                             'fileNamePrefix': fileNamePrefix,
                             'crs': outProj,
                             'crsTransform': affine,
                             'maxPixels': 1e13
                             }
    
            task = ee.batch.Export.image.toDrive(allStack.clip(aoi), fileNamePrefix, config=export_params);
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

    # Download
    # Figure out how to give authentication key without manual intervention
    #   prbably some way to authenticate at the beginning of the script and keep oauth token

    

if __name__ == '__main__':
    
    sys.path.append(os.path.dirname(sys.argv[1]))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore') # raises SyntaxWarning because of *
        from ee_test_params import * # brings all variables into namespace
    sys.exit(main())#'''
    