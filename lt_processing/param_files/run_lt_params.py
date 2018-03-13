featureValue = 'r0801';
featureCol = '1Xk_FlopYbVU1fiW44rjrT0KjtNTh51AK7qtq6dH0';
featureKey = 'epa_code';
startYear = 1984; # what year do you want to start the time series 
endYear = 2016; # what year do you want to end the time series
startDay = ['04-15', '06-15', '08-15']; # what is the beginning of date filter | month-day
endDay =   ['07-15', '09-15', '11-15']; # what is the end of date filter | month-day
mosaicType = "medoid"; # how to make annual mosaic - options: "medoid", "targetDay" 
targetDay = None ; # if running "targetDay" mosaic, what day of year should be the target
outProj = 'EPSG:5070'; # what should the output projection be? 'EPSG:5070' is North American Albers
gDriveFolder = 'ltee_stem_trials_' + featureValue; # what is the name of the Google Drive folder that you want the outputs placed in
affine = [30.0, 0, 15.0, 0, -30.0, 15.0];
aoiBuffer = 300;

indexDict = {'NBR': -1,
             'B5': 1,
             'NDVI': -1, 
             'TCB': 1,
             'NDSI': -1,
             'TCG': -1,
             'B3': 1}

# define the segmentation parameters - see paper (NEED CITATION)
run_params = { 
    'maxSegments': 6,
    'spikeThreshold': 0.9,
    'vertexCountOvershoot': 3,
    'preventOneYearRecovery': True,
    'recoveryThreshold': 0.25,
    'pvalThreshold': 0.05,
    'bestModelProportion': 0.75,
     'minObservationsNeeded': 6
     }
  
