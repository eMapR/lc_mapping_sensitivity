# -*- coding: utf-8 -*-
"""
Created on Sun Apr 02 09:15:58 2017

@author: braatenj

https://googledrive.github.io/PyDrive/docs/build/html/index.html
https://pypi.python.org/pypi/PyDrive
"""

from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive
import os
import sys
import time
import multiprocessing
from functools import partial
import pandas as pd


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
    
# example of inputs
#gDirName = "ltgee_nbr_ltr_tile28_06200910_prevOneYearOffTest_ltr_stack"
#outDirPath = "/vol/v1/proj/ltee_vs_ltidl/raster/region_28/ltgee/seg_prevOneYearOff/prep"

def main(gDirName, outDirPath, njobs=4):
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

if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:]))