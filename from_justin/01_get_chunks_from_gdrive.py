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
def download_files(fileName, outDirPath):
  print(fileName['title'])
  getFile = drive.CreateFile({'id': fileName['id']})
  #getFile.GetContentFile(outDirPath+fileName['title']) # Download file
  ################### added by SMH 2/22/18 ##################################
  try:
      getFile.GetContentFile(outDirPath+fileName['title']) # Download file
  except Exception as e:
      print('Problem with file %s' % fileName['title'])
      logFile = os.path.join(outDirPath, 'error_log.txt')
      with open(logFile, 'a') as f:
          f.write('%s\nfile: %s\n\n' % (e, fileName['title']))
  ###########################################################################

# get the arguments
args = sys.argv
gDirName = args[1]
outDirPath = args[2]
############### added by SMH 2/22/18 ##############
if len(args)>=3:
    njobs = int(args[3])
else:
    njobs = 4
###################################################
# example of inputs
#gDirName = "ltgee_nbr_ltr_tile28_06200910_prevOneYearOffTest_ltr_stack"
#outDirPath = "/vol/v1/proj/ltee_vs_ltidl/raster/region_28/ltgee/seg_prevOneYearOff/prep"


# make sure the paths end in '/'
if outDirPath[-1] != '/':
  outDirPath += '/'

os.chdir('/vol/v1/general_files/script_library/earth_engine/') #GoogleAuth looks in here for an authorization file - could pass the file as an argument and the get the os.path.dirname

# authenticate gDrive application and request access to gDrive account
gauth = GoogleAuth()
gauth.LocalWebserverAuth() # creates local webserver and auto handles authentication.
drive = GoogleDrive(gauth)     
                        
# find files in the specified gDrive folder
gDir = drive.ListFile({'q': "mimeType='application/vnd.google-apps.folder' and title contains '"+gDirName+"'"}).GetList()
if len(gDir) == 1:
  fileList = drive.ListFile({'q': "'"+gDir[0]['id']+"' in parents and title contains '.tif'"}).GetList()
  
  # create the output folder if it does not already exist
  if not os.path.isdir(outDirPath):
    os.mkdir(outDirPath)
  
  ############ added 2/1/18 by SMH ##########################################
  # check which files have already been downloaded fully
  fileInfo = pd.DataFrame(fileList)
  existing = fileInfo[fileInfo.title.apply(
              lambda z: os.path.exists(os.path.join(outDirPath, z)))]
  unfinished = existing[existing.fileSize.astype(int) !=\
                  existing.title.apply(lambda z: os.stat(outDirPath + z).st_size)]
  # find files that either didn't finish or that don't exist locally
  fileList = pd.concat([unfinished, fileInfo[~fileInfo.title.isin(existing.title)]])
  fileList = [dict(r) for i, r in fileList.iterrows()]#back to what download_files expects
  
  ############################################################################
  # wait 10 seconds to start - if the folder is created in the line above
  # then the download won't start, rerunning the script will get it to start
  # could be that the folder is not fully registered before pool.map(func, fileList) 
  # is called
  time.sleep(10)

  # loop through downloading the files in parallel
  '''for fn in fileList:
    download_files(fn, outDirPath)'''

  if __name__ == '__main__':
      pool = multiprocessing.Pool(processes=njobs) 
      func = partial(download_files, outDirPath=outDirPath)
      pool.map(func, fileList)  
      pool.close()  
      pool.join()

else:
    raise RuntimeError('Could not find %s' % gDir)