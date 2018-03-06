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


# define function to download tif files found in gDrive folder - called by multiprocessing.Pool().map()
def download_files(fileName, outDirPath):
  print('   file: '+fileName['title'])
  getFile = drive.CreateFile({'id': fileName['id']})
  getFile.GetContentFile(outDirPath+fileName['title']) # Download file

# get the arguments
#args = sys.argv
#gDirName = args[1]
#outDirPath = args[2]

# example of inputs
gDirName = "pecora_time_test_17"
outDirPath = "/vol/v1/proj/stem_improv_paper/test_regions/tiles/staging/download17"

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
#if len(gDir) == 1: # TODO else print problem and exit
fileList = drive.ListFile({'q': "'"+gDir[0]['id']+"' in parents and title contains '.tif'"}).GetList()

# create the output folder if it does not already exist
if not os.path.isdir(outDirPath):
  os.mkdir(outDirPath)

# wait 10 seconds to start - if the folder is created in the line above
# then the download won't start, rerunning the script will get it to start
# could be that the folder is not fully registered before pool.map(func, fileList) 
# is called
time.sleep(10)

for i, thisFile in enumerate(fileList):
  print("i: "+str(i))
  download_files(thisFile, outDirPath)


# loop through downloading the files in parallel
#pool = multiprocessing.Pool(processes=3) 
#func = partial(download_files, outDirPath=outDirPath)
#pool.map(func, fileList)  
#pool.close()  
  
