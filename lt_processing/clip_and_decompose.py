#!/usr/bin/env python
"""
Utility function to call clipAndDecompose() from lt_gee_processing.py.

Usage:
    clip_and_decompose.py <chunkDir> <outDir> <clipFile> [--searchStr=<str>] [--njobs=<int>] [--tileIdField=<str>] [--proj=<str>] [--returnOutDirs=<bool>] [-o | --overwrite] [--downloadLog=<str>] [--maxWait=<int>]
    clip_and_decompose.py -h | --help


Options:
    -h --help               Show this screen.
    --searchStr=<str>       Glob-style pattern to find run_info.csv files. Script searches for csvs with searchStr + "*run_info.csv". Default is "".
    --njobs=<int>           Integer number of cores to use[default: 1]
    --tileIdField=<str>     Field in clipFile containing tile IDs. [default: name]
    --proj=<str>            OGR-interpretable projection string [default: EPSG:5070]
    --returnOutDirs=<bool>  Boolean indicating whether to return lowest level output subdirectories [default:'']
    -o --overwrite          Boolean indicating whether to write over existing files
    --downloadLog=<str>     CSV file containing information about a current LT GEE run. If either a path of 'True' specified, the script will wait for tasks to finish downloading frwhile lt_gee_processing.listenAndDownload() is running. This file is automatically created by . Either a path to a log file or
    --maxWait=<int>         Number of seconds to wait for a update from lt_gee_processing.listenAndDownload(). If exceeded, an exception is raised
"""


import sys
import os
import re
import time
import multiprocessing
import pandas as pd
import docopt

import lt_gee_processing as ltee



def main(chunkDir, outDir, clipFile, searchStr='', njobs=1, tileIdField='name', proj='EPSG:5070', returnOutDirs=False, overwrite=False, downloadLog=None, maxWait=14400):
    
    # Set up a pool and a queue to wait for downloads to finish
    if not os.path.isdir(outDir):
        os.mkdir(outDir)
    
    if not os.path.exists(chunkDir):
        raise IOError('chunkDir does not exist: %s' % chunkDir)
    if downloadLog:
        maxWait = int(maxWait)
        
        mgr = multiprocessing.Manager()
        queue = mgr.Queue()
        pool = multiprocessing.Pool(njobs)
        #pool.apply_async(ltee._callTranslateFromQueue, (queue, ))
        
        # Check if the log file exists
        if os.path.isfile(downloadLog):
            log = pd.read_csv(downloadLog)
        else:#assume the log is in chunkDir with the default name
            downloadLog = os.path.join(chunkDir, 'download_log.csv')
            if not os.path.isfile(downloadLog): 
                raise IOError('download log does not exist: %s' % downloadLog)
            log = pd.read_csv(downloadLog)
        imgTasks = log.loc[(log.task_type == 'EXPORT_IMAGE')]
        nTasks = len(imgTasks)
        
        # read the download log
        FAILED_STATES = ['FAILED', 'CANCELLED', 'CANCEL_REQUESTED']
        inQueue = []
        t0 = time.time()
        nCommands = 0
        while not imgTasks.description.isin(inQueue).all():
            # read the log again
            try:
                log = pd.read_csv(downloadLog)
            except:
                continue
            #print log
            # Get all images that finished since the last iteration
            #   They will:
            #   1. be and image
            #   3. not already be in the queue
            imgTasks = log.loc[(log.task_type == 'EXPORT_IMAGE') &\
                                ~log.description.isin(inQueue)]
            csvTasks = log.loc[(log.task_type == 'EXPORT_FEATURES')]
            #import pdb; pdb.set_trace()
            
            for i, task in imgTasks.loc[(imgTasks.state == 'COMPLETED') & imgTasks.downloadDone].iterrows():
                # Find the csvs from this run
                
                csvInfo = csvTasks[[task.description in d for d in csvTasks.description]]
                # There should be 2: one with run info and one with band info
                #import pdb; pdb.set_trace()
                if len(csvInfo) == 2:
                    runInfoTask = log.loc[log.description == (task.description + '-run_info'), 'description'].iloc[0]
                    bandInfoTask = log.loc[log.description == (task.description + '-band_info'), 'description'].iloc[0]
                    runPath = os.path.join(chunkDir, runInfoTask + '.csv')
                    bandPath = os.path.join(chunkDir, bandInfoTask + '.csv')
                    commands = ltee.getDecomposeCommands(chunkDir, runPath, bandPath, outDir, clipFile, tileIdField=tileIdField, proj=proj, nTasks=len(imgTasks), startTime=t0, overwrite=overwrite)
                    for i, (outFile, cmd) in enumerate(commands):
                        queue.put([outFile, cmd, i + nCommands, nCommands, t0])
                        #print queue.qsize()
                        pool.apply_async(ltee._callTranslateFromQueue, (queue, ))
                    #print queue.qsize()
                    nCommands += len(commands)
                    inQueue.append(task.description)
                    #sys.stdout.write('\rSubmitted commands for %s of %s tasks' % (len(inQueue), nTasks))
                    #sys.stdout.flush()
                    #import pdb; pdb.set_trace()
                
                # If any of them failed
                elif csvInfo.state.isin(FAILED_STATES).any():
                    # Say that it's in the queue
                    inQueue.append(task.description)
            
            # If any of the image tasks failed, say it's in the queue
            for i, task in log.loc[(log.task_type == 'EXPORT_IMAGE') & log.state.isin(FAILED_STATES)].iterrows():
                if not task.description in inQueue:
                    inQueue.append(task.description)
                
            # check the mtime of the log file and if it's too old, break
            last_modified = os.stat(downloadLog).st_mtime
            if (time.time() - last_modified) > maxWait:
                raise RuntimeError('Last task update exceeded maxWait time of %.1f hours' % (maxWait/3600.))
        import pdb; pdb.set_trace()
        # Close all processes
        for _ in range(njobs): 
            queue.put(None) # add kill switch to break of while loop in downloadFrom
        pool.close()
        pool.join()
        
    # Presumably, all files are done downloading, so process all of them at once
    else:
        ltee.clipAndDecompose(**args)
        

if __name__ == '__main__':
    
    # Any args that don't have a default value and weren't specifed will be None
    cl_args = {k: v for k, v in docopt.docopt(__doc__).iteritems() if v is not None}

    # get rid of extra characters from doc string and 'help' entry
    args = {re.sub('[<>-]*', '', k): v for k, v in cl_args.iteritems()
            if k != '--help' and k != '-h'}  
    
    try:
        args['njobs'] = int(args['njobs'])
    except:
        raise ValueError('njobs not understood: ' + args['njobs'])

    if 'returnOutDirs' in args:        
        if args['returnOutDirs'].lower() == 'false':
            args['returnOutDirs'] = False
    
    if 'o' in args or 'overwite' in args:
        args['overwrite'] = True
    
    sys.exit(main(**args))