# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 12:35:56 2017

@author: braatenj
"""

import os, fnmatch

searchDir = '/vol/v1/proj/stem_improv_paper/test_regions/tiles/tiles_seg/'

files = []
for root, dirnames, filenames in os.walk(searchDir):
  for filename in fnmatch.filter(filenames, '*'):
    files.append(os.path.join(root, filename))

for thisFile in files:
  dname = os.path.dirname(thisFile)
  bname = os.path.basename(thisFile)
  if 'nbr' not in bname[0:38]:
    newFile = os.path.join(dname, '{0}_nbr_{1}'.format(bname[:28], bname[29:]))
    os.rename(thisFile, newFile)