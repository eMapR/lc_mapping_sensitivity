
import os
import fnmatch
import subprocess



def run_cmd(cmd):
  print(cmd)  
  return subprocess.call(cmd, shell=True)



# make a list of tile tifs
def make_vrt(inputFiles, vrtFile):
  #vrtFile = chunkDir+name+'_'+'_'.join(os.path.basename(tileFile).split('_')[4:7])+'.vrt'
  print(vrtFile)  
  inputListFile = vrtFile.replace('.vrt', '_filelist.txt')
  inputList = open(inputListFile, 'w')
  for inputFile in inputFiles:
    inputList.write(inputFile+'\n')
  inputList.close()
  
  # create vrt
  cmd = 'gdalbuildvrt -input_file_list '+inputListFile+' '+vrtFile
  subprocess.call(cmd, shell=True)

"""
# get the arguments
chunkDir = sys.argv[1]
outDir = sys.argv[2]  
tileFile = sys.argv[3]
name = sys.argv[4]  
indexID = sys.argv[5]
nVert = int(sys.argv[6])
startYear = int(sys.argv[7])
endYear = int(sys.argv[8])
proj = sys.argv[9]
delete = sys.argv[9]
"""

chunkDir = '/vol/v1/proj/stem_improv_paper/test_regions/tiles/staging/download17'
outDir = '/vol/v1/proj/stem_improv_paper/test_regions/tiles/staging/vrt' #'/vol/v2/conus_tiles/tiles'  #'/vol/v1/general_files/user_files/justin/for_others/fs_karen/data/landtrendr/ycd_test/'



# make sure path parts are right
if chunkDir[-1] != '/':
  chunkDir += '/'
if outDir[-1] != '/':
  outDir += '/'


# find the tif chunks
rasters = []
for root, dirnames, filenames in os.walk(chunkDir):
  for filename in fnmatch.filter(filenames, '*.tif'):
    rasters.append(os.path.join(root, filename))

eeLTfileParts = [thisOne.split('-')[0] for thisOne in rasters]
eeLTfilePartsUni = list(set(eeLTfileParts))
for thisPart in eeLTfilePartsUni:
  theseFiles = [rasters[i] for i, j in enumerate(eeLTfileParts) if j == thisPart]
  for thisFile in theseFiles:
    print(os.path.basename(thisFile))
  print('')
  vrtFile = outDir+os.path.basename(thisPart)+'.vrt'
  make_vrt(theseFiles, vrtFile)
  
  
  







