#!/usr/bin/python                                                                                                                                                                           
#pipeline_template.py                                                                                                                                                                       

'''                                                                                                                                                                                         
The MIT License (MIT)                                                                                                                                                                       
                                                                                                                                                                                            
Copyright (c) 2015 Charles Lin                                                                                                                                                              
                                                                                                                                                                                            
Permission is hereby granted, free of charge, to any person obtaining a copy                                                                                                                
of this software and associated documentation files (the "Software"), to deal                                                                                                               
in the Software without restriction, including without limitation the rights                                                                                                                
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell                                                                                                                   
copies of the Software, and to permit persons to whom the Software is                                                                                                                       
furnished to do so, subject to the following conditions:                                                                                                                                    
                                                                                                                                                                                            
The above copyright notice and this permission notice shall be included in                                                                                                                  
all copies or substantial portions of the Software.                                                                                                                                         
                                                                                                                                                                                            
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR                                                                                                                  
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,                                                                                                                    
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE                                                                                                                 
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER                                                                                                                      
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,                                                                                                               
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN                                                                                                                   
THE SOFTWARE.                                                                                                                                                                               
'''

#generic pipeline template for human data                                                                                                                                                   


#==========================================================================                                                                                                                 
#=============================DEPENDENCIES=================================                                                                                                                 
#==========================================================================                                                                                                                 


import sys
sys.path.append('/storage/cylin/bin/pipeline/')

import pipeline_dfci
import utils
import string
import time
import datetime

from collections import defaultdict

#==========================================================================                                                                                                                 
#============================PARAMETERS====================================                                                                                                                 
#==========================================================================                                                                                                                 



projectName = 'proving_ground'
dataFile = '/storage/cylin/grail/projects/%s/chordoma_rna_data_table.txt' % (projectName)

#project folders                                                                                                                                                                            
projectFolder = '/storage/cylin/grail/projects/%s/' % (projectName) #PATH TO YOUR PROJECT FOLDER                                                                                            

#standard folder names                                                                                                                                                                      
gffFolder ='%sgff/' % (projectFolder)
macsFolder = '%smacsFolder/' % (projectFolder)
macsEnrichedFolder = '%smacsEnriched/' % (projectFolder)
mappedEnrichedFolder = '%smappedEnriched/' % (projectFolder)
mappedFolder = '%smappedFolder/' % (projectFolder)
wiggleFolder = '%swiggles/' % (projectFolder)
metaFolder = '%smeta/' % (projectFolder)
genome = 'hg19_ercc'
annotFile = '/storage/cylin/bin/pipeline/annotation/hg19_refseq.ucsc'

#making folders                                                                                                                                                                             
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder]

#for folder in folderList:
#    pipeline_dfci.formatFolder(folder,True)

#==========================================================================                                             
#=======================LOADING DATA ANNOTATION============================                                             
#==========================================================================                                             

##THIS SECTION LOADS A DATA TABLE.  MUST BE UNCOMMENTED FOR REST OF CODE TO WORK                                        


#LOADING THE DATA TABLE                                                                                                 
dataDict = pipeline_dfci.loadDataTable(dataFile)

print(dataDict.keys())

pipeline_dfci.summary(dataFile)

#print(dataDict)

#==========================================================================                                                                                                                       
#=======================HISAT2 ALIGNMENT FASTQ TO BAM======================
#========================================================================== 



#pipeline_dfci.mapHisat(dataFile,namesList=[],useSRA=False,pCount=16,Launch=True)

#==========================================================================                                                                                                                                 
#============================Cufflinks=====================================                                                                                                                                 
#========================================================================== 

gtfFile='/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes_ercc.gtf'
analysisName='chordoma_rnaseq'
cufflinksFolder='/storage/cylin/grail/projects/proving_ground/cufflinks/'

pipeline_dfci.makeCuffTableSlurm(dataFile,analysisName,gtfFile,cufflinksFolder,groupList=[['THZ1_00H_1','THZ1_00H_2','THZ1_00H_3'],['THZ1_04H_1','THZ1_04H_2','THZ1_04H_3'],['THZ1_08H_1','THZ1_08H_2','THZ1_08H_3'],['THZ1_12H_1','THZ1_12H_2','THZ1_12H_3'],['THZ1_24H_1','THZ1_24H_2','THZ1_24H_3']],bashFileName = '')
