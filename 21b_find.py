#This program takes input of accension number for a gene 
#searches for the chromsome#, (+/-) strand, and location of 3' UTR 
#from a tab delimited refGene.txt file 

import sys
sys.path.append('/storage/cylin/bin/pipeline')
import utils
import os
import re
import itertools
import difflib 
from difflib import SequenceMatcher 
#parses the data table
table = utils.parseTable('/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/refGene.txt', '\t')


for line in table[:]:
        if line[1] == sys.argv[1]:
                                
                #if (+) strand make sure you use the 3' end
                if line[3] == '+':
                        
                        print line[1]
                        print line[2]
                        print line[3]
                        line[7] = int(line[7])
                        
                        #add a postion after end of coding sequence within the transcript 
                        UTR_start = (line[7] + 1)
                        print UTR_start
                        print '^ is the start postion of the 3\' UTR'
                        print line[5] 
                        print '^ is the stop postion of the 3\' UTR' 
			
                        #print line[4] + '<- this is where the transcript starts'
                        #print line[6] + '<- this is where the coding seq starts'
                        
                        
                        #turns values into ints for use with fetchSeq function
                        UTR_start = int(UTR_start)
                        line[5] = int(line[5])
                        UTR_stop = line[5]
                        dna = utils.fetchSeq('/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg38/Sequence/Chromosomes/', line[2], UTR_start, UTR_stop)
                        line[4] = int(line[4])
                        transcript1 = utils.fetchSeq('/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg38/Sequence/Chromosomes/', line[2], line[4], line[5])
                        transcript1 = transcript1.upper()
                        dna1 = dna
                  
                        
                                
                                        
                  
                  
 
                #if the strand was - make sure to grab the 3' UTR from the correct side
                elif line[3] == '-':
                        print line[2]
                        print line[3]
                        line[6] = int(line[6])
                        UTR_start2 = (line[6] - 1)
                        print UTR_start2
                        print '^ is the start position of the (-) strand 3\' UTR'
                        print line[4]
                        print '^ is the stop postion of the (-) strand 3\' UTR'

                        UTR_start = int(UTR_start2)
                        line[4] = int(line[4])
                        UTR_stop = line[4]
                        dna = utils.fetchSeq('/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg38/Sequence/Chromosomes/', line[2], UTR_stop, UTR_start)
                        
                        dna = utils.revComp(dna)
                        print dna

               
        elif line[1] == sys.argv[2]:

                #if (+) strand make sure you use the 3' end                                                                                                                                                 
                if line[3] == '+':

                        print line[1]
                        print line[2]
                        print line[3]
                        line[7] = int(line[7])

                        #add a postion after end of coding sequence within the transcript                                                                                                                   
                        UTR_start = (line[7] + 1)
                        print UTR_start
                        print '^ is the start postion of the 3\' UTR'
                        print line[5]
                        print '^ is the stop postion of the 3\' UTR'

                        #print line[4] + '<- this is where the transcript starts'
                        #print line[6] + '<- this is where the coding seq starts'

                        #turns values into ints for use with fetchSeq function                                                                                                                              
                        UTR_start = int(UTR_start)
                        line[5] = int(line[5])
                        UTR_stop = line[5]
                        dna = utils.fetchSeq('/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg38/Sequence/Chromosomes/', line[2], UTR_start, UTR_stop)
                        line[4] = int(line[4])
                        transcript2 = utils.fetchSeq('/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg38/Sequence/Chromosomes/', line[2], line[4], line[5])
                        transcript2 = transcript2.upper()
                        dna2 = dna
                 





                #if the strand was - make sure to grab the 3' UTR from the correct side                                                                                                                     
                elif line[3] == '-':
                        print line[2]
                        print line[3]
                        line[6] = int(line[6])
                        UTR_start2 = (line[6] - 1)
                        print UTR_start2
                        print '^ is the start position of the (-) strand 3\' UTR'
                        print line[4]
                        print '^ is the stop postion of the (-) strand 3\' UTR'

                        UTR_start = int(UTR_start2)
                        line[4] = int(line[4])
                        UTR_stop = line[4]
                        dna = utils.fetchSeq('/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg38/Sequence/Chromosomes/', line[2], UTR_stop, UTR_start)
                        dna2 = dna
                        dna2 = utils.revComp(dna)
                        print dna2
dna1 = dna1.upper()
dna2 = dna2.upper()

i = (transcript1)
j = (transcript2)
match = SequenceMatcher(None, dna1, dna2, autojunk=False).find_longest_match(0, len(dna1), 0, len(dna2))                                                                                                               

print match                                                                                                                                                                                                
print (dna1[match.a: match.a + match.size])

if match.size != 21 or match.size <= 21:
        print 'The match is not of optimum length. \nNow searching for a longer match in the entire transcript.'


        match = SequenceMatcher(None, i, j, autojunk=False).find_longest_match(0, len(i), 0, len(j))  
                                       
        print match                                                                                                           
        print (i[match.a: match.a + match.size])
else:
        print 'Match is of optimum length'


