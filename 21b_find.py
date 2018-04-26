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
                        line[7] = int(line[7])
                        UTR_start = (line[7]+1)
                        print '( ' + line[2] + ', ' + line[3] + ' ,', UTR_start, line[5] + ' )'
                        print '^ represents the 3 prime UTR region'
                         
                        UTR_start = int(UTR_start)
                        line[5] = int(line[5])
                        UTR_stop = line[5]
                        dna1 = utils.fetchSeq('/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg38/Sequence/Chromosomes/', line[2], UTR_start, UTR_stop)
                        dna1 = dna1.upper()
                        line[4] = int(line[4])
                        transcript1 = utils.fetchSeq('/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg38/Sequence/Chromosomes/', line[2], line[4], line[5])
                        transcript1 = transcript1.upper()
                        
                #if the strand was - make sure to grab the 3' UTR from the correct side
                elif line[3] == '-':
                        
                        print line[1]
                        line[6] = int(line[6])
                        UTR_start = (line[6]-1)
                        print '( ' + line[2] + ', ' + line[3] + ' , ' + line[4],',', (line[6]-1),')'
                        print '^ represents the 3 prime region'
                        UTR_start = int(UTR_start)
                        line[4] = int(line[4])
                        line[5] = int(line[5])
                        UTR_stop = line[4]
                        dna1 = utils.fetchSeq('/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg38/Sequence/Chromosomes/', line[2], UTR_stop, UTR_start)
                        dna1 = dna1.upper()
                        dna1 = utils.revComp(dna1)
                        transcript1 = utils.fetchSeq('/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg38/Sequence/Chromosomes/', line[2], UTR_stop, line[5])
                        transcript1 = transcript1.upper()
                        transcript1 = utils.revComp(transcript1)                        
                      
        elif line[1] == sys.argv[2]:

                #if (+) strand make sure you use the 3' end                                                                                                                                                 
                if line[3] == '+':

                        print line[1]
                        line[7] = int(line[7])
                        UTR_start2 = (line[7]+1)
                        print '( ' + line[2] + ', ' + line[3] + ' ,', UTR_start2,',', line[5] + ' )'
                        print '^ represents the 3 prime UTR region'
                                                                                              
                        UTR_start2 = int(UTR_start2)
                        line[5] = int(line[5])
                        UTR_stop2 = line[5]
                        dna2 = utils.fetchSeq('/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg38/Sequence/Chromosomes/', line[2], UTR_start2, UTR_stop2)
                        dna2 = dna2.upper()
                        line[4] = int(line[4])
                        transcript2 = utils.fetchSeq('/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg38/Sequence/Chromosomes/', line[2], line[4], line[5])
                        transcript2 = transcript2.upper()
                        
                #if the strand was - make sure to grab the 3' UTR from the correct side                                                                                                                     
                elif line[3] == '-':
                        
                        print line[1]
                        line[6] = int(line[6])
                        UTR_start = (line[6]-1)
                        print '( ' + line[2] + ', ' + line[3] + ' , ' + line[4],',', (line[6]-1),')'
                        print '^ represents the 3 prime region'
                        UTR_start = int(UTR_start)
                        line[4] = int(line[4])
                        line[5] = int(line[5])
                        UTR_stop = line[4]
                        dna2 = utils.fetchSeq('/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg38/Sequence/Chromosomes/', line[2], UTR_stop, UTR_start)
                        dna2 = dna2.upper()
                        dna2 = utils.revComp(dna2)
                        transcript2 = utils.fetchSeq('/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg38/Sequence/Chromosomes/', line[2], UTR_stop, line[5])
                        transcript2 = transcript2.upper()
                        transcript2 = utils.revComp(transcript2)

print 'Now searching for match in 3 prime UTR'

match = SequenceMatcher(None, dna1, dna2, autojunk=False).find_longest_match(0, len(dna1), 0, len(dna2))                                                                                              
print match                                                                                                                                                                                           
print (dna1[match.a: match.a + match.size])

if match.size != 21 or match.size <= 21:
        print 'The match is not of optimum length. \nNow searching for a longer match in the entire transcript.'

        match = SequenceMatcher(None, transcript1, transcript2, autojunk=False).find_longest_match(0, len(transcript1), 0, len(transcript2))  
                                       
        print match                                                                                                           
        print (transcript1[match.a: match.a + match.size])
else:
        print 'Match is of optimum length'


