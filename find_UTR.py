#This program takes input of accension number for a gene 
#searches for the chromsome#, (+/-) strand, and location of 3' UTR 
#from a tab delimited refGene.txt file 

import sys
sys.path.append('/home/harrison/Downloads/pipeline-master')

import utils

table = utils.parseTable('refGene.txt', '\t')

for line in table[:]:

	if line[1] == sys.argv[1]:

		if line[3] == '+':
			print line[2]
			print line[3]
			line[7] = int(line[7])
			UTR_start = (line[7] + 1)
			print UTR_start
			print '^ is the start postion of the 3\' UTR'
			print line[5] 
			print '^ is the stop postion of the 3\' UTR' 
			

		elif line[3] == '-':
			print line[2]
			print line[4]
			print '^ is the start position of the (-) strand 3\' UTR'
			line[6] = int(line[6])
			UTR_start2 = (line[6] - 1)
			print UTR_start2
			print '^ is the stop postion of the (-) strand 3\' UTR'

		




