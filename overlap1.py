#Overlap test

#Sample information of two Chromosomes: name, start, and stop location.
dna1 = ['chr1', 10, 30]
dna2 = ['chr1', 18, 19]

#If sequence is similar and any overlap exists, state whether this is overlap
if dna1[0] == dna2[0]:
	if ((dna1[1] >= dna2[1] and dna1[1] <= dna2[2]) or (dna2[1] >= dna1[1] and dna2[1] <= dna1[2])) : #Checks if dna1 start is in range/touching second sequence and vice versa for dna2 start
		print ("overlap exists") 
	else:
		print('regions do not overlap')
else:
	print("Regions are not on the same chromosome")    #If the starting chromosomes arent the same the program doesn't bother checking for overlap. 
