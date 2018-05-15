#This is the Cas9 Blocker! 

#Determine sequence match
#Locate pam site frame
#Alter the overlapping PAM sequence without changing the codon-AA

import re                             
import collections
from collections import defaultdict
import os
import sys
sys.path.append('/home/harrison/Downloads/pipeline-master')
import utils
import math


#parses the tab delimited table
table=utils.parseTable('/home/harrison/PPscripts/HS_tools/codon_hg19.txt', '\t')

#make the first row the key and the next 2 rows the values for that key
hg19_bias = defaultdict()
for line in table[1:]:
    hg19_bias[line[0]]=[line[1],float(line[2])]

#This checks if the input is a text
if os.path.isfile(sys.argv[1]) and os.path.isfile(sys.argv[2]):

#This opens the input files
    with open(sys.argv[1]) as f, open(sys.argv[2]) as g:

#This reads the file and makes the text a usuable string                          
        dna1 = f.read().rstrip()                                                
        gRNA = g.read().rstrip()

#If the input was not a file than continue with the code
else:                                                                                
    dna1 = sys.argv[1]
    gRNA = sys.argv[2]  
    pass

#This will replace all 'T's with 'U's if DNA was inputted otherwise    
dna1 = dna1.replace('T', 'U')                                                
gRNA = gRNA.replace('T', 'U')
dna1 = dna1.upper()
gRNA = gRNA.upper() 

def reverse_complement(dna):
                complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A',}
                return ''.join([complement[base] for base in dna[::-1]])

#this gives you a reverse version of the target and guide
rev_gRNA = gRNA[::-1]

revcom_dna1 = reverse_complement(dna1)


#Checks if the match is unique in the original sequence
if dna1.count(gRNA) == 1 and revcom_dna1.count(rev_gRNA) == 0:
    
    print 'This is the target:\n'+ sys.argv[1]
    print 'This is the guide:\n' + sys.argv[2]     

    #Finds the location of the match
    for a in re.finditer(gRNA, dna1):    

        #This makes sure the adjacent 3 bases make a PAM site
        if dna1[a.end() + 1] == 'G' and dna1[a.end() + 2] == 'G':
            
            #This makes the PAM site a PAM variable
            PAM = dna1[a.end()], dna1[a.end() + 1], dna1[a.end() + 2]

            PAM = ''.join(PAM)

            print 'The guide matches the sequence once at: '
            print (a.start(), a.end())

            #This print's out the same codon/PAM in DNA if DNA was given or leaves it as RNA if not.
            if sys.argv[1].find('U') == -1:
                t = PAM.replace('U', 'T')
            elif sys.argv[1].find('T') == -1:
                t = PAM    
            print t + ' is the PAM sequence'
                
            #remainder indirectly gives postion for base 1 in PAM telling us the frame its in
            if (a.end()%3) == 0:
                print'The PAM is in a 123 frame'
                
                #if 123 frame, PAM is the same as the codon that needs to be changed
                print t + ' Is the overlapping codon'
                if hg19_bias[PAM][0] == 'W':
                    print 'Cannot alter tryptophan due to only one codon'
                    quit()
                #print out the amino acid the codon makes
                print t + ' produces ' + hg19_bias[PAM][0]
                    
                #make a list of other codons that have the same amino acid
                codon_list = []
                for i in hg19_bias:
                    if hg19_bias[i][0]==hg19_bias[PAM][0]:
                        codon_list.append(i)

                #finds a different codon with the highest codon bias       
    
                
                for b in codon_list:
                    if b == t:
                        codon_list.remove(b)
                                       
                        f = hg19_bias[PAM]
                        del hg19_bias[PAM]
                        

                        freq_diff = 0
                        freq_diff_list = [] 
                    
                        for i in codon_list:
                            

                            freq_diff = f[1] - hg19_bias[i][1]
                            freq_diff = math.fabs(freq_diff)
                            freq_diff_list.append(freq_diff)
                                      
                            
                        freq_diff_list.sort()
                        x = freq_diff_list[0] 
                        for i in codon_list:
                            
                            if math.fabs(f[1] - hg19_bias[i][1]) == x and i != 'AGG':
                               

                                if sys.argv[1].find('U') == -1:
                                    i = i.replace('U', 'T')
                                elif sys.argv[1].find('T') == -1:
                                    i = i

                                print 'Now swapping ' + t + ' for ' + i + ' because its frequencies are close and both produce ' + f[0]
                                
                                #So it is easier to read in the seqeunce
                                i = i.lower()

                                dna2 = list(dna1)
                                
                                #replaces old codon overlapping PAM with new similar codon
                                dna2[a.end()], dna2[a.end() + 1], dna2[a.end() + 2] = i
                                    
                                #list becomes a string again for legibility 
                                dna2 = ''.join(dna2)
                                
                                #print out the new DNA string
                                if sys.argv[1].find('U') == -1:
                                    dna2 = dna2.replace('U', 'T')
                                    print dna2
                                    quit()
                                else:
                                    print dna2   
                                    quit()
                   
            #This asks whether the end position is in a different frame based on the positions remainder
            elif (a.end()%3) == 1:
                        
                #everything after this point is the same as above for a different frame or checking the reverse complement of the target
                codon_PAM = dna1[a.end() - 1], dna1[a.end()], dna1[a.end() + 1]

                codon_PAM = ''.join(codon_PAM)
                        
                print'The PAM is in a 231 frame'
                
                if sys.argv[1].find('T') == -1:
                    t = codon_PAM.replace('T', 'U')
                elif sys.argv[1].find('U') == -1:
                    t = codon_PAM
                      
                print t + ' Is the overlapping codon'

                print t +' produces ' + hg19_bias[codon_PAM][0]
                        
                codon_list = []
                for i in hg19_bias:
                    if hg19_bias[i][0]==hg19_bias[codon_PAM][0]:
                        codon_list.append(i)

                
                for b in codon_list:
                    if b == t:
                        codon_list.remove(b)
              
                        f = hg19_bias[codon_PAM]
                        del hg19_bias[codon_PAM]
                        

                        freq_diff = 0
                        freq_diff_list = [] 
                    
                        for i in codon_list:
                            

                            freq_diff = f[1] - hg19_bias[i][1]
                            freq_diff = math.fabs(freq_diff)
                            freq_diff_list.append(freq_diff)
                                      
                            
                        freq_diff_list.sort()
                        x = freq_diff_list[0] 
                        for i in codon_list:
                            
                            if math.fabs(f[1] - hg19_bias[i][1]) == x:

                                if sys.argv[1].find('U') == -1:
                                    i = i.replace('U', 'T')
                                elif sys.argv[1].find('T') == -1:
                                    i = i
                            
                                print 'Now swapping ' + t + ' for ' + i + ' because its frequencies are close and both produce ' + f[0]
                                
                                i = i.lower()

                                dna2 = list(dna1)

                                dna2[a.end() - 1], dna2[a.end()], dna2[a.end() + 1] = i
                                
                                dna2 = ''.join(dna2)
                            
                                if sys.argv[1].find('U') == -1:
                                    dna2 = dna2.replace('U', 'T')
                                print dna2
                
            else:
                print 'PAM is in 321 frame. You cannot alter the PAM without altering the future peptide'
        else:
            print 'Target not adjacent to a PAM site'

elif dna1.count(gRNA) == 0 and revcom_dna1.count(rev_gRNA) == 1:
    
    print 'Reversing the guide and reverse-complementing the target:'
    
    if sys.argv[1].find('U') == -1:
        rD = revcom_dna1.replace('U', 'T')
        rG = rev_gRNA.replace('U', 'T')
    elif sys.argv[1].find('T') == -1:
        rD = revcom_dna1
        rG = rev_gRNA
        
    print 'This is the target:\n'+ rD
    print 'This is the guide:\n' + rG    

    dna1 = revcom_dna1
    gRNA = rev_gRNA
    
    for a in re.finditer(gRNA, dna1):    

        
        if dna1[a.end() + 1] == 'G' and dna1[a.end() + 2] == 'G':
            
            
            PAM = dna1[a.end()], dna1[a.end() + 1], dna1[a.end() + 2]

            PAM = ''.join(PAM)

            print 'The guide matches the sequence once at: '
            print (a.start(), a.end())
            if sys.argv[1].find('U') == -1:
                t = PAM.replace('U', 'T')
            elif sys.argv[1].find('T') == -1:
                t = PAM 
   
            print t + ' is the PAM sequence'
                
            if (a.end()%3) == 0:
                print'The PAM is in a 123 frame'
                
                print t + ' Is the overlapping codon'

                if hg19_bias[PAM][0] == 'W':
                    print 'Cannot alter tryptophan due to only one codon'
                else:

                    print t + ' produces ' + hg19_bias[PAM][0]
                        

                    codon_list = []
                    for i in hg19_bias:
                        if hg19_bias[i][0]==hg19_bias[PAM][0]:
                            codon_list.append(i)
                           
                    for b in codon_list:
                        if b == t:
                            codon_list.remove(b)
                                       
                            f = hg19_bias[PAM]
                            del hg19_bias[PAM]
                        

                            freq_diff = 0
                            freq_diff_list = [] 
                    
                            for i in codon_list:
                            

                                freq_diff = f[1] - hg19_bias[i][1]
                                freq_diff = math.fabs(freq_diff)
                                freq_diff_list.append(freq_diff)
                                      
                            
                            freq_diff_list.sort()
                            x = freq_diff_list[0] 
                            for i in codon_list:
                            
                                if math.fabs(f[1] - hg19_bias[i][1]) == x and i != 'AGG':

                                    if sys.argv[1].find('U') == -1:
                                        i = i.replace('U', 'T')
                                    elif sys.argv[1].find('T') == -1:
                                        i = i

                                    print 'Now swapping ' + t + ' for ' + i + ' because its frequencies are close and both produce ' + f[0]
                                        
                                    i = i.lower()

                                    dna2 = list(dna1)
                                            
                                    dna2[a.end()], dna2[a.end() + 1], dna2[a.end() + 2] = i
                                    
                                    dna2 = ''.join(dna2)
                                    
                                    if sys.argv[1].find('U') == -1:
                                        dna2 = dna2.replace('U', 'T')
                                        print dna2
                                    else:
                                        print dna2
                
            elif (a.end()%3) == 1:
                        
                codon_PAM = dna1[a.end() - 1], dna1[a.end()], dna1[a.end() + 1]

                codon_PAM = ''.join(codon_PAM)
                        
                print'The PAM is in a 231 frame'
                
                if sys.argv[1].find('T') == -1:
                    t = codon_PAM.replace('T', 'U')
                elif sys.argv[1].find('U') == -1:
                    t = codon_PAM
                      
                print t + ' Is the overlapping codon'

                print t +' produces ' + hg19_bias[codon_PAM][0]
                        
                codon_list = []
                for i in hg19_bias:
                    if hg19_bias[i][0]==hg19_bias[codon_PAM][0]:
                        codon_list.append(i)
                       
                for b in codon_list:
                    if b == t:
                        codon_list.remove(b)
              
                        f = hg19_bias[codon_PAM]
                        del hg19_bias[codon_PAM]
                        

                        freq_diff = 0
                        freq_diff_list = [] 
                    
                        for i in codon_list:
                            

                            freq_diff = f[1] - hg19_bias[i][1]
                            freq_diff = math.fabs(freq_diff)
                            freq_diff_list.append(freq_diff)
                                      
                            
                        freq_diff_list.sort()
                        x = freq_diff_list[0] 
                        for i in codon_list:
                            
                            if math.fabs(f[1] - hg19_bias[i][1]) == x:

                                if sys.argv[1].find('U') == -1:
                                    i = i.replace('U', 'T')
                                elif sys.argv[1].find('T') == -1:
                                    i = i
                            
                                print 'Now swapping ' + t + ' for ' + i + ' because its frequencies are close and both produce ' + f[0]
                                
                                i = i.lower()

                                dna2 = list(dna1)
                                    
                                dna2[a.end() - 1], dna2[a.end()], dna2[a.end() + 1] = i
                                                    
                                dna2 = ''.join(dna2)
                            
                                if sys.argv[1].find('U') == -1:
                                    dna2 = dna2.replace('U', 'T')
                                print dna2
                
            else:
                print 'PAM is in 321 frame. You cannot alter the PAM without altering the future peptide'
        else:
            print 'Target not adjacent to a PAM site'
else:
    print 'Guide does not match or is not unique'
