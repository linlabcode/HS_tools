#This is the Cas9 Blocker! 

#Determine sequence match
#Locate pam site frame
#Alter the overlapping PAM sequence without changing the codon-AA

import re                             
import collections
import sys
sys.path.append('/home/harrison/Downloads/pipeline-master')
import os
import utils

from collections import defaultdict

table=utils.parseTable('CodonUsageTables.txt', '\t')

hg19_bias = defaultdict()
for line in table[1:]:
    hg19_bias[line[0]]=[line[1],float(line[2])]



#This checks if the input is a text
if os.path.isfile(sys.argv[1]) and os.path.isfile(sys.argv[2]):
#This opens the input files
    with open(sys.argv[1]) as f, open(sys.argv[2]) as g:
#This reads the file and makes the text a usuable string                          
        dna_input = f.read().rstrip()                                                
        gRNA_input = g.read().rstrip()

#If the input was not a file than continue with the code
else:                                                                                
    dna_input = sys.argv[1]
    gRNA_input = sys.argv[2]  
    pass



#This will replace all 'T's with 'U's if DNA was inputted otherwise    
dna_input = dna_input.replace('T', 'U')                                                
gRNA_input = gRNA_input.replace('T', 'U')
dna1_input = dna_input.upper()
gRNA_input = gRNA_input.upper()

#This will make the input a variable that the code uses.
dna1 = dna_input                                                                    
gRNA = gRNA_input 

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
            if sys.argv[1].find('U') == -1:
                t = PAM.replace('U', 'T')
            elif sys.argv[1].find('T') == -1:
                t = PAM    
            print t + ' is the PAM sequence'
                
            #the remainder indirectly tells us the postion of the first base in the PAM site which tells us what frame the codon is in
            if (a.end()%3) == 0:
                print'The PAM is in a 123 frame'
                
                print t + ' Is the overlapping codon'
                if hg19_bias[PAM][0] == 'W':
                    print 'Cannot alter tryptophan due to only one codon'
                    
                print t + ' produces ' + hg19_bias[PAM][0]
                    

                codon_list = []
                for i in hg19_bias:
                    if hg19_bias[i][0]==hg19_bias[PAM][0]:
                        codon_list.append(i)
                       
                freq = 0
                if i != PAM:
                    for i in codon_list:
                        if freq == 0:
                            freq = hg19_bias[i][1]
                            best_codon = i
                        if hg19_bias[i][1] > freq:
                            freq = hg19_bias[i][1]

                    if sys.argv[1].find('U') == -1:
                        i = i.replace('U', 'T')
                    elif sys.argv[1].find('T') == -1:
                        i = i

                    print 'Now swapping ' + t + ' for ' + i + ' because it has the highest bias for ' + hg19_bias[PAM][0]
                        
                    i = i.lower()

                    dna2 = list(dna1)
                            #replaces the orginal codon in the sequence to a new desirable codon that produces the same amino acid with highest bias
                    dna2[a.end()], dna2[a.end() + 1], dna2[a.end() + 2] = i
                        
                            #list becomes a string again for legibility 
                    dna2 = ''.join(dna2)
                    
                    
                    if sys.argv[1].find('U') == -1:
                        dna2 = dna2.replace('U', 'T')
                        print dna2
                    else:
                        print dna2
                

                   #This asks whether the end position is in a different frame based on the positions remainder
            elif (a.end()%3) == 1:
                        
                        #everything after this point is the same as above for a different frame or checking the reverse complement of the target
                codon_PAM = dna1[a.end() - 1], dna1[a.end()], dna1[a.end() + 1]

                codon_PAM = ''.join(codon_PAM)
                        
                print'The PAM is in a 231 frame'
                
                if sys.argv[1].find('U') == -1:
                    t = codon_PAM.replace('U', 'T')
                elif sys.argv[1].find('T') == -1:
                    t = codon_PAM
                      
                print t + ' Is the overlapping codon'

                print t +' produces ' + hg19_bias[codon_PAM][0]
                        
                codon_list = []
                for i in hg19_bias:
                    if hg19_bias[i][0]==hg19_bias[codon_PAM][0]:
                        codon_list.append(i)
                       
                freq = 0
                
                    
                for i in codon_list:
                    if i != codon_PAM:
                        
                        if freq == 0:
                            freq = hg19_bias[i][1]
                            best_codon = i
                        if hg19_bias[i][1] > freq:
                            freq = hg19_bias[i][1]

                        if sys.argv[1].find('U') == -1:
                            i = i.replace('U', 'T')
                        elif sys.argv[1].find('T') == -1:
                            i = i
                    
                        print 'Now swapping ' + t + ' for ' + i + ' because it has the next highest bias for ' + hg19_bias[codon_PAM][0]
                        
                        i = i.lower()

                        dna2 = list(dna1)
                            #replaces the orginal codon in the sequence to a new desirable codon that produces the same amino acid with highest bias
                        dna2[a.end() - 1], dna2[a.end()], dna2[a.end() + 1] = i
                        
                            #list becomes a string again for legibility 
                        dna2 = ''.join(dna2)
                    
                        #print 'Now swapping ' + PAM + ' for ' + i + ' because it has the highest bias for ' + hg19_bias[PAM][0]
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
    #Finds the location of the match

    dna1 = revcom_dna1
    gRNA = rev_gRNA
    
    for a in re.finditer(gRNA, dna1):    

        #This makes sure the adjacent 3 bases make a PAM site
        if dna1[a.end() + 1] == 'G' and dna1[a.end() + 2] == 'G':
            
            #This makes the PAM site a PAM variable
            PAM = dna1[a.end()], dna1[a.end() + 1], dna1[a.end() + 2]

            PAM = ''.join(PAM)

            print 'The guide matches the sequence once at: '
            print (a.start(), a.end())
            if sys.argv[1].find('U') == -1:
                t = PAM.replace('U', 'T')
            elif sys.argv[1].find('T') == -1:
                t = PAM 
   
            print t + ' is the PAM sequence'
                
            #the remainder indirectly tells us the postion of the first base in the PAM site which tells us what frame the codon is in
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
                           
                    freq = 0
                    if i != PAM:
                        for i in codon_list:
                            if freq == 0:
                                freq = hg19_bias[i][1]
                                best_codon = i
                            if hg19_bias[i][1] > freq:
                                freq = hg19_bias[i][1]

                        if sys.argv[1].find('U') == -1:
                            i = i.replace('U', 'T')
                        elif sys.argv[1].find('T') == -1:
                            i = i

                        print 'Now swapping ' + t + ' for ' + i + ' because it has the highest bias for ' + hg19_bias[PAM][0]
                            
                        i = i.lower()

                        dna2 = list(dna1)
                                #replaces the orginal codon in the sequence to a new desirable codon that produces the same amino acid with highest bias
                        dna2[a.end()], dna2[a.end() + 1], dna2[a.end() + 2] = i
                            
                                #list becomes a string again for legibility 
                        dna2 = ''.join(dna2)
                        
                        
                        if sys.argv[1].find('U') == -1:
                            dna2 = dna2.replace('U', 'T')
                            print dna2
                        else:
                            print dna2
                

                   #This asks whether the end position is in a different frame based on the positions remainder
            elif (a.end()%3) == 1:
                        
                        #everything after this point is the same as above for a different frame or checking the reverse complement of the target
                codon_PAM = dna1[a.end() - 1], dna1[a.end()], dna1[a.end() + 1]

                codon_PAM = ''.join(codon_PAM)
                        
                print'The PAM is in a 231 frame'
                
                if sys.argv[1].find('U') == -1:
                    t = codon_PAM.replace('U', 'T')
                elif sys.argv[1].find('T') == -1:
                    t = codon_PAM
                      
                print t + ' Is the overlapping codon'

                print t +' produces ' + hg19_bias[codon_PAM][0]
                        
                codon_list = []
                for i in hg19_bias:
                    if hg19_bias[i][0]==hg19_bias[codon_PAM][0]:
                        codon_list.append(i)
                       
                freq = 0
                if i != codon_PAM:
                    for i in codon_list:
                        if freq == 0:
                            freq = hg19_bias[i][1]
                            best_codon = i
                        if hg19_bias[i][1] > freq:
                            freq = hg19_bias[i][1]

                    if sys.argv[1].find('U') == -1:
                        i = i.replace('U', 'T')
                    elif sys.argv[1].find('T') == -1:
                        i = i
                
                    print 'Now swapping ' + t + ' for ' + i + ' because it has the highest bias for ' + hg19_bias[codon_PAM][0]
                    
                    i = i.lower()

                    dna2 = list(dna1)
                        #replaces the orginal codon in the sequence to a new desirable codon that produces the same amino acid with highest bias
                    dna2[a.end() - 1], dna2[a.end()], dna2[a.end() + 1] = i
                    
                        #list becomes a string again for legibility 
                    dna2 = ''.join(dna2)
                
                    #print 'Now swapping ' + PAM + ' for ' + i + ' because it has the highest bias for ' + hg19_bias[PAM][0]
                    if sys.argv[1].find('U') == -1:
                        dna2 = dna2.replace('U', 'T')
                    print dna2
                
            else:
                print 'PAM is in 321 frame. You cannot alter the PAM without altering the future peptide'
        else:
            print 'Target not adjacent to a PAM site'
else:
    print 'Guide does not match or is not unique'


#     print 'Reversing the guide and reverse-complementing the target:'

#     if sys.argv[1].find('U') == -1:
#         revcom_dna1 = revcom_dna1.replace('U', 'T')
#         rev_gRNA = rev_gRNA.replace('U', 'T')
            
#     print 'This is the target:\n'+ revcom_dna1
#     print 'This is the guide:\n' + rev_gRNA      
    
#     for a in re.finditer(rev_gRNA, revcom_dna1):    
#         dna1 = revcom_dna1
#         gRNA = rev_gRNA
#         #This makes sure the adjacent 3 bases make a PAM site
#         if dna1[a.end() + 1] == 'G' and dna1[a.end() + 2] == 'G':
            
#             #This makes the PAM site a PAM variable
#             PAM = dna1[a.end()], dna1[a.end() + 1], dna1[a.end() + 2]

#             PAM = ''.join(PAM)

#             print 'The guide matches the sequence once at: '
#             print (a.start(), a.end())
#             print PAM + ' is the PAM sequence'
#             print hg19_bias[PAM]
#             #the remainder indirectly tells us the postion of the first base in the PAM site which tells us what frame the codon is in
#             if (a.end()%3) == 0:
#                 print'The PAM is in a 123 frame'
#                 print PAM +' produces ' + hg19_bias[PAM][0]

#                 codon_list = []
#                 for i in hg19_bias:
#                     if hg19_bias[i][0]==hg19_bias[PAM][0]:
#                         codon_list.append(i)
                   
#                 freq = 0
#                 for i in codon_list:
#                     if freq == 0:
#                         freq = hg19_bias[i][1]
#                         best_codon = i
#                     if hg19_bias[i][1] > freq:
#                         freq = hg19_bias[i][1]
                
#                 print 'Now swapping ' + PAM + ' for ' + i + ' because it has the highest bias for ' + hg19_bias[PAM][0]
                
#                 i = i.lower()

#                 dna2 = list(dna1)
#                     #replaces the orginal codon in the sequence to a new desirable codon that produces the same amino acid with highest bias
#                 dna2[a.end()], dna2[a.end() + 1], dna2[a.end() + 2] = i
                
#                     #list becomes a string again for legibility 
#                 dna2 = ''.join(dna2)
            
#                 #print 'Now swapping ' + PAM + ' for ' + i + ' because it has the highest bias for ' + hg19_bias[PAM][0]
#                 print dna2


#             #This asks whether the end position is in a different frame based on the positions remainder
#             elif (a.end()%3) == 1:
                
#                 #everything after this point is the same as above for a different frame or checking the reverse complement of the target
#                 codon_PAM = dna1[a.end() - 1], dna1[a.end()], dna1[a.end() + 1]

#                 codon_PAM = ''.join(codon_PAM)
                
#                 print'The PAM is in a 231 frame'
#                 print codon_PAM + ' Is the overlapping codon'
#                 print codon_PAM +' produces ' + hg19_bias[codon_PAM][0]

#                 codon_list = []
#                 for i in hg19_bias:
#                     if hg19_bias[i][0]==hg19_bias[codon_PAM][0]:
#                         codon_list.append(i)
                   
#                 freq = 0
#                 for i in codon_list:
                
#                     if freq == 0:
#                         freq = hg19_bias[i][1]
#                         best_codon = i
#                     if hg19_bias[i][1] > freq:
#                         freq = hg19_bias[i][1]
            
#                 print 'Now swapping ' + codon_PAM + ' for ' + i + ' because it has the highest bias for ' + hg19_bias[codon_PAM][0]
                
#                 i = i.lower()

#                 dna2 = list(dna1)
#                     #replaces the orginal codon in the sequence to a new desirable codon that produces the same amino acid with highest bias
#                 dna2[a.end() -1], dna2[a.end()], dna2[a.end() + 1] = i
                
#                     #list becomes a string again for legibility 
#                 dna2 = ''.join(dna2)
            
#                 #print 'Now swapping ' + PAM + ' for ' + i + ' because it has the highest bias for ' + hg19_bias[PAM][0]
                
#                 print dna2
            
#             else:
#                 print 'PAM is in 321 frame. You cannot alter the PAM without altering the future peptide'
#         else:
#             print 'Target not adjacent to a PAM site'
# else:
#     print 'Guide does not match'
