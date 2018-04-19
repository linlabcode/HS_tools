#This is the Cas9 Blocker! 

#Determine sequence match
#Locate pam site frame
#Alter the overlapping PAM sequence without changing the codon-AA

import re                             
import collections
import random
import sys
import os

#This checks if the input is a text
if os.path.isfile(sys.argv[1]) and os.path.isfile(sys.argv[2]):
#This opens the input files
    with open(sys.argv[1]) as f, open(sys.argv[2]) as g:
#This reads the file and makes the text a usuable string                          
        dna_input = f.read().rstrip()                                                
        gRNA_input = g.read().rstrip()

#If the input was not a file than continue witht the code
else:                                                                                
    dna_input = sys.argv[1]
    gRNA_input = sys.argv[2]  
    pass

print 'This is the target:\n'+ dna_input
print 'This is the guide:\n' + gRNA_input  

#This will replace all 'T's with 'U's if DNA was inputted otherwise    
dna_input = dna_input.replace('T', 'U')                                                
gRNA_input = gRNA_input.replace('T', 'U')
dna1_input = dna_input.upper()
gRNA_input = gRNA_input.upper()

#This will make the input a variable that the code uses.
dna1 = dna_input                                                                    
gRNA = gRNA_input 

codon_dict = {
    'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
    'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
    'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
    'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
    'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
    'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
    'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
    'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
    'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_',
    'UGC':'C', 'UGU':'C', 'UGA':'_', 'UGG':'W'}

#this makes the table into a grouped list and swaps the key and value from the previous table
codon_dict2 = collections.defaultdict(list)                                         
for k, v in codon_dict.iteritems():
    codon_dict2[v].append(k)



def reverse_complement(dna):
                complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A',}
                return ''.join([complement[base] for base in dna[::-1]])

#this gives you a reverse version of the target and guide
rev_gRNA = gRNA[::-1]

revcom_dna1 = reverse_complement(dna1)


#Checks if the match is unique in the original sequence
if dna1.count(gRNA) == 1 and revcom_dna1.count(rev_gRNA) == 0:
        
    #Finds the location of the match
    for a in re.finditer(gRNA, dna1):    

        #This makes sure the adjacent 3 bases make a PAM site
        if dna1[a.end() + 1] == 'G' and dna1[a.end() + 2] == 'G':
            
            #This makes the PAM site a PAM variable
            PAM = dna1[a.end()], dna1[a.end() + 1], dna1[a.end() + 2]

            PAM = ''.join(PAM)

            print 'The guide matches the sequence once at: '
            print (a.start(), a.end())
            print PAM + ' is the PAM sequence'

            #the remainder indirectly tells us the postion of the first base in the PAM site which tells us what frame the codon is in
            if (a.end()%3) == 0:
                print'The PAM is in a 123 frame'
                print PAM +' produces ' + codon_dict[PAM]
                
                #variable for value of the PAM codon key 
                aa = codon_dict[PAM]
                #variable for the values that refer to the amino acid key found based on the codon given
                codon_list = codon_dict2[aa]

                #variable for selection of a codon from the given amino acid key
                codon = (random.choice(codon_list))

                b = codon

                #This selects for a new codon in the list group and avoids undesirable replacements
                if b != PAM and b != 'AGG' and b != 'CGG':
                    print PAM + ' has been changed to ' + b
                    
                    #So its easier to find on the new sequence
                    b = b.lower()    
                    #Makes the string a list for manipulation
                    dna2 = list(dna1)
                    #replaces the orginal codon in the sequence to a new desirable codon that produces the same amino acid as the previous one
                    dna2[a.end()], dna2[a.end() + 1], dna2[a.end() + 2] = b
                    #list becomes a string again for legibility 
                    dna2 = ''.join(dna2)
                    
                    #if DNA was given then output DNA as well otherwise RNA will remain
                    if dna_input.find('T') == -1:
                        dna2 = dna2.replace('U', 'T')
                    
                    print 'The new dna sequence; ' + dna2 + ', and contains no PAM site'

                else:
                    print 'Try again'
            
            #This asks whether the end position is in a different frame based on the positions remainder
            elif (a.end()%3) == 1:
                
                #everything after this point is the same as above for a different frame or checking the reverse complement of the target
                nearby_codon = dna1[a.end() - 1], dna1[a.end()], dna1[a.end() + 1]

                nearby_codon = ''.join(nearby_codon)
                
                print'The PAM is in a 231 frame'
                print nearby_codon + ' overlaps with PAM and produces ' + codon_dict[nearby_codon]
                
                aa = codon_dict[nearby_codon]
                codon_list = codon_dict2[aa]
                
                codon = (random.choice(codon_list))

                b = codon 

                if b != nearby_codon and b != 'AGG' and b != 'UGG':
                    print nearby_codon + ' has been changed to ' + b
                    
                    b = b.lower()
                    dna2 = list(dna1)
                    dna2[a.end() - 1], dna2[a.end()], dna2[a.end() + 1] = b
                    dna2 = ''.join(dna2)

                    if dna_input.find('T') == -1:
                        dna2 = dna2.replace('U', 'T')

                    print 'The new dna sequence; ' + dna2 + ', and contains no PAM site'

                else:
                    print 'Try again'
            else:
                print 'PAM is in a 312 frame. \nAltering PAM would cause a change in the peptide produced.' 
        else:
            print 'Target not adjacent to a PAM site'
            
#Same as everything above but for the opposite strand of DNA/RNA given
elif dna1.count(gRNA) == 0 and revcom_dna1.count(rev_gRNA) == 1:                
    print 'Reverse compliment match'
    for a in re.finditer(rev_gRNA, revcom_dna1):    

        if revcom_dna1[a.end() + 1] == 'G' and revcom_dna1[a.end() + 2] == 'G':
            
            PAM = revcom_dna1[a.end()], revcom_dna1[a.end() + 1], revcom_dna1[a.end() + 2]

            PAM = ''.join(PAM)

            print 'We have a pam site at: '
            print (a.start(), a.end())
            print PAM

            if (a.end()%3) == 0:
                print'The PAM is in a 123 frame'
                print PAM +' produces ' + codon_dict[PAM]
                
                aa = codon_dict[PAM]
                codon_list = codon_dict2[aa]
                
                codon = (random.choice(codon_list))

                b = codon 

                if b != PAM and b != 'AGG' and b != 'CGG':
                    print PAM + ' has been changed to ' + b
                    
                    b = b.lower()
                    dna2 = list(revcom_dna1)
                    dna2[a.end()], dna2[a.end() + 1], dna2[a.end() + 2] = b
                    dna2 = ''.join(dna2)

                    if dna_input.find('T') == -1:
                        dna2 = dna2.replace('U', 'T')                        

                    print 'The new dna sequence; ' + dna2 + ', and contains no PAM site'

                else:
                    print 'Try again'
            
            elif (a.end()%3) == 1:
                
                nearby_codon = revcom_dna1[a.end() - 1], revcom_dna1[a.end()], revcom_dna1[a.end() + 1]

                nearby_codon = ''.join(nearby_codon)
                
                print'The PAM is in a 231 frame'
                print nearby_codon + ' overlaps with PAM and produces ' + codon_dict[nearby_codon]
                
                aa = codon_dict[nearby_codon]
                codon_list = codon_dict2[aa]
                
                codon = (random.choice(codon_list))

                b = codon 

                if b != nearby_codon and b != 'AGG' and b != 'UGG' and 'CGG':
                    print nearby_codon + ' has been changed to ' + b
                    
                    b = b.lower()
                    dna2 = list(revcom_dna1)
                    dna2[a.end() - 1], dna2[a.end()], dna2[a.end() + 1] = b
                    dna2 = ''.join(dna2)
                    
                    if dna_input.find('T') == -1:
                        dna2 = dna2.replace('U', 'T')
                        
                    print 'The new dna sequence; ' + dna2 + ', and contains no PAM site'

                else:
                    print 'Try again'
            else:
                print 'PAM is in a 312 frame. \nAltering PAM would cause a change in the peptide produced.' 
        else:
            print 'Target not adjacent to a PAM site'
    else:
        print('Target not unique')
else:
    print('Target not unique')  