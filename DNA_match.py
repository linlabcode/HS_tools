import difflib
from difflib import SequenceMatcher

dna1 = 'GTGTGACTGGGTCAAAAAAAAA'
dna2 = 'GTGTGAGTTAATTCTAACTTAG'


#searches for the same substring within the two sequences and provides location and length

match = SequenceMatcher(None, dna1, dna2).find_longest_match(0, len(dna1), 0, len(dna2))

print match
print (dna1[match.a: match.a + match.size])



