'''
EG and AZF 2024 DragonRNA
This code will identify the most common DragonRNA sequences starting with the DNA input
and extending with sequences complementary to the desired length

Syntax: python3 DragonMostCommon.py <input_fastq_file> <input_DNA_oligo_sequence> 
Option to get the output saved in a file: python3 DragonMostCommon.py <input_fastq_file> <input_DNA_oligo_sequence> > <output_name>

Note: The output sequence ID is hardcoded at the end and 
    an if statement at the end could have a clause added to specify the length of the sequence read
    which can be changed or removed (we tested out of all the most common without length requirements and the expected lengths based on gel assay)
'''

import sys
from collections import Counter

input_file=sys.argv[1]
template=sys.argv[2]

def complement(s):
    '''This function will provide the complement of a given sequence; 
    to get the antisense of a given sequence, run antisense(sequence)[::-1]
    '''
    return s.upper().replace('A','t').replace('T','a').replace('G','c').replace('C','g').upper()

F = open(input_file,mode='rt').read().splitlines()

Sequences1 = F[1::4]
StartDNA1 = str(template)
StartDNAComplement1 = antisense(StartDNA1)[::-1]

klen1 = 10
tlen1 = len(StartDNAComplement1)

WordsInStartDNA1 = Counter([StartDNA1[i:i+klen1] for i in range(tlen1-klen1+1)])
WordsInStartComplement1 = Counter([StartDNAComplement1[i:i+klen1] for i in range(tlen1-klen1+1)])
DNAWords1 = WordsInStartDNA1+WordsInStartComplement1

IlluminaTetritis1 = open('<path_to>/illuminatetritis1223.fa', mode='rt').read()
IlluminaTetritis1+=complement(IlluminaTetritis1)
WordsInIllunimaTetritis1 = Counter([IlluminaTetritis1[i:i+klen1] for i in range(len(IlluminaTetritis1)-klen1+1)])

rawCounter1 = Counter(Sequences1)
filteredCounter1 = Counter()

for seq1 in rawCounter1:
    OM1 = False ## Oligo  Match
    IM1 = False ## Illumina  Match
    for i1 in range(len(seq1)-klen1+1):
        w1 = seq1[i1:i1+klen1]
        if w1 in DNAWords1:
            OM1 = True
        if w1 in WordsInIllunimaTetritis1:
            IM1 = True
    #have the option here to edit the length desired for the sequences depending on the size of template and the size of expected extension
    #this if statement could have a size range added to restrict output to the most common sequence reads of a desired size range based on gel assays
    #syntax of that if statement would be: if OM1 and not(IM1) and 60<=len(seq1)<=80:
    if OM1 and not(IM1):
        filteredCounter1[seq1] = rawCounter1[seq1]

#output will be to the terminal or you can end the syntax to save the output in a file: > <output file> 
for i1,(s1,c1) in enumerate(filteredCounter1.most_common(10)):
    print('>hmtRNAP_vlongtemp-0h_MostCommon_'+str(i1)+'_Counts:'+str(c1))
    print(s1)      
