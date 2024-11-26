import sys
import gzip

input_file=sys.argv[1]
input_DNA=sys.argv[2].upper()

path=input_file.split('.')[0]
sample=path.split('/')[-1]
output_filename='DragonMostCommon'+sample.split('.')[0]+'.fa'

#add function/output to give total number of reads that start with DNA primer 
#can add the proportion of counting reads as a proportion/percent/fraction of that DNA primed number as fraction/percent 

def antisense(s):
    return s.upper().replace('A','t').replace('T','a').replace('G','c').replace('C','g').upper()
from collections import Counter
F = gzip.open(input_file,mode='rt').read().splitlines()

Sequences1 = F[1::4]
StartDNA1 = input_DNA
StartDNAComplement1 = antisense(StartDNA1).upper()

klen1 = 10
tlen1 = len(StartDNAComplement1)

WordsInStartDNA1 = Counter([StartDNA1[i:i+klen1] for i in range(tlen1-klen1+1)])
WordsInStartComplement1 = Counter([StartDNAComplement1[i:i+klen1] for i in range(tlen1-klen1+1)])
DNAWords1 = WordsInStartDNA1+WordsInStartComplement1

IlluminaTetritis1 = open('/Users/emilygreenwald/Documents/0Stanford/Fire_lab/Polymerase/DragonRNA/mtRNAP_templated_experiment/sequencing_adapter_trimmed/AndyDropFiles_040924/illuminatetritis1223.fa', mode='rt').read()
IlluminaTetritis1+=antisense(IlluminaTetritis1)
WordsInIllunimaTetritis1 = Counter([IlluminaTetritis1[i:i+klen1] for i in range(len(IlluminaTetritis1)-klen1+1)])

rawCounter1 = Counter(Sequences1)
filteredCounter1 = Counter()
total_with_primer=0


for seq1 in rawCounter1:
    OM1 = False ## Oligo  Match
    IM1 = False ## Illumina  Match
     # Check if the sequence starts with input_DNA or its complement
    if seq1.startswith(StartDNA1) or seq1.startswith(StartDNAComplement1):
        OM1 = True
        total_with_primer += rawCounter1[seq1]

    for i1 in range(len(seq1)-klen1+1):
        w1 = seq1[i1:i1+klen1]
        if w1 in DNAWords1:
            OM1 = True
        if w1 in WordsInIllunimaTetritis1:
            IM1 = True
        
    if OM1 and not(IM1):
        filteredCounter1[seq1] = rawCounter1[seq1]

output = open(output_filename,mode='wt')

output.write('>total_counts_startwith_primer: ' + str(total_with_primer))
output.write('\n')


for i1,(s1,c1) in enumerate(filteredCounter1.most_common(10)):
    print('>total_counts_startwith_primer: ' + str(total_with_primer))
    print('>'+sample+' '+str(i1)+' Counts:'+str(c1))
    print(s1)
    output.write('>'+sample+' '+str(i1)+' Counts:'+str(c1)+' Proportion of reads that start with primer: '+str(c1/total_with_primer))
    output.write('\n')
    output.write(s1)
    output.write('\n')





