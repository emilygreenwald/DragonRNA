'''
Emily Greenwald
Sept 2024
this script will identify sequences that match an input sequence with extension using the reverse complement of the primer as a template
and count how frequently they appear in a fastq file. the input contains a region of unspecified bases "NN"s
'''

from collections import Counter
import sys
import re
import gzip
import argparse
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import seaborn as sns
import pandas as pd 

LinkerStart='TGGAATTCTC'

def complement(s):
    '''this function will give the complement version of an inputed sequence (note NOT reverse complement, only complement)
    '''
    return s.upper().replace('A', 't').replace('T', 'a').replace('G', 'c').replace('C', 'g').upper()

def generate_primers(input_file):
    '''this will generate a list of primers given that there are two unknown bases
    using the string name of the input file which indicates which primer was used
    '''
    bases = 'ATCG'
    
    doublets=[]
    primer_set_73=[]
    primer_set_74=[]

    for base1 in bases:
        for base2 in bases:
            doublets.append(base1 + base2)
    
    for doublet in doublets:
        primer_73='AGA'+doublet+'ATTATTACGTGCTTTTGTTCAA'
        primer_set_73.append(primer_73)
        primer_74='TT'+doublet+'ACGTCAACGATATAAGTTTTGAC'
        primer_set_74.append(primer_74)
    
    if '73' in input_file:
        return primer_set_73
    elif '74' in input_file:
        return primer_set_74
    elif 'nn' in input_file:
        return doublets
    else:
        print("Oligo number 73 and 74 not found in input_file name.")
    

def search_and_count_sequences(F, primers):
    '''
    this function will search sequences in the FASTQ file data that match the primer with some extension
    and count the occurrences of each matching sequence.

    F=List of lines from the FASTQ file.
    primer_sequence=The primer sequence to search for, hardcoded into "main" because there's only two options for this NN analysis
    return=A Counter object with sequence counts, total number of reads with the primer at the 5' end at all
    '''
    #read only every 4th line because that contains sequence information
    Sequences1 = F[1::4]
    # print(f"First few sequences: {Sequences1[:5]}") #works identifies sequences

    #counter object to keep track of each primer: counter object of each sequence read from that primer
    primer_counters = {primer: Counter() for primer in primers}
    #dictionary of lists to keep track of primers and sequences
    # primer_reads_dict = {primer:[] for primer in primers}

    #dictionaries to keep track of how many total counts of reads there are for each primer and for each extension product of minimum length
    primer_total_counts = {primer: 0 for primer in primers}
    primer_min_length_counts = {primer: 0 for primer in primers}

    #need minimum length of primer sequence
    min_length = len(primers[0])

    #have to subtract last 4 bases of primer for template portion because they are the purported tetraloop
    # template_sequence=complement(primer_sequence)[::-1][len(primer_sequence)-4:] #this was for if have a single primer sequence
    #loop through each sequence
    for seq in Sequences1: 
        #loop through each primer option 
        for primer_sequence in primers:
            
            if primer_sequence in seq and LinkerStart in seq:
                pos=seq.rfind(LinkerStart)
                seq_trim=seq[:pos]

                primer_total_counts[primer_sequence] += 1

                if len(seq_trim) >= min_length:
                    # primer_reads_dict[primer_sequence].append(seq)
                    primer_min_length_counts[primer_sequence] += 1

                    #add to the counter object for the given sequence read for the given primer
                    primer_counters[primer_sequence][seq_trim] += 1
                    #if this primer is in the sequence, move to next sequence 
                    break

    # print(primer_reads_dict)
    #also save total counts across primers
    total_with_primer = sum(primer_total_counts.values())
    total_with_primer_of_min_length = sum(primer_min_length_counts.values())

    print('total counts with primer', total_with_primer)
    print('total counts with primer of minimum length of slightly less than double primer length: ', total_with_primer_of_min_length)
    return primer_counters, primer_total_counts, primer_min_length_counts, total_with_primer, total_with_primer_of_min_length


def kmer_count(sequence):
    '''this function will analyze the input oligo sequence
     for kmer sequences of length 10 and return 
    a Counter object of k-mers from that sequence.'''
    klen1 = 10

    return Counter([sequence[i:i+klen1] for i in range(len(sequence) - klen1 + 1)])


def analyze_nn_complementarity(primer_counters, primers, nn_start_pos):
    '''this function will analyze the counter object generated in search_and_count_sequences
    for whether the portions that are extended were extended based on the same molecule as 
    the start of the read, or a different molecule (testing for inter vs intra molecular priming and 
    templating). this function will use kmer counting to test for the same sequence vs different sequences.
    
    inputs are a dictionary of counter objects for primers to sequence reads and the counts of the reads
    and list of primer options and the position where the unspecified bases start
    
    output is a dictionary of counts of whether, for each input primer sequence, the extension is 
    complementary to the primer portion of the read or to a different primer'''

    #dictionary of matches vs non matches and total sequences with counts
    nn_complementarity_counts={}

    #dictionary of {primer_nn_sequence:{extension_sequence:count,etc},etc}
    nn_regions_counts={}

    doublet_options=generate_primers('nn')

    #keep track of whether extension NN is antisense or not to primer-region NN in same read
    total_match = 0
    total_not_match = 0

    #make key in dictionary of NN region for all options for NNs
    for nn in doublet_options:
        nn_regions_counts[nn]=Counter()

        for unspec in doublet_options:
            nn_regions_counts[nn][unspec]=0

    for primer_sequence in primers:
        # nn_complementary_count = 0
        # nn_non_complementary_count = 0
        # nn_match_details = {'Match':Counter(), 'Non-match':Counter()}

        sequences = primer_counters[primer_sequence]

        #go through counter object's sequences one sequence at a time
        for seq in sequences:

            primer_start = seq.find(primer_sequence)
            
            if primer_start == -1:
                continue
            elif primer_start != 0:
                continue
            
            primer_region = seq[primer_start:primer_start + len(primer_sequence)]

            extension_region = seq[primer_start + len(primer_sequence):]

            #determine the region with unspecified bases using the start position of NN (will be different for diff oligo inputs)
            nn_region = primer_region[nn_start_pos:nn_start_pos+2]
            
            known_region_before=primer_sequence[0 : nn_start_pos]
            known_region_after = primer_sequence[nn_start_pos+2:nn_start_pos+7]

            #check if NN in extension region is reverse complement of primer NN region
            #but need to know where NN is in extension region so need to look for sequences nearby with reverse complement
            #search for NN region using region before and after NN
            for i in range(len(extension_region) - 8):
                # Extract the sequence segment to check
                segment = extension_region[i:]

                #check if known region is matches reverse complement of  the segment, we know we are in the right place
                if segment[0:5] == complement(known_region_after)[::-1] and segment[7:]==complement(known_region_before)[::-1]:

                    extended_nn_region = segment[5:7]
                    # print(extended_nn_region)
                    
                    nn_regions_counts[nn_region][extended_nn_region] += primer_counters[primer_sequence][seq]
                    # print(primer_counters[primer_sequence][seq])

                    if extended_nn_region == complement(nn_region)[::-1]:
                        total_match += primer_counters[primer_sequence][seq]


                    #     nn_complementary_count += 1
                    #     nn_match_details['Match'][seq] += 1
                    #     break

                    else:
                        total_not_match += primer_counters[primer_sequence][seq]

                    #     nn_non_complementary_count += 1
                    #     nn_match_details['Non-match'][seq] += 1

        # nn_complementarity_counts[primer_sequence] = {
        #     'complementary_count': nn_complementary_count,
        #     'non_complementary_count': nn_non_complementary_count,
        #     'match_details': nn_match_details
        # }
        #confirmed works from this output with agnialign and checking whether NN region is actually complementary or not

    return nn_regions_counts, total_match, total_not_match


def plot_nn_complementarity_heatmapseaborn(nn_regions_counts, nn_start_pos, sample, total_match, total_not_match, save_as_svg=False, save_as_png=False, filename="nn_complementarity_heatmap.jpg"):
    '''Generates a heatmap showing how frequently each NN primer sequence is extended to each NN extension sequence.'''
    
    # Get all possible "NN" combinations for primer and extension
    bases = 'ATCG'
    nn_combinations = [base1 + base2 for base1 in bases for base2 in bases]
    nn_revcom = [complement(combo)[::-1] for combo in nn_combinations]

    # Initialize DataFrame for heatmap
    heatmap_data = pd.DataFrame(0, index=nn_combinations, columns=nn_revcom)

    for primer_sequence, details in nn_regions_counts.items():
        for extension_seq, count in details.items():
            heatmap_data.loc[extension_seq, primer_sequence] += count

    # Create heatmap
    cmap = sns.color_palette("YlGnBu", as_cmap=True)
    cmap.set_under('white')

    plt.figure(figsize=(10, 8))
    sns.heatmap(heatmap_data, annot=True, cmap=cmap, fmt='g', linewidths=.5, cbar_kws={'label': 'Count'}, square=True, vmin=1, vmax=heatmap_data.max().max())

    # Add labels and title
    plt.title(sample)
    plt.xlabel('Unspecified bases in input DNA oligonucleotide sequence')
    plt.ylabel('Unspecified bases in extension sequence')

    # Add total match and not match text
    plt.text(-0.1, -0.1, f"Total match primer (cis priming): {total_match}, Total not match primer (trans priming): {total_not_match}", 
         fontsize=10, ha='left', va='center', transform=plt.gca().transAxes)

    # Save as SVG or PNG if requested
    if save_as_svg:
        plt.savefig(f"{filename}.svg", format='svg', bbox_inches='tight')
    elif save_as_png:
        plt.savefig(f"{filename}.png", format='png', bbox_inches='tight')
    else:
        plt.show()

def main():
    #take in .fastq.gz input file of sequences and primer sequence in command line
    parser = argparse.ArgumentParser(description="Named parameters in script")
    parser.add_argument('--input_file', required=True, help="path to the input fastq file")
    args = parser.parse_args()
    
    input_file = args.input_file

    primers=generate_primers(input_file)
    #print('primer set: ', primers) #confirmed works and gets correct sequence

    if '73' in input_file.split('/')[-1]:
        unspecified_pos=3
    elif '74' in input_file.split('/')[-1]:
        unspecified_pos=2

    with gzip.open(input_file, mode='rt') as f:
        F = f.read().splitlines()
    # print(f"Number of lines read: {len(F)}") #confirmed reads file

    #this object [0] is the dictionary of primers to sequence counts counter object per primer, 
    #[1] is the dictionary of total counts of reads per primer
    #[2] is the dictionary of total counts of reads of minimum length per primer 
    #[3] is the total counts of every read that contains any primer sequence from the input
    #[4] is the total counts of everything of min length with any primer sequence from the input
    primer_counters, primer_total_counts, primer_min_length_counts, total_with_primer, total_with_primer_of_min_length = search_and_count_sequences(F, primers)

    nn_dict, total_match, total_not_match = analyze_nn_complementarity(primer_counters, primers, unspecified_pos)

    sample = input_file.split('/')[-1].split('.')[0]
    sample_title=sample.split('_')[:3]
    sample_title = ' '.join(sample_title)
    sample_title=re.sub('min ', ' minute timepoint with ', sample_title)

    plot_nn_complementarity_heatmapseaborn(nn_dict, unspecified_pos, sample_title, total_match, total_not_match, save_as_svg=True, filename=f'NNcomplementarity_heatmap_{sample}_{datetime.now().strftime("%Y%m%d-%H%M%S")}.jpg')
    # plot_nn_complementarity_seaborn(nn_dict, unspecified_pos, sample, save_as_jpg=True, filename=f'NNcomplementarity_{sample}_{datetime.now().strftime("%Y%m%d-%H%M%S")}.jpg')

    #make output file with total number of sequences that start with the primer and top 10 most common sequences and their counts and proportion
    output_filename = f'NNcomplementarity_heatmap_{sample}_{datetime.now().strftime("%Y%m%d-%H%M%S")}.jpg'

    num_entries_to_print = 3
        
    print(f"Results drawn to {output_filename}")


    


if __name__ == "__main__":
    main()
