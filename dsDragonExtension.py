'''
Emily Greenwald
August 2024
this script will identify sequences that match an input primer and optional template (or use the reverse complement of the primer as the template)
and count how frequently they appear in a fastq file and return how frequently the input primer is observed at all in the sequencing reads
and the top 10 most frequent sequence reads that start with the input primer and the extension sequence that matches it
'''

from collections import Counter
import sys
import re
import gzip
import argparse
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd 

LinkerStart='TGGAATTCTC'

def complement(s):
    '''this function will give the complement version of an inputed sequence (note NOT reverse complement, only complement)
    '''
    return s.upper().replace('A', 't').replace('T', 'a').replace('G', 'c').replace('C', 'g').upper()

def search_and_count_sequences(F, primer_sequence, template_input):
    '''
    this function will search sequences in the FASTQ file data that match the primer extended with the template sequence
    and count the occurrences of each matching sequence.

    F=List of lines from the FASTQ file.
    primer_sequence=The primer sequence to search for.
    template_sequence=The template sequence to extend the primer with.
    return=A Counter object with sequence counts, total number of reads with the primer at the 5' end at all
    '''
    #read only every 4th line because that contains sequence information
    Sequences1 = F[1::4]

    #template sequence input will be in form of template DNA input which will be reverse complement and include sequence of primer at 3' end
    #primer has 5bp overhang on 5' end and then region of complementarity with template from 3' end
    templated_sequence=complement(template_input)[::-1][len(primer_sequence)-5:]

    sequence_counts = Counter()
    total_with_primer=0
    extension_lengths = Counter()  

    #look for counts that only go to gapped/nicked point
    
    #loop through each sequence
    for seq in Sequences1:  
        if (primer_sequence in seq) and (LinkerStart in seq):
            pos=seq.rfind(LinkerStart)
            seq_trim=seq[:pos]

            if len(seq_trim)<20:
                continue
            # print('seq_trim: ', seq_trim)

            total_with_primer+=1
            sequence_counts[seq_trim]+=1

            start = seq_trim.find(primer_sequence)
            if start != 0:
                continue

            extension_portion=seq_trim[start+len(primer_sequence):]

            for i in range(len(templated_sequence),-1,-1):
                if extension_portion==templated_sequence[:i]:
                    extension_lengths[len(extension_portion)] += 1
                    break
                # extension_portion_pos=extension_portion.find(templated_sequence[:i])
                # if extension_portion_pos != -1:
                    

                # print('extension_length: ', extn_length)
                # print('templated_region: ', templated_region)

        #     if start != -1:
        #         end = start + len(primer_sequence) + r'[ATGC]*' + len(template_sequence)
        #         if end <= len(seq_trim) and seq_trim[start:end].endswith(template_sequence:
        #             sequence_counts[seq_trim] += 1

        #             extn_length = len(template_sequence)
        #             extension_lengths[extn_length] += 1
        #             print('extension_seq: ', )



        # #check if it's matching a sequence based on the primer/template extension
        # for i in range(len(primer_sequence+template_sequence)*2):
        #     partial_extended_sequence = primer_sequence + r'[ATGC]*' + template_sequence[:i] + r'[ATGC]*'

        #     options_for_extensions = re.compile(partial_extended_sequence)
            
        #     #check if extended sequence is found in the sequence
        #     if options_for_extensions.search(seq_trim):
        #         sequence_counts[seq_trim] += 1

        #         #find templated extension region after primer+potential stutter bases
        #         extension_portion = re.search(template_sequence, seq_trim)
        #         if extension_portion:

        #             #match contains templated extension region, need to know length
        #             extn_length = len(extension_portion.group())

        #             extension_lengths[extn_length] += 1  
                    
        #             break
        #         break

        # if len(seq)>(len(primer_sequence+template_sequence)+1) and primer_sequence in seq:
        #     print('**Longer sequence than templating alone: ', seq)
        #should find a way to determine if these regions are templated/rolling hairpin or circle 


    max_length = max(len(seq) for seq in sequence_counts)
    return sequence_counts, total_with_primer, extension_lengths

def plot_extension_seaborn(extension_lengths, primer_sequence, template_sequence, gap_vs_nick, sample, total_with_primer, show_legend=False, save_as_svg=False, save_as_png=False, filename="ds_extension_histogram.svg"):
    '''Plot the length of extension as a percentage of total_with_primer compared to where the gap/nick was in the template.
    
    extension_lengths: Counter object of extension lengths
    total_with_primer: Total number of reads that had the primer sequence (used to calculate percentage)
    gap_vs_nick: integer indicating if the input has a 3-base gap (3) or a 0-base nick (0)
    sample: sample name for labeling
    '''
    # Create a DataFrame for plotting
    df = pd.DataFrame(list(extension_lengths.items()), columns=['Extension Length', 'Count'])
    
    # Calculate percentage of total_with_primer for each count
    df['Percentage'] = (df['Count'] / total_with_primer) * 100

    min_length = min(df['Extension Length'].min(), 0)
    max_length = max(df['Extension Length'].max(), len(template_sequence))

    max_percentage = df['Percentage'].max()

    # Dynamic height adjustment (scaling factor for height based on the maximum percentage)
    base_height = 2  # Minimum height
    scaling_factor = 0.05  # Adjust this for more or less sensitivity
    dynamic_height = base_height + scaling_factor * max_percentage  # Final height

    # Set width dynamically based on number of extension lengths, using preset widths for 18, 21, 40, and 43
    if '83' in sample:
        figsize = (5, dynamic_height)  # Smallest width
    elif '82' in sample:
        figsize = (5.5, dynamic_height)  # Slightly wider
    elif '81' in sample:
        figsize = (9, dynamic_height)  # Medium width
    else:
        figsize = (9.5, dynamic_height)  # Largest width

    plt.figure(figsize=figsize)
    sns.histplot(
        data=df, 
        x='Extension Length', 
        weights='Percentage',  # Plot the percentages instead of counts
        bins=range(min_length, max_length + 1), 
        discrete=True, 
        color='#D48C8C', 
        edgecolor='#B00C0C'
    )
    
    # Add dashed lines for expected length and end of template
    expected_length = gap_vs_nick
    plt.axvline(expected_length, color='#FFD700', linestyle='--', label=f'Reach end of gap/start of double-stranded region ({expected_length} bases)')
    
    longest_expected_length = len(template_sequence) - 16
    plt.axvline(longest_expected_length, color='black', linestyle='--', label=f'Length to end of template (minus primer) ({longest_expected_length} bases)')

    # Set x-axis ticks and add padding on either side
    plt.xticks(np.arange(df['Extension Length'].min(), df['Extension Length'].max() + 1, 3))
    padding = 0.5  # Amount of padding to add
    plt.xlim(df['Extension Length'].min() - padding, df['Extension Length'].max() + padding)

    # Update the y-axis label for percentages
    plt.ylabel('Percentage of Total With Primer (%)')

    plt.title(f'Extension Lengths Distribution \n {sample}', fontsize=10)
    plt.xlabel('Length of Extension')

    # Conditionally show the legend
    if show_legend:
        # Position the legend below the plot
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=1, frameon=False)  # Adjust bbox_to_anchor and ncol for positioning

    # Save as SVG or PNG if requested
    if save_as_svg:
        plt.savefig(f"{filename}.svg", format='svg', bbox_inches='tight')
    elif save_as_png:
        plt.savefig(f"{filename}.png", format='png', bbox_inches='tight')
    else:
        plt.show()
def main():
    #take in .fastq.gz input file of sequences, primer sequence, and optional template sequence in the command line
    parser = argparse.ArgumentParser(description="Named parameters in script")
    
    parser.add_argument('--input_file', required=True, help="path to the input fastq file")
    parser.add_argument('--primer', required=True, help="sequence of primer")
    parser.add_argument('--template', required=False, help="sequence of template, optional")

    args = parser.parse_args()

    input_file = args.input_file
    #need to determine if nick or gapped template--gap or nick is an integer variable representing length of single-stranded region after primer before ds region of template
    #input tempalte 80 and 82 are gapped, 81 and 83 are nicked
    if '80' in input_file or '82' in input_file:
        gap_or_nick=3
    elif '81' in input_file or '83' in input_file:
        gap_or_nick=0  
    else:
        print('cannot find template number in input_file')      
    

    primer = str(args.primer).strip().upper()

    if 'template' in args:
        template=str(args.template).strip().upper()
    else:
        template=complement(primer).upper()

    with gzip.open(input_file, mode='rt') as f:
        F = f.read().splitlines()

    sequence_counts, total_with_primer, extension_lengths = search_and_count_sequences(F, primer, template)


    #make output file with total number of sequences that start with the primer and top 10 most common sequences and their counts and proportion
    sample = input_file.split('/')[-1].split('.')[0]

    plot_extension_seaborn(extension_lengths, primer, template, gap_or_nick, sample, total_with_primer, save_as_svg=True, show_legend=False, filename=f'dsextension_histogram_WITHLEGEND_{sample}_{datetime.now().strftime("%Y%m%d-%H%M%S")}.svg')

    output_filename = f'dsDragonMostCommon_{sample}_{datetime.now().strftime("%Y%m%d-%H%M%S")}.fa'

    with open(output_filename, mode='wt') as output:
        output.write('>total_counts_startwith_primer: ' + str(total_with_primer))
        output.write('\n')
        for i, (seq, count) in enumerate(sequence_counts.most_common(10)):
            # print(i, seq, count, 'i, seq, count\n')
            output.write(f'>{sample} {i} Counts:{count} Proportion of primer-start reads:{str(round(count/total_with_primer, 3))}\n')
            output.write(f'{seq}\n')

    print(f"Results written to {output_filename}")

if __name__ == "__main__":
    main()
