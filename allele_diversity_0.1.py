#!/usr/bin/env python

from Bio import SeqIO
import math

def hamming_distance(s1, s2):
    """Return the Hamming distance between equal-length sequences."""
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length.")
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))
    
working_folder = "/Users/peter/Documents/Sequence_data/MAUI_Sara_2016/results/MAUIcount_output_nodA/"
sequence_file = working_folder + "accepted_sequences.fas"
table_file = working_folder + "nodA1_sequences.tab"
output_file = working_folder + "nodA1_allele_diversity.tab"
distance_matrix = working_folder + "nodA1_distance_matrix.tab"

distance_dict = {}
seq_list = list(SeqIO.parse(sequence_file, "fasta"))

seq_len = len(seq_list[0])
#NB: this assumes that the alignment is ungapped (true for our data)

with open(distance_matrix, "w") as outfile2:

    for j in range(len(seq_list)):            
        sj = seq_list[j].id
        ssj = sj[:sj.find("_",4)]
        outfile2.write("\t"+ssj)
    outfile2.write("\n")
        
    for i in range(len(seq_list)):

        si = seq_list[i].id
        ssi = si[:si.find("_",4)]
        outfile2.write(ssi)

        distance_dict[ssi] = {}
        for j in range(len(seq_list)):
            hd = hamming_distance(seq_list[i],seq_list[j])
            sj = seq_list[j].id
            ssj = sj[:sj.find("_",4)]
            outfile2.write('\t{:d}'.format(hd))
            distance_dict[ssi][ssj] = hd
        outfile2.write("\n")
        
output_table = []
with open(table_file) as infile:
    for line in infile:
        sample_data = line.rstrip("\n\t").split("\t")
        sample_name = sample_data.pop(0)
        if sample_name == "":
            seq_name = sample_data
        elif sample_name == "total": 
            break
        else:
            sample_data = list(map(int, sample_data))
            sample_freq = []
            max_dist = 0
            total_counts = sum(sample_data)
            if total_counts == 0:
                mpd = max_dist = pi = math.nan
                allele_count = 0
            else:
                num_seqs = len(sample_data)
                mpd = 0
                for count in sample_data:
                    sample_freq.append(count/total_counts)
                allele_count = num_seqs - sample_data.count(0)
                for i in range(num_seqs):
                    for j in range(num_seqs):
                        mpd += sample_freq[i]*sample_freq[j]*distance_dict[seq_name[i]][seq_name[j]]
                        if sample_freq[i]*sample_freq[j] > 0:
                            max_dist = max(max_dist,distance_dict[seq_name[i]][seq_name[j]])
                pi = mpd/seq_len            
                
            output_table.append((sample_name, total_counts, allele_count, mpd, max_dist, pi))
    
with open(output_file, "w") as outfile:
    outfile.write("sample\tcounts\talleles\tmpd\tmax_dist\tpi\n") 
    for (sample_name, total_counts, allele_count, mpd, max_dist, pi) in output_table:
        outfile.write('{:10}\t{:6d}\t{:d}\t{:>5.4F}\t{:F}\t{:>7.6F}\n'.format(sample_name, total_counts, allele_count, mpd, max_dist, pi))  
        
        
