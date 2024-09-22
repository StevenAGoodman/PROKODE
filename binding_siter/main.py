import numpy as np
import pandas as pd
import re

def get_deltaG_score(motif, site):
    return None
def binding_siter(motifs_loc, fasta_seq_loc, match_thresh):
    # process fasta_seq into dict
    seq_dict = {}
    with open(fasta_seq_loc, 'r') as fasta_seq_file:
        for _, line in enumerate(fasta_seq_file):
            is_header = line.find('>')

            if is_header != -1: # if line is header line
                current_key = line[1:].replace('\n','')

            elif re.search('[A,T,G,C]+\n', line) != None: # if line is seq
                seq_dict[current_key] = line.replace('\n','')

            elif line == '\n': # if line is blank
                continue

            else:
                raise Exception('line in fasta seq not recognized')
            
    # preprocess motifs
    
    # scan over each dict seq and  write results to csv
    with open('./results.csv', 'a') as results_file:
        results_file.write('binder,bindee,kd')
        for key, seq in seq_dict.items:

            for i in range(len(seq) - len(seq) + 1):
                matches = []
                substring = seq[i:i+len(motif)]
                
                binding_energy = get_deltaG_score(motif, substring)

                if binding_energy >= match_thresh:
                    matches.append(binding_energy)
                    seq_bound = True
                
                if seq_bound:
                    results_file.write(f'{key}')
                
    #
    




binding_siter("TAAAACTACG", '../src/preprocessing/promoters.fa', 0.5)