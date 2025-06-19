
# Generate primers for the selected target region
# Forward primer are selcted from upstream of the target
# Reverse primers are selected from downstream of the target
# Primer size should be from 18 to 25 bp

import itertools
import pandas as pd
from insilico_pcr.primer_validator import validate_primers, validate_primer_pairs
from Bio.Seq import Seq

def generate_candidate_primers(dna, target_start, target_end, flank):

    # Get Upstream flank

    upstream_start = max(0, target_start-flank-1)
    upstream_end = target_start - 1
    upstream_seq = dna[upstream_start : upstream_end]

    # Get downstream flank

    downstream_start = min(target_end+1, (len(dna)-1))
    downstream_end = target_end + flank + 1
    downstream_seq = dna[downstream_start : downstream_end]

    # Generate all possible forward and reverse primers

    forward_primer_candidates = []
    reverse_primer_candidates = []

    for length in range (18, 26):
        for i in range (0, len(upstream_seq)-length+1):
            primer = upstream_seq[i : i+length]
            forward_primer_candidates.append(primer)
    
        for j in range (0, len(downstream_seq)-length+1):
            primer = downstream_seq[j : j+length]
            primer_reverse_complement = Seq(primer).reverse_complement()
            reverse_primer_candidates.append(str(primer_reverse_complement))
    
       
    return forward_primer_candidates, reverse_primer_candidates


# Perform primer validation and keep only the valid ones

def filter_valid_primer_candidates(forward_primer_candidates, reverse_primer_candidates):
   
    valid_forward_primer_candidates = []
    valid_reverse_primer_candidates = []
    forward_primers_data = []
    reverse_primers_data = []

    for candidate in forward_primer_candidates:
        boolean_value, data_dict = validate_primers(candidate)
        forward_primers_data.append(data_dict)
        if boolean_value == True:
            valid_forward_primer_candidates.append(candidate)

    for candidate in reverse_primer_candidates:
        boolean_value, data_dict = validate_primers(candidate)
        reverse_primers_data.append(data_dict)
        if boolean_value == True:
            valid_reverse_primer_candidates.append(candidate)

    fp_df = pd.DataFrame(forward_primers_data)
    rp_df = pd.DataFrame(reverse_primers_data)

    fp_df.to_csv("Insilico-PCR-Tool/results/CandidateForwardPrimers_Data.csv", index=False)
    rp_df.to_csv("Insilico-PCR-Tool/results/CandidateReversePrimers_Data.csv", index=False)

    return valid_forward_primer_candidates, valid_reverse_primer_candidates



# Filter all valid primer pairs

def filter_valid_primer_pair_candidates(primer_pair_candidates):

    valid_primer_pair_candidates_list = []
    primer_pairs_data = []
    valid_data_list = []

    for primer_pair in primer_pair_candidates:
        fp = primer_pair[0]
        rp = primer_pair[1]
        boolean_value, data_dict = validate_primer_pairs(fp, rp)
        primer_pairs_data.append(data_dict)
        if boolean_value == True:
            valid_pair = (fp, rp)
            valid_primer_pair_candidates_list.append(valid_pair)
            valid_data_list.append(data_dict)

    # Data for all primer pairs
    primer_pair_data_df = pd.json_normalize(primer_pairs_data)
    primer_pair_data_df.to_csv("Insilico-PCR-Tool/results/CandidatePrimerPairs_Data.csv", index=False)

    # Data for only valid primer pairs
    valid_primer_pair_data_table = pd.json_normalize(valid_data_list)


    return valid_primer_pair_candidates_list, valid_primer_pair_data_table





