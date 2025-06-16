
# Generate primers for the selected target region
# Forward primer are selcted from upstream of the target
# Reverse primers are selected from downstream of the target
# Primer size should be from 18 to 25 bp

from insilico_pcr.primer_validator import validate_primers
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

def generate_valid_primer_candidates(forward_primer_candidates, reverse_primer_candidates):
   
    valid_forward_primer_candidates = []
    valid_reverse_primer_candidates = []

    for candidate in forward_primer_candidates:
        if validate_primers(candidate) == True:
            valid_forward_primer_candidates.append(candidate)

    for candidate in reverse_primer_candidates:
        if validate_primers(candidate) == True:
            valid_reverse_primer_candidates.append(candidate)

    #print(f"Valid Candidate Primers count: {len(valid_forward_primer_candidates)}, {len(valid_forward_primer_candidates)}")
    
    return valid_forward_primer_candidates, valid_reverse_primer_candidates