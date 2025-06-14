
# Generate primers for the selected target region
# Forward primer are selcted from upstream of the target
# Reverse primers are selected from downstream of the target
# Primer size should be from 18 to 25 bp

from insilico_pcr.target_extractor import handle_fasta, extract_target
from Bio.Seq import Seq

dna = handle_fasta(r"C:\Users\DELL\Documents\Personal\PORTFOLIO\Insilico-PCR-Tool\data\ACTB_cds.fasta")
target = extract_target(dna)

def generate_candidate_primers(dna, target_start, target_end, flank=50):

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
            reverse_primer_candidates.append(primer_reverse_complement)
    
       
    return forward_primer_candidates, reverse_primer_candidates




generate_candidate_primers(dna,281,385,50)