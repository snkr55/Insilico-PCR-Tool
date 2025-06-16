
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


from insilico_pcr.target_extractor import *
from insilico_pcr.primer_generator import *
from insilico_pcr.primer_validator import *
#from insilico_pcr.pcr_simulator import *

dna = handle_fasta(r"C:\Users\DELL\Documents\Personal\PORTFOLIO\Insilico-PCR-Tool\data\ACTB_cds.fasta")
target = extract_target(dna)
TARGET_START = target[1]
TARGET_END = target[2]
FLANK = target[3]

primer_candidates = generate_candidate_primers(dna=dna, target_start=TARGET_START, target_end=TARGET_END, flank=FLANK)
forward_primer_candidates = primer_candidates[0]
reverse_primer_candidates = primer_candidates[1]
print(f"\nCandidate Primers count: {len(forward_primer_candidates)}, {len(reverse_primer_candidates)}")

valid_primer_candidates = generate_valid_primer_candidates(forward_primer_candidates, reverse_primer_candidates)
valid_forward_primer_candidates = valid_primer_candidates[0]
valid_reverse_primer_candidates = valid_primer_candidates[1]
print(f"\nValid Candidate Primers count: {len(valid_forward_primer_candidates)}, {len(valid_forward_primer_candidates)}")
