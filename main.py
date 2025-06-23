
import pandas as pd
from insilico_pcr.target_extractor import *
from insilico_pcr.primer_generator import *
from insilico_pcr.primer_validator import *
from insilico_pcr.primer_pair_scoring import score_primer_pairs
from insilico_pcr.utilitis import *
from insilico_pcr.pcr_simulator import run_insilico_pcr
# from insilico_pcr.plot import *


dna = handle_fasta(r"data\ACTB_cds.fasta")
TARGET_SEQ, TARGET_START, TARGET_END, FLANK = extract_target(dna)


# Generation of valid forward and reverse primer candiates
forward_primer_candidates, reverse_primer_candidates = generate_candidate_primers(dna=dna, target_start=TARGET_START, target_end=TARGET_END, flank=FLANK)
print(f"\nNo. of forward and reverse primer candidates: {len(forward_primer_candidates)}, {len(reverse_primer_candidates)}")

valid_forward_primer_candidates, valid_reverse_primer_candidates = filter_valid_primer_candidates(forward_primer_candidates, reverse_primer_candidates)
print(f"No. of valid forward and reverse primers: {len(valid_forward_primer_candidates)}, {len(valid_reverse_primer_candidates)}")


# Generation of valid primer pair candidates
primer_pair_candidates = list(itertools.product(valid_forward_primer_candidates, valid_reverse_primer_candidates))
print(f"\nNo. of primer pair candidates: {len(primer_pair_candidates)}")

valid_primer_pairs, valid_pair_data = filter_valid_primer_pair_candidates(primer_pair_candidates)
print(f"No. of valid primer pairs: {len(valid_primer_pairs)}")


# Add amplicon info to the valid primer pair data table
primerPair_Amplicon_data = add_amplicon_columns(valid_pair_data, dna)


# Scoring the valid primer pairs
scored_primer_pairs = score_primer_pairs(valid_primer_pairs, primerPair_Amplicon_data)


# Print top primer pairs
num_to_print = int(input("\nEnter the number of top primer pairs you wish to view: "))
top_5_pairs = scored_primer_pairs[['fp.primer','rp.primer','fp.length','rp.length','delta_tm','amplicon_length','score']].head(num_to_print)
print(f'The top {num_to_print} primer pairs are:\n{top_5_pairs}')


# Select the primer pair to perform PCR simulation
rank, fp, rp = select_primer_pair(scored_primer_pairs)

# Simulate PCR
pcr_products = run_insilico_pcr(dna,fp,rp)
print("\nPCR Product(s):")
for index, product in enumerate(pcr_products, 1):
    print(f"\nResult {index}")
    print(f"Start Position: {product['start']} | End Position: {product['end']} | Amplicon length: {product['length']}")
    print(f"PCR Product: {product['amplicon_seq']}")
    print(f"-"*100)