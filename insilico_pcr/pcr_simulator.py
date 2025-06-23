
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt 


# Parameters

tm_range = (50,65)
min_amplicon_size = 100
max_amplicon_size = 250
MAX_MISMATCH_ALLOWED = 1


# Primer matching function

def primer_match_positions(primer,dna,max_mismatch_allowed=MAX_MISMATCH_ALLOWED):
    match_positions = []
    for i in range(len(dna) - len(primer) + 1):
        window = dna[i : i + len(primer)]
        #print(i, window)
        #print(list(zip(primer, window)))
        mismatch = 0
        for a,b in zip(primer,window):
            if a!=b:
                #print(a,b)
                mismatch+=1
        mismatch_count = mismatch
        #print(i,mismatch)
        if mismatch_count <= max_mismatch_allowed:
            match_positions.append(i)
    # print(matches)

    return match_positions


# Running In-silico PCR simulation function

def run_insilico_pcr(dna,fp,rp):

    rp_rev_comp = str(Seq(rp).reverse_complement())

    forward_match_positions = primer_match_positions(fp,dna,MAX_MISMATCH_ALLOWED)
    # print(f"\nForward Match Positions: {forward_match_positions}")
    reverse_match_positions = primer_match_positions(rp_rev_comp,dna,MAX_MISMATCH_ALLOWED)
    # print(f"Reverse Match Positions: {reverse_match_positions}")

    results = []
    for forward_start in forward_match_positions:
        for reverse_start in reverse_match_positions:
            if reverse_start > forward_start:
                amplicon_length = (reverse_start + len(rp)) - forward_start + 1
                # print(amplicon_length)
                if min_amplicon_size <= amplicon_length <= max_amplicon_size:
                    amplicon = dna[forward_start : (reverse_start + len(rp))]
                    # print(amplicon)
                    results.append({'start': forward_start,
                                    'end': reverse_start + len(rp),
                                    'length': amplicon_length,
                                    'amplicon_seq': str(amplicon)})
    return results