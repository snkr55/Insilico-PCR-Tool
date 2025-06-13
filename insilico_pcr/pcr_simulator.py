# To perform in-silico PCR

print("In-silico PCR")

# Importing Libraries

from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt 


# User input


dna = input("Enter template sequence: ").strip().upper().replace(" ", "")
fp = input("Enter forward primer (5' to 3'): ").strip().upper()
rp = input("Enter reverse primer (5' to 3'): ").strip().upper()

#dna = Seq("AATGCTACCGGGCCGCTATCACGATGCGGGCGCTCGCTAGCGTTGATCTACCGGGCCGCTATCAGGTACGATCGTAGCTGCGGGCGCTCGCTAGCAGCTTACCGGGCCGCTATCAGGCGGGCGCTCGCTAGCATCG")
#fp = "TACCGGGCCGCTATCA"
#rp = "GCTAGCGAGCGCCCGC"


# Parameters


tm_range = (50,65)
min_amplicon_size = 20
max_amplicon_size = 50
max_mismatch_allowed = 1


# Primer matching function


def primer_match_positions(primer,dna,max_mismatch_allowed):
    matches = []
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
            matches.append(i)
    #print(matches)

    return matches


# Running In-silico PCR simulation function


def run_insilico_pcr(dna,fp,rp):

    fp_tm = round(mt.Tm_NN(Seq(fp)), 3)
    rp_tm = round(mt.Tm_NN(Seq(rp)), 3)
    print(f"Forward Primer Tm: {fp_tm}°C")
    print(f"Reverse Primer Tm: {rp_tm}°C")
    if not (tm_range[0] <= fp_tm <= tm_range[1]) or not (tm_range[0] <= rp_tm <= tm_range[1]):
        print(f"The melting temperature of one or both primers are not within the compatible Tm range ({tm_range[0]}°C-{tm_range[1]}°C)")
        print(f"Forward Primer Tm: {fp_tm}°C")
        print(f"Reverse Primer Tm: {rp_tm}°C")
        return None
    else:
        print(f"The melting temperature of both primers are within the compatible Tm range ({tm_range[0]}°C-{tm_range[1]}°C)")

    rp_rev_comp = Seq(rp).reverse_complement()

    forward_match_positions = primer_match_positions(fp,dna,max_mismatch_allowed)
    #print(forward_match_positions)
    reverse_match_positions = primer_match_positions(rp_rev_comp,dna,max_mismatch_allowed)
    #print(reverse_match_positions)

    results = []
    for forward_start in forward_match_positions:
        for reverse_start in reverse_match_positions:
            if reverse_start > forward_start:
                amplicon_length = (reverse_start + len(rp)) - forward_start
                #print(amplicon_length)
                if min_amplicon_size <= amplicon_length <= max_amplicon_size:
                    amplicon = dna[forward_start : (reverse_start + len(rp))]
                    #print(amplicon)
                    results.append({'start': forward_start,
                                    'end': reverse_start + len(rp),
                                    'length': amplicon_length,
                                    'amplicon_seq': str(amplicon)})
    return results



pcr_products = run_insilico_pcr(dna,fp,rp)
print("PCR Product(s):")
print(f"-"*80)
for index, product in enumerate(pcr_products, 1):
    print(f"Result {index}")
    print(f"Start Position: {product['start']} | End Position: {product['end']} | Amplicon length: {product['length']}")
    print(f"PCR Product: {product['amplicon_seq']}")
    print(f"-"*80)



