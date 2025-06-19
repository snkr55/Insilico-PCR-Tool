
from Bio.Seq import Seq




# Updating primer pair data table with amplicon information

def add_amplicon_columns(primer_pair_data_table, dna):

    df = primer_pair_data_table.copy()
    amplicon_starts = []
    amplicon_ends = []
    amplicon_lengths = []

    for index, row in df.iterrows():
        fp_seq = row['fp.primer']
        rp = row['rp.primer']
        rp_seq = str(Seq(rp).reverse_complement())

        fwd_index = dna.find(fp_seq)
        rev_index = dna.find(rp_seq)

        if fwd_index != -1 and rev_index != -1:
            amplicon_start = fwd_index
            amplicon_end = rev_index + len(rp_seq)
            amplicon_length = amplicon_end - amplicon_start + 1
        else:
            amplicon_start = None
            amplicon_end = None
            amplicon_length = None

        amplicon_starts.append(amplicon_start)
        amplicon_ends.append(amplicon_end)
        amplicon_lengths.append(amplicon_length)

    df['amplicon_start'] = amplicon_starts
    df['amplicon_end'] = amplicon_ends
    df['amplicon_length'] = amplicon_lengths

    # df.to_csv("Insilico-PCR-Tool/results/ValidPrimerPairs_Data.csv", index=False)

    return df




# Getting user input for selecting primer pair for PCR simulation

def select_primer_pair(scoredPrimerPairs_df):

    while True:
        
        selected_rank = int(input("\nEnter the Rank number of the primer pair you wish to use for PCR simulation: "))

        if selected_rank not in scoredPrimerPairs_df.index:
            print(f"\nInvalid Rank number. Please enter a valid Rank from 1 to {len(scoredPrimerPairs_df)}")
            continue
        
        else:
            print("\nSelected Primer Pair:")
            fp = scoredPrimerPairs_df.loc[selected_rank,'fp.primer']
            rp = scoredPrimerPairs_df.loc[selected_rank,'rp.primer']
            print(f"Forward Primer: {fp}")
            print(f"Primer Length: {scoredPrimerPairs_df.loc[selected_rank,'fp.length']} | Tm: {scoredPrimerPairs_df.loc[selected_rank,'fp.tm']} | GC content: {scoredPrimerPairs_df.loc[selected_rank,'fp.gc']}")
            print(f"Reverse Primer: {rp}")
            print(f"Primer Length: {scoredPrimerPairs_df.loc[selected_rank,'rp.length']} | Tm: {scoredPrimerPairs_df.loc[selected_rank,'rp.tm']} | GC content: {scoredPrimerPairs_df.loc[selected_rank,'rp.gc']}")
            print(f"Amplicon Length: {scoredPrimerPairs_df.loc[selected_rank,'amplicon_length']}")
            print(f"SCORE: {scoredPrimerPairs_df.loc[selected_rank,'score']}")
            return selected_rank, fp, rp


    
    