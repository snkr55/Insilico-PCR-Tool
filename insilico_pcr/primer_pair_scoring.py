def score_primer_pairs(valid_primer_pairs_list, valid_primer_pair_data):

    for i,pair in enumerate(valid_primer_pairs_list):

        fp = pair[0]
        rp = pair[1]
        fp_length = valid_primer_pair_data.loc[i,'fp.length']
        rp_length = valid_primer_pair_data.loc[i,'rp.length']
        delta_tm = valid_primer_pair_data.loc[i,'delta_tm']
        amplicon_length = valid_primer_pair_data.loc[i,'amplicon_length']

        score = 0

        # 1. Primer length
        if 18 <= fp_length <= 21:
            score += 2
        if 22 <= fp_length <= 25:
            score += 1
        if 18 <= rp_length <= 21:
            score += 2
        if 22 <= rp_length <= 25:
            score += 1
       
        # 2. Melting temperature difference
        if delta_tm < 1:
            score += 2
        if 1 <= delta_tm <= 2:
            score += 1

        # 3. GC Clamp
        if fp[-1] in ['G', 'C']:
            score += 1
        if rp[-1] in ['G', 'C']:
            score += 1
        
        # 4. Amplicon Length
        if amplicon_length < 100:
            score -= 4
        if 100 <= amplicon_length <= 200:
            score += 4
        if 201 <= amplicon_length <= 220:
            score += 3
        if 221 <= amplicon_length <= 240:
            score += 2
        if 241 <= amplicon_length <= 260:
            score += 1
        if amplicon_length > 260:
            score -= 4

        # Add score to the table
        valid_primer_pair_data.loc[i, 'score'] = score

    # Sort the updated table by score
    scored_primer_pair_data = valid_primer_pair_data.sort_values(by=['score','delta_tm','amplicon_length'], ascending=[False,True,True])
    scored_primer_pair_data = scored_primer_pair_data.reset_index(drop=True)
    scored_primer_pair_data.index = scored_primer_pair_data.index + 1
    scored_primer_pair_data.index.name = 'Rank'

    scored_primer_pair_data.to_csv("results/ScoredPrimerPairs_Data.csv")

    return scored_primer_pair_data

        
