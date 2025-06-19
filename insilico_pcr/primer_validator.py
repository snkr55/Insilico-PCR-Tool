# Validate all the primers based on certain criteria
# Criteria to be checked:
# (1) Primer Length - 18 to 25 bp
# (2) GC content - 40 to 60 %
# (3) Melting Temperature - 50 to 65 °C 
# (5) Homopolymers - Avoid >3 repeated bases (e.g. AAAA)
# (8) Hairpins - Avoid significant self-binding
# (6) Self-dimers - Avoid strong self-dimers
# (4) Tm difference - ΔTm ≤ 2°C between FWD and REV
# (7) Cross-dimers - Avoid between FWD and REV


import primer3
from Bio.SeqUtils import gc_fraction
import re

# Parameters
MIN_LEN = 18
MAX_LEN = 25
MIN_GC = 40
MAX_GC = 60
MIN_TM = 50
MAX_TM = 65
MAX_HOMOPOLYMER = 3
MAX_TM_DIFF = 2.0
MAX_DIMER_DG = -6.0
MAX_HAIRPIN_DG = -2.0

def has_homopolymer(primer, limit=MAX_HOMOPOLYMER):
    bases = ['A','T','G','C']
    for base in bases:
        if base*(limit+1) in primer:
            return True # Found homopolymer
    return False # No homopolymer


def calculate_primer_characteristics(primer):
    length = len(primer)
    gc = round(gc_fraction(primer) * 100, 3)
    tm = round(primer3.calc_tm(primer), 3)
    homopolymer = has_homopolymer(primer)
    hairpin = primer3.calc_hairpin(primer)
    homodimer = primer3.calc_homodimer(primer)

    results = {
        'primer': primer,
        'length': length,
        'gc': gc,
        'tm': tm,
        'has_homopolymer': homopolymer,
        'hairpin_dG': round(hairpin.dg, 3) if hairpin.structure_found else 0,
        'hairpin_found': hairpin.structure_found,
        'homodimer_dG': round(homodimer.dg, 3) if homodimer.structure_found else 0,
        'homodimer_found': homodimer.structure_found,
    }

    return results


def validate_primers(primer):
    primer_data = calculate_primer_characteristics(primer)
    if not(MIN_LEN <= primer_data['length'] <= MAX_LEN) or\
        not(MIN_GC <= primer_data['gc'] <= MAX_GC) or \
        not(MIN_TM <= primer_data['tm'] <= MAX_TM) or \
        primer_data['has_homopolymer'] == True or \
        primer_data['hairpin_found'] == True:
        #primer_data['homodimer_found'] == True:
        return False, primer_data
    
    return True, primer_data


def calculate_primer_pair_characteristics(fp ,rp):
    fwdp_data = calculate_primer_characteristics(fp)
    revp_data = calculate_primer_characteristics(rp)
    
    delta_tm = abs(fwdp_data['tm'] - revp_data['tm'])
    cross_dimer = primer3.calc_heterodimer(fp, rp)

    results = {
        'fp': fwdp_data,
        'rp': revp_data,
        'cross_dimer_found': cross_dimer.structure_found,
        'cross_dimer_dG': cross_dimer.dg if cross_dimer.structure_found else 0,
        'delta_tm': delta_tm
    }

    return results

def validate_primer_pairs(fp, rp):
    primer_pair_data = calculate_primer_pair_characteristics(fp, rp)
    if primer_pair_data['delta_tm'] > MAX_TM_DIFF: \
        #or primer_pair_data['cross_dimer_found'] == True:
        return False, primer_pair_data
    
    return True, primer_pair_data

