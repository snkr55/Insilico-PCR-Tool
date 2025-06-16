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

    return {
        'length': length,
        'gc': gc,
        'tm': tm,
        'has_homopolymer': homopolymer,
        'hairpin_dG': round(hairpin.dg, 3) if hairpin.structure_found else 0,
        'hairpin_found': hairpin.structure_found,
        'homodimer_dG': round(homodimer.dg, 3) if homodimer.structure_found else 0,
        'homodimer_found': homodimer.structure_found,
    }


def validate_primers(primer):
    primer_data = calculate_primer_characteristics(primer)
    if primer_data['length'] not in range (MIN_LEN, MAX_LEN+1):
        return False
    if primer_data['gc'] not in range (MIN_GC, MAX_GC+1):
        return False
    if primer_data['tm'] not in range (MIN_TM, MAX_TM+1):
        return False
    if primer_data['has_homopolymer'] == True:
        return False
    if primer_data['hairpin_found'] == True:
        return False
    if primer_data['homodimer_found'] == True:
        return False
    
    return True


def calculate_primer_pair_characteristics(fwdp,revp):
    fwdp_data = calculate_primer_characteristics(fwdp)
    revp_data = calculate_primer_characteristics(revp)
    
    delta_tm = abs(fwdp_data['tm'] - revp_data['tm'])
    cross_dimer = primer3.calcHeterodimer(fwdp, revp)
    cross_dimer_dG = cross_dimer.dg if cross_dimer.structure_found else 0

    results = {
        'forward': fwdp_data,
        'reverse': revp_data,
        'cross_dimer_dG': cross_dimer_dG,
        'delta_tm': delta_tm,
        'passes': True
    }

    return results



#primer_data = calculate_primer_characteristics("GAGAAAATCTGGCACCACACCTTCTACAATG")
#print(primer_data)
#validate_primers("GAGAAAATCTGGCACCACACCTTCTACAATG")

