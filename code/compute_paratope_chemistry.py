#!/usr/bin/env python3
"""Analyze paratope chemical composition and create hydrophobic patch map."""

import csv
from collections import Counter

# Define paratope residues from contact analysis
paratope_residues = {
    'A': [31, 32, 33, 50, 52, 55, 57, 59, 99, 102, 103, 104, 105, 106, 107, 109,
          2, 27, 28, 30, 58, 98, 100, 101, 108, 110],
    'B': [2, 27, 31, 98, 99, 100, 102, 1, 26, 28, 29, 30, 32, 33, 34, 35, 36, 38,
          55, 56, 96, 97, 101]
}

# Conserved interface residues (bind both antigens)
conserved_residues = {
    'A': [31, 32, 33, 50, 52, 55, 57, 59, 99, 102, 103, 104, 105, 106, 107, 109],
    'B': [2, 27, 31, 98, 99, 100, 102]
}

# Residue classifications
hydrophobic = {'PHE', 'TYR', 'TRP', 'LEU', 'ILE', 'VAL', 'ALA', 'MET'}
aromatic = {'PHE', 'TYR', 'TRP', 'HIS'}
basic = {'ARG', 'LYS', 'HIS'}
acidic = {'ASP', 'GLU'}
polar = {'SER', 'THR', 'ASN', 'GLN', 'CYS'}

def parse_pdb_residues(pdb_file):
    """Extract residue types from PDB."""
    residues = {}
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                chain = line[21].strip()
                resnum = int(line[22:26].strip())
                resname = line[17:20].strip()
                
                if chain not in residues:
                    residues[chain] = {}
                residues[chain][resnum] = resname
    return residues

def analyze_composition(residues, res_dict, label):
    """Analyze chemical composition of residue set."""
    res_types = []
    for chain, positions in res_dict.items():
        for pos in positions:
            if chain in residues and pos in residues[chain]:
                res_types.append(residues[chain][pos])
    
    counts = Counter(res_types)
    total = len(res_types)
    
    # Calculate fractions
    hydro_count = sum(1 for r in res_types if r in hydrophobic)
    arom_count = sum(1 for r in res_types if r in aromatic)
    basic_count = sum(1 for r in res_types if r in basic)
    acidic_count = sum(1 for r in res_types if r in acidic)
    polar_count = sum(1 for r in res_types if r in polar)
    
    print(f"\n=== {label} Composition (n={total}) ===")
    print(f"Hydrophobic: {hydro_count} ({100*hydro_count/total:.1f}%)")
    print(f"Aromatic:    {arom_count} ({100*arom_count/total:.1f}%)")
    print(f"Basic:       {basic_count} ({100*basic_count/total:.1f}%)")
    print(f"Acidic:      {acidic_count} ({100*acidic_count/total:.1f}%)")
    print(f"Polar:       {polar_count} ({100*polar_count/total:.1f}%)")
    
    print(f"\nTop residue types:")
    for res, count in counts.most_common(10):
        print(f"  {res}: {count} ({100*count/total:.1f}%)")
    
    return {
        'total': total,
        'hydrophobic': hydro_count,
        'aromatic': arom_count,
        'basic': basic_count,
        'acidic': acidic_count,
        'polar': polar_count,
        'residue_counts': dict(counts)
    }

def create_hydrophobic_map(residues, res_dict):
    """Create CSV mapping hydrophobic patches."""
    rows = []
    for chain in ['A', 'B']:
        if chain in res_dict:
            for pos in sorted(res_dict[chain]):
                if chain in residues and pos in residues[chain]:
                    resname = residues[chain][pos]
                    is_hydro = resname in hydrophobic
                    is_arom = resname in aromatic
                    is_basic = resname in basic
                    is_acidic = resname in acidic
                    is_polar = resname in polar
                    
                    rows.append({
                        'Chain': chain,
                        'Position': pos,
                        'ResName': resname,
                        'Hydrophobic': is_hydro,
                        'Aromatic': is_arom,
                        'Basic': is_basic,
                        'Acidic': is_acidic,
                        'Polar': is_polar,
                        'In_Conserved': chain in conserved_residues and pos in conserved_residues[chain]
                    })
    return rows

def main():
    # Parse both complex structures
    print("Parsing CL1-Ab4_p53.pdb...")
    residues_p53 = parse_pdb_residues('/Users/aloksinha/Desktop/CL1-Ab4_p53.pdb')
    
    print("Parsing CL1-Ab4_CALR.pdb...")
    residues_calr = parse_pdb_residues('/Users/aloksinha/Desktop/CL1-Ab4_CALR.pdb')
    
    # Analyze paratope composition
    print("\n" + "="*60)
    print("PARATOPE CHEMICAL COMPOSITION ANALYSIS")
    print("="*60)
    
    para_stats_p53 = analyze_composition(residues_p53, paratope_residues, "Full Paratope (p53 complex)")
    cons_stats_p53 = analyze_composition(residues_p53, conserved_residues, "Conserved Interface (p53 complex)")
    
    para_stats_calr = analyze_composition(residues_calr, paratope_residues, "Full Paratope (CALR complex)")
    cons_stats_calr = analyze_composition(residues_calr, conserved_residues, "Conserved Interface (CALR complex)")
    
    # Create hydrophobic map CSV
    hydro_map = create_hydrophobic_map(residues_p53, paratope_residues)
    
    with open('/Users/aloksinha/Desktop/paratope_chemistry_map.csv', 'w', newline='') as csvfile:
        fieldnames = ['Chain', 'Position', 'ResName', 'Hydrophobic', 'Aromatic', 'Basic', 'Acidic', 'Polar', 'In_Conserved']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in hydro_map:
            writer.writerow(row)
    
    print(f"\n✓ Generated paratope_chemistry_map.csv with {len(hydro_map)} residues")
    
    # Key findings
    print("\n" + "="*60)
    print("KEY FINDINGS")
    print("="*60)
    
    print("\n1. PARATOPE ENRICHMENT:")
    hydro_frac = cons_stats_p53['hydrophobic'] / cons_stats_p53['total']
    arom_frac = cons_stats_p53['aromatic'] / cons_stats_p53['total']
    basic_frac = cons_stats_p53['basic'] / cons_stats_p53['total']
    
    print(f"   - Hydrophobic enrichment: {100*hydro_frac:.1f}%")
    print(f"   - Aromatic enrichment: {100*arom_frac:.1f}%")
    print(f"   - Basic enrichment: {100*basic_frac:.1f}%")
    
    if arom_frac > 0.25:
        print("\n   → HIGH AROMATIC CONTENT suggests π-π and cation-π interactions")
    if basic_frac > 0.20:
        print("   → HIGH BASIC CONTENT suggests electrostatic/salt bridge promiscuity")
    if hydro_frac > 0.50:
        print("   → HIGH HYDROPHOBIC CONTENT suggests hydrophobic patch-mediated binding")
    
    print("\n2. POLYREACTIVITY INDICATORS:")
    if cons_stats_p53['aromatic'] >= 3 and cons_stats_p53['basic'] >= 2:
        print("   ✓ Multiple aromatics + basics → typical polyreactive signature")
    
    print("\n3. CONSERVED RESIDUE CHEMISTRY:")
    print("   Residues that bind both antigens are enriched in:")
    for res, count in sorted(cons_stats_p53['residue_counts'].items(), key=lambda x: x[1], reverse=True)[:5]:
        print(f"   - {res}: {count} positions")

if __name__ == '__main__':
    main()

