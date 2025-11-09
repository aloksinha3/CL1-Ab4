#!/usr/bin/env python3
"""Analyze B-factors to assess paratope flexibility."""

import csv
import numpy as np

# Define paratope residues
paratope_residues = {
    'A': [31, 32, 33, 50, 52, 55, 57, 59, 99, 102, 103, 104, 105, 106, 107, 109,
          2, 27, 28, 30, 58, 98, 100, 101, 108, 110],
    'B': [2, 27, 31, 98, 99, 100, 102, 1, 26, 28, 29, 30, 32, 33, 34, 35, 36, 38,
          55, 56, 96, 97, 101]
}

# Conserved interface residues
conserved_residues = {
    'A': [31, 32, 33, 50, 52, 55, 57, 59, 99, 102, 103, 104, 105, 106, 107, 109],
    'B': [2, 27, 31, 98, 99, 100, 102]
}

def parse_bfactors(pdb_file):
    """Extract B-factors from PDB file."""
    bfactors = {}
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                chain = line[21].strip()
                resnum = int(line[22:26].strip())
                resname = line[17:20].strip()
                atom = line[12:16].strip()
                bfactor = float(line[60:66].strip())
                
                if chain not in bfactors:
                    bfactors[chain] = {}
                if resnum not in bfactors[chain]:
                    bfactors[chain][resnum] = {'resname': resname, 'atoms': []}
                
                bfactors[chain][resnum]['atoms'].append((atom, bfactor))
    
    # Calculate average B-factor per residue
    avg_bfactors = {}
    for chain in bfactors:
        avg_bfactors[chain] = {}
        for resnum in bfactors[chain]:
            atom_bfactors = [b for a, b in bfactors[chain][resnum]['atoms']]
            avg_bfactors[chain][resnum] = {
                'resname': bfactors[chain][resnum]['resname'],
                'avg_bfactor': np.mean(atom_bfactors),
                'std_bfactor': np.std(atom_bfactors),
                'min_bfactor': np.min(atom_bfactors),
                'max_bfactor': np.max(atom_bfactors),
                'n_atoms': len(atom_bfactors)
            }
    
    return avg_bfactors

def analyze_bfactor_distribution(bfactors, res_dict, label):
    """Analyze B-factor distribution for a residue set."""
    values = []
    for chain, positions in res_dict.items():
        for pos in positions:
            if chain in bfactors and pos in bfactors[chain]:
                values.append(bfactors[chain][pos]['avg_bfactor'])
    
    if not values:
        print(f"No B-factors found for {label}")
        return None
    
    print(f"\n=== {label} B-factor Statistics ===")
    print(f"Mean:   {np.mean(values):.2f} Ų")
    print(f"Median: {np.median(values):.2f} Ų")
    print(f"StdDev: {np.std(values):.2f} Ų")
    print(f"Min:    {np.min(values):.2f} Ų")
    print(f"Max:    {np.max(values):.2f} Ų")
    
    # Flexibility classification
    high_flex = sum(1 for v in values if v > 50)
    med_flex = sum(1 for v in values if 30 <= v <= 50)
    low_flex = sum(1 for v in values if v < 30)
    
    print(f"\nFlexibility distribution:")
    print(f"  High (>50):  {high_flex} ({100*high_flex/len(values):.1f}%)")
    print(f"  Medium (30-50): {med_flex} ({100*med_flex/len(values):.1f}%)")
    print(f"  Low (<30):   {low_flex} ({100*low_flex/len(values):.1f}%)")
    
    return {
        'mean': np.mean(values),
        'median': np.median(values),
        'std': np.std(values),
        'min': np.min(values),
        'max': np.max(values),
        'high_flex_count': high_flex,
        'med_flex_count': med_flex,
        'low_flex_count': low_flex
    }

def create_bfactor_csv(bfactors, res_dict, output_file):
    """Create CSV with B-factor data for paratope."""
    rows = []
    for chain in ['A', 'B']:
        if chain in res_dict:
            for pos in sorted(res_dict[chain]):
                if chain in bfactors and pos in bfactors[chain]:
                    data = bfactors[chain][pos]
                    is_conserved = chain in conserved_residues and pos in conserved_residues[chain]
                    
                    # Classify flexibility
                    if data['avg_bfactor'] > 50:
                        flex_class = 'High'
                    elif data['avg_bfactor'] >= 30:
                        flex_class = 'Medium'
                    else:
                        flex_class = 'Low'
                    
                    rows.append({
                        'Chain': chain,
                        'Position': pos,
                        'ResName': data['resname'],
                        'Avg_Bfactor': round(data['avg_bfactor'], 2),
                        'Std_Bfactor': round(data['std_bfactor'], 2),
                        'Min_Bfactor': round(data['min_bfactor'], 2),
                        'Max_Bfactor': round(data['max_bfactor'], 2),
                        'Flexibility': flex_class,
                        'In_Conserved': is_conserved
                    })
    
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['Chain', 'Position', 'ResName', 'Avg_Bfactor', 'Std_Bfactor', 
                      'Min_Bfactor', 'Max_Bfactor', 'Flexibility', 'In_Conserved']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    
    print(f"\n✓ Generated {output_file} with {len(rows)} residues")

def main():
    print("="*60)
    print("B-FACTOR FLEXIBILITY ANALYSIS")
    print("="*60)
    
    # Parse both structures
    print("\nParsing CL1-Ab4_p53.pdb...")
    bfactors_p53 = parse_bfactors('/Users/aloksinha/Desktop/CL1-Ab4_p53.pdb')
    
    print("Parsing CL1-Ab4_CALR.pdb...")
    bfactors_calr = parse_bfactors('/Users/aloksinha/Desktop/CL1-Ab4_CALR.pdb')
    
    # Analyze B-factor distributions
    print("\n" + "="*60)
    print("p53 COMPLEX")
    print("="*60)
    
    para_stats_p53 = analyze_bfactor_distribution(bfactors_p53, paratope_residues, "Full Paratope")
    cons_stats_p53 = analyze_bfactor_distribution(bfactors_p53, conserved_residues, "Conserved Interface")
    
    print("\n" + "="*60)
    print("CALR COMPLEX")
    print("="*60)
    
    para_stats_calr = analyze_bfactor_distribution(bfactors_calr, paratope_residues, "Full Paratope")
    cons_stats_calr = analyze_bfactor_distribution(bfactors_calr, conserved_residues, "Conserved Interface")
    
    # Create CSV files
    create_bfactor_csv(bfactors_p53, paratope_residues, '/Users/aloksinha/Desktop/bfactors_p53.csv')
    create_bfactor_csv(bfactors_calr, paratope_residues, '/Users/aloksinha/Desktop/bfactors_calr.csv')
    
    # Key findings
    print("\n" + "="*60)
    print("KEY FINDINGS")
    print("="*60)
    
    if cons_stats_p53 and cons_stats_calr:
        avg_flex = (cons_stats_p53['mean'] + cons_stats_calr['mean']) / 2
        
        print(f"\n1. OVERALL FLEXIBILITY:")
        print(f"   Average conserved interface B-factor: {avg_flex:.2f} Ų")
        
        if avg_flex > 50:
            print("   → HIGH FLEXIBILITY: Suggests induced-fit binding mechanism")
        elif avg_flex > 30:
            print("   → MODERATE FLEXIBILITY: Some conformational adaptability")
        else:
            print("   → LOW FLEXIBILITY: Rigid interface, less likely induced-fit")
        
        high_flex_frac = (cons_stats_p53['high_flex_count'] + cons_stats_calr['high_flex_count']) / \
                         (len([p for c in conserved_residues.values() for p in c]) * 2)
        
        print(f"\n2. FLEXIBLE HOTSPOTS:")
        print(f"   Fraction of highly flexible residues (B>50): {100*high_flex_frac:.1f}%")
        
        if high_flex_frac > 0.3:
            print("   → Many flexible residues suggest conformational promiscuity")
        
        print("\n3. POLYREACTIVITY MECHANISM:")
        if avg_flex > 40:
            print("   ✓ High flexibility supports INDUCED-FIT polyreactivity")
        else:
            print("   ✓ Lower flexibility suggests RIGID/CHEMICAL-PATCH polyreactivity")
    
    # Create PyMOL visualization script
    with open('/Users/aloksinha/Desktop/visualize_bfactors.pml', 'w') as f:
        f.write("""# PyMOL script to visualize B-factors
load /Users/aloksinha/Desktop/CL1-Ab4_p53.pdb, p53_complex
load /Users/aloksinha/Desktop/CL1-Ab4_CALR.pdb, calr_complex

# Color by B-factor (blue=low, white=medium, red=high)
spectrum b, blue_white_red, minimum=20, maximum=80

# Show paratope as sticks
select paratope, chain A and resi 31+32+33+50+52+55+57+59+99+102+103+104+105+106+107+109 or chain B and resi 2+27+31+98+99+100+102
show sticks, paratope

# Surface colored by B-factor
show surface, chain A+B
set transparency, 0.3

# Label highly flexible residues (B>50)
select flexible, paratope and b > 50
label flexible and name CA, "%s%s" % (resn,resi)

print "Visualization complete. Red regions = high B-factor (flexible)"
print "Blue regions = low B-factor (rigid)"
""")
    
    print("\n✓ Generated visualize_bfactors.pml for PyMOL visualization")
    print("\nTo visualize in PyMOL:")
    print("  @/Users/aloksinha/Desktop/visualize_bfactors.pml")

if __name__ == '__main__':
    main()

