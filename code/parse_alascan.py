#!/usr/bin/env python3
"""Parse FoldX AlaScan output and generate ΔΔG analysis for antibody interface residues."""

import csv
import re

# Antibody interface residues (from previous contact analysis)
ab_residues = {
    'A': [31, 32, 33, 50, 52, 55, 57, 59, 99, 102, 103, 104, 105, 106, 107, 109,
          2, 27, 28, 30, 58, 98, 100, 101, 108, 110],  # All heavy chain residues
    'B': [2, 27, 31, 98, 99, 100, 102, 1, 26, 28, 29, 30, 32, 33, 34, 35, 36, 38,
          55, 56, 96, 97, 101]  # All light chain residues
}

def parse_alascan(filepath):
    """Parse FoldX AlaScan output file."""
    results = {}
    with open(filepath, 'r') as f:
        for line in f:
            # Format: "RES NUM to ALA energy change is VALUE"
            match = re.match(r'([A-Z0-9]{3})\s+(\d+)\s+to ALA energy change is\s+([-\d.]+)', line)
            if match:
                res_name = match.group(1)
                res_num = int(match.group(2))
                ddg = float(match.group(3))
                results[res_num] = {'res_name': res_name, 'ddg': ddg}
    return results

def main():
    p53_data = parse_alascan('/Users/aloksinha/Desktop/foldx_out/CL1-Ab4_p53_Repair_AS.fxout')
    calr_data = parse_alascan('/Users/aloksinha/Desktop/foldx_out/CL1-Ab4_CALR_Repair_AS.fxout')
    
    # Create comprehensive CSV
    rows = []
    
    # Process heavy chain (A)
    for res_num in sorted(set(ab_residues['A'])):
        if res_num in p53_data or res_num in calr_data:
            p53_info = p53_data.get(res_num, {'res_name': '?', 'ddg': None})
            calr_info = calr_data.get(res_num, {'res_name': '?', 'ddg': None})
            
            rows.append({
                'Chain': 'A',
                'Position': res_num,
                'ResName': p53_info['res_name'],
                'ddG_p53': p53_info['ddg'],
                'ddG_CALR': calr_info['ddg'],
                'ddG_avg': (p53_info['ddg'] + calr_info['ddg']) / 2 if p53_info['ddg'] is not None and calr_info['ddg'] is not None else None,
                'ddG_diff': abs(p53_info['ddg'] - calr_info['ddg']) if p53_info['ddg'] is not None and calr_info['ddg'] is not None else None
            })
    
    # Process light chain (B) - need to find chain B data in the files
    # Note: The AlaScan output shows all residues sequentially, so we need to identify chain B residues
    # For now, let's focus on chain A (heavy) which has most of the conserved contacts
    
    # Write CSV
    with open('/Users/aloksinha/Desktop/foldx_ddg_analysis.csv', 'w', newline='') as csvfile:
        fieldnames = ['Chain', 'Position', 'ResName', 'ddG_p53', 'ddG_CALR', 'ddG_avg', 'ddG_diff']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    
    print(f"✓ Generated foldx_ddg_analysis.csv with {len(rows)} residues")
    
    # Print summary statistics
    print("\n=== ΔΔG Summary (kcal/mol) ===")
    print(f"{'Residue':<12} {'p53':<8} {'CALR':<8} {'Avg':<8} {'Diff':<8}")
    print("-" * 50)
    for row in rows:
        if row['ddG_p53'] is not None:
            print(f"{row['Chain']}{row['Position']:<11} {row['ddG_p53']:<8.2f} {row['ddG_CALR']:<8.2f} {row['ddG_avg']:<8.2f} {row['ddG_diff']:<8.2f}")
    
    # Identify key findings
    print("\n=== Key Findings ===")
    conserved_residues = [r for r in rows if r['ddG_avg'] is not None and r['ddG_avg'] > 1.0 and r['ddG_diff'] < 2.0]
    print(f"Conserved hotspots (ΔΔG > 1.0, similar in both): {len(conserved_residues)}")
    for r in sorted(conserved_residues, key=lambda x: x['ddG_avg'], reverse=True)[:10]:
        print(f"  {r['Chain']}{r['Position']} ({r['ResName']}): p53={r['ddG_p53']:.2f}, CALR={r['ddG_CALR']:.2f}")
    
    differential_residues = [r for r in rows if r['ddG_diff'] is not None and r['ddG_diff'] > 2.0]
    print(f"\nDifferential residues (ΔΔG diff > 2.0): {len(differential_residues)}")
    for r in sorted(differential_residues, key=lambda x: x['ddG_diff'], reverse=True)[:5]:
        print(f"  {r['Chain']}{r['Position']} ({r['ResName']}): p53={r['ddG_p53']:.2f}, CALR={r['ddG_CALR']:.2f}, diff={r['ddG_diff']:.2f}")

if __name__ == '__main__':
    main()

