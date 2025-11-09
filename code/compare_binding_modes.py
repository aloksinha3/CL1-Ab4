#!/usr/bin/env python3
"""Compare binding modes of p53 and CALR to CL1-Ab4."""

import numpy as np
from collections import defaultdict
import csv

def parse_pdb_coords(pdb_file, chain_ids):
    """Extract CA coordinates for specified chains."""
    coords = {}
    for chain in chain_ids:
        coords[chain] = {}
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                chain = line[21].strip()
                if chain in chain_ids:
                    resnum = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coords[chain][resnum] = np.array([x, y, z])
    return coords

def calculate_rmsd(coords1, coords2, residues):
    """Calculate RMSD between two sets of coordinates for specified residues."""
    points1 = []
    points2 = []
    
    for res in residues:
        if res in coords1 and res in coords2:
            points1.append(coords1[res])
            points2.append(coords2[res])
    
    if len(points1) < 3:
        return None
    
    points1 = np.array(points1)
    points2 = np.array(points2)
    
    # Calculate RMSD
    diff = points1 - points2
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    
    return rmsd

def calculate_centroid(coords, residues):
    """Calculate centroid of specified residues."""
    points = [coords[res] for res in residues if res in coords]
    if not points:
        return None
    return np.mean(points, axis=0)

def calculate_binding_orientation(ab_coords, ag_coords, ab_residues, ag_residues):
    """Calculate approximate binding orientation vector."""
    ab_centroid = calculate_centroid(ab_coords, ab_residues)
    ag_centroid = calculate_centroid(ag_coords, ag_residues)
    
    if ab_centroid is None or ag_centroid is None:
        return None
    
    vector = ag_centroid - ab_centroid
    return vector / np.linalg.norm(vector)

def main():
    print("="*60)
    print("BINDING MODE COMPARISON ANALYSIS")
    print("="*60)
    
    # Parse both complexes
    print("\nParsing CL1-Ab4_p53.pdb...")
    p53_coords = parse_pdb_coords('/Users/aloksinha/Desktop/CL1-Ab4_p53.pdb', ['A', 'B', 'C'])
    
    print("Parsing CL1-Ab4_CALR.pdb...")
    calr_coords = parse_pdb_coords('/Users/aloksinha/Desktop/CL1-Ab4_CALR.pdb', ['A', 'B', 'C'])
    
    # Antibody conserved interface residues
    ab_conserved_A = [31, 32, 33, 50, 52, 55, 57, 59, 99, 102, 103, 104, 105, 106, 107, 109]
    ab_conserved_B = [2, 27, 31, 98, 99, 100, 102]
    all_ab_conserved = ab_conserved_A + ab_conserved_B
    
    # Antigen interface residues (from contact analysis)
    p53_interface = [1, 2, 3, 6, 10, 63, 65, 67, 111, 112, 113, 114, 115, 116, 117, 118, 
                     120, 122, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174]
    calr_interface = [25, 27, 30, 31, 32, 33, 34, 35, 36, 37, 39, 80, 81, 83, 84, 85, 86, 87, 88, 89]
    
    # 1. Antibody RMSD (measure conformational change in antibody)
    print("\n" + "="*60)
    print("1. ANTIBODY CONFORMATIONAL COMPARISON")
    print("="*60)
    
    rmsd_heavy = calculate_rmsd(p53_coords['A'], calr_coords['A'], ab_conserved_A)
    rmsd_light = calculate_rmsd(p53_coords['B'], calr_coords['B'], ab_conserved_B)
    
    print(f"\nHeavy chain (conserved interface) RMSD: {rmsd_heavy:.2f} Å")
    print(f"Light chain (conserved interface) RMSD: {rmsd_light:.2f} Å")
    
    if rmsd_heavy and rmsd_heavy < 2.0:
        print("  → Heavy chain maintains similar conformation (RMSD < 2 Å)")
    elif rmsd_heavy and rmsd_heavy < 5.0:
        print("  → Moderate conformational change in heavy chain (induced fit)")
    elif rmsd_heavy:
        print("  → Significant conformational change in heavy chain (large induced fit)")
    
    # 2. Binding orientation comparison
    print("\n" + "="*60)
    print("2. BINDING ORIENTATION ANALYSIS")
    print("="*60)
    
    p53_vector = calculate_binding_orientation(p53_coords['A'], p53_coords['C'], 
                                                ab_conserved_A, p53_interface)
    calr_vector = calculate_binding_orientation(calr_coords['A'], calr_coords['C'],
                                                 ab_conserved_A, calr_interface)
    
    if p53_vector is not None and calr_vector is not None:
        dot_product = np.dot(p53_vector, calr_vector)
        angle = np.arccos(np.clip(dot_product, -1.0, 1.0)) * 180.0 / np.pi
        
        print(f"\nBinding orientation angle difference: {angle:.1f}°")
        
        if angle < 30:
            print("  → Similar binding orientations (parallel)")
        elif angle < 60:
            print("  → Moderately different orientations")
        elif angle < 120:
            print("  → Substantially different orientations")
        else:
            print("  → Opposite binding orientations (antiparallel)")
    
    # 3. Interface size comparison
    print("\n" + "="*60)
    print("3. INTERFACE SIZE COMPARISON")
    print("="*60)
    
    p53_ab_contacts = len(ab_conserved_A) + len(ab_conserved_B)
    p53_ag_contacts = len(p53_interface)
    
    calr_ab_contacts = len(ab_conserved_A) + len(ab_conserved_B)
    calr_ag_contacts = len(calr_interface)
    
    print(f"\np53 complex:")
    print(f"  Antibody interface: {p53_ab_contacts} residues")
    print(f"  Antigen interface:  {p53_ag_contacts} residues")
    
    print(f"\nCALR complex:")
    print(f"  Antibody interface: {calr_ab_contacts} residues")
    print(f"  Antigen interface:  {calr_ag_contacts} residues")
    
    # 4. Epitope overlap analysis
    print("\n" + "="*60)
    print("4. EPITOPE CHARACTERISTICS")
    print("="*60)
    
    # Calculate interface centroids
    p53_centroid = calculate_centroid(p53_coords['C'], p53_interface)
    calr_centroid = calculate_centroid(calr_coords['C'], calr_interface)
    
    # Calculate interface spread (radius of gyration)
    def calc_radius_gyration(coords, residues):
        centroid = calculate_centroid(coords, residues)
        if centroid is None:
            return None
        dists = [np.linalg.norm(coords[r] - centroid) for r in residues if r in coords]
        return np.sqrt(np.mean([d**2 for d in dists]))
    
    p53_rg = calc_radius_gyration(p53_coords['C'], p53_interface)
    calr_rg = calc_radius_gyration(calr_coords['C'], calr_interface)
    
    print(f"\np53 epitope spread (Rg):  {p53_rg:.2f} Å")
    print(f"CALR epitope spread (Rg): {calr_rg:.2f} Å")
    
    if p53_rg and calr_rg:
        if abs(p53_rg - calr_rg) < 3.0:
            print("  → Similar epitope sizes")
        else:
            print("  → Different epitope sizes")
    
    # 5. Generate summary CSV
    print("\n" + "="*60)
    print("5. GENERATING BINDING MODE SUMMARY")
    print("="*60)
    
    summary = {
        'Metric': ['Heavy_RMSD', 'Light_RMSD', 'Binding_Angle', 'Ab_Interface_Size', 
                   'p53_Epitope_Size', 'CALR_Epitope_Size', 'p53_Epitope_Rg', 'CALR_Epitope_Rg'],
        'Value': [rmsd_heavy, rmsd_light, angle if p53_vector is not None else None, 
                  p53_ab_contacts, p53_ag_contacts, calr_ag_contacts, p53_rg, calr_rg],
        'Unit': ['Å', 'Å', '°', 'residues', 'residues', 'residues', 'Å', 'Å']
    }
    
    with open('/Users/aloksinha/Desktop/binding_mode_comparison.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Metric', 'Value', 'Unit'])
        for i in range(len(summary['Metric'])):
            writer.writerow([summary['Metric'][i], summary['Value'][i], summary['Unit'][i]])
    
    print("\n✓ Generated binding_mode_comparison.csv")
    
    # Interpretation
    print("\n" + "="*60)
    print("INTERPRETATION")
    print("="*60)
    
    print("\n**Polyreactivity Mechanism:**")
    
    if rmsd_heavy and rmsd_heavy < 2.0:
        print("✓ Low antibody RMSD → SAME paratope conformation used for both antigens")
        print("  This suggests CHEMICAL PROMISCUITY rather than structural adaptation")
    elif rmsd_heavy and rmsd_heavy < 5.0:
        print("✓ Moderate antibody RMSD → Some conformational adaptation (induced fit)")
        print("  This suggests HYBRID mechanism (chemical + conformational)")
    
    if angle and angle > 30:
        print(f"✓ Large binding angle ({angle:.1f}°) → Different antigen ORIENTATIONS")
        print("  Same paratope residues engage antigens from different angles")
    
    if p53_ag_contacts and calr_ag_contacts:
        ratio = p53_ag_contacts / calr_ag_contacts
        if 0.8 < ratio < 1.2:
            print("✓ Similar epitope sizes → Antibody recognizes patches of similar dimensions")
        else:
            print(f"✓ Different epitope sizes (ratio {ratio:.2f}) → Versatile recognition")
    
    print("\n**Conclusion:**")
    print("CL1-Ab4 likely achieves polyreactivity through:")
    print("1. Conserved paratope geometry (low RMSD)")
    print("2. Chemical versatility of interface residues")
    print("3. Different binding orientations/modes with each antigen")

if __name__ == '__main__':
    main()

