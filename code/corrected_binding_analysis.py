#!/usr/bin/env python3
"""Corrected binding mode analysis with proper antibody superposition."""

import numpy as np
from scipy.spatial.transform import Rotation
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

def kabsch_superposition(coords1, coords2, common_residues):
    """Perform Kabsch superposition to align two sets of coordinates."""
    # Get common residues
    points1 = []
    points2 = []
    for res in common_residues:
        if res in coords1 and res in coords2:
            points1.append(coords1[res])
            points2.append(coords2[res])
    
    if len(points1) < 3:
        return None, None, None
    
    points1 = np.array(points1)
    points2 = np.array(points2)
    
    # Center coordinates
    centroid1 = np.mean(points1, axis=0)
    centroid2 = np.mean(points2, axis=0)
    
    centered1 = points1 - centroid1
    centered2 = points2 - centroid2
    
    # Calculate rotation matrix using SVD
    H = centered1.T @ centered2
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    
    # Ensure proper rotation (det(R) = 1)
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    
    # Calculate translation
    t = centroid2 - R @ centroid1
    
    return R, t, centroid1

def apply_transformation(coords, R, t):
    """Apply rotation and translation to coordinates."""
    transformed = {}
    for res, coord in coords.items():
        transformed[res] = R @ coord + t
    return transformed

def calculate_rmsd(coords1, coords2, residues):
    """Calculate RMSD between two sets of coordinates."""
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
    
    diff = points1 - points2
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    
    return rmsd

def calculate_binding_orientation(ab_coords, ag_coords, ab_residues, ag_residues):
    """Calculate binding orientation vector."""
    ab_points = [ab_coords[res] for res in ab_residues if res in ab_coords]
    ag_points = [ag_coords[res] for res in ag_residues if res in ag_coords]
    
    if not ab_points or not ag_points:
        return None
    
    ab_centroid = np.mean(ab_points, axis=0)
    ag_centroid = np.mean(ag_points, axis=0)
    
    vector = ag_centroid - ab_centroid
    return vector / np.linalg.norm(vector)

def main():
    print("="*60)
    print("CORRECTED BINDING MODE ANALYSIS")
    print("="*60)
    
    # Parse both complexes
    print("\nParsing structures...")
    p53_coords = parse_pdb_coords('/Users/aloksinha/Desktop/CL1-Ab4_p53.pdb', ['A', 'B', 'C'])
    calr_coords = parse_pdb_coords('/Users/aloksinha/Desktop/CL1-Ab4_CALR.pdb', ['A', 'B', 'C'])
    
    # Conserved antibody interface residues (for superposition)
    ab_conserved_A = [31, 32, 33, 50, 52, 55, 57, 59, 99, 102, 103, 104, 105, 106, 107, 109]
    ab_conserved_B = [2, 27, 31, 98, 99, 100, 102]
    all_ab_conserved = ab_conserved_A + ab_conserved_B
    
    # Antigen interface residues
    p53_interface = [1, 2, 3, 6, 10, 63, 65, 67, 111, 112, 113, 114, 115, 116, 117, 118, 
                     120, 122, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174]
    calr_interface = [25, 27, 30, 31, 32, 33, 34, 35, 36, 37, 39, 80, 81, 83, 84, 85, 86, 87, 88, 89]
    
    # Perform superposition using conserved antibody residues
    print("\nPerforming Kabsch superposition on conserved antibody residues...")
    
    # Combine heavy and light chain coordinates for superposition
    ab_p53_combined = {**p53_coords['A'], **{f"B{r}": p53_coords['B'][r] for r in p53_coords['B']}}
    ab_calr_combined = {**calr_coords['A'], **{f"B{r}": calr_coords['B'][r] for r in calr_coords['B']}}
    
    common_residues = [r for r in all_ab_conserved if r in ab_p53_combined and r in ab_calr_combined]
    common_residues.extend([f"B{r}" for r in ab_conserved_B if r in p53_coords['B'] and r in calr_coords['B']])
    
    R, t, centroid = kabsch_superposition(ab_p53_combined, ab_calr_combined, common_residues)
    
    if R is None:
        print("Error: Could not perform superposition")
        return
    
    print(f"Superposition performed on {len(common_residues)} conserved residues")
    
    # Apply transformation to CALR complex
    calr_coords_aligned = {}
    for chain in ['A', 'B', 'C']:
        calr_coords_aligned[chain] = {}
        for res, coord in calr_coords[chain].items():
            if chain in ['A', 'B']:  # Transform antibody
                calr_coords_aligned[chain][res] = R @ coord + t
            else:  # Don't transform antigen
                calr_coords_aligned[chain][res] = coord
    
    # Calculate corrected RMSDs
    print("\n" + "="*60)
    print("CORRECTED RMSD ANALYSIS")
    print("="*60)
    
    rmsd_heavy = calculate_rmsd(p53_coords['A'], calr_coords_aligned['A'], ab_conserved_A)
    rmsd_light = calculate_rmsd(p53_coords['B'], calr_coords_aligned['B'], ab_conserved_B)
    rmsd_all = calculate_rmsd({**p53_coords['A'], **{f"B{r}": p53_coords['B'][r] for r in p53_coords['B']}}, 
                             {**calr_coords_aligned['A'], **{f"B{r}": calr_coords_aligned['B'][r] for r in calr_coords_aligned['B']}}, 
                             common_residues)
    
    print(f"\nHeavy chain (conserved interface) RMSD: {rmsd_heavy:.3f} Å")
    print(f"Light chain (conserved interface) RMSD: {rmsd_light:.3f} Å")
    print(f"Combined conserved interface RMSD: {rmsd_all:.3f} Å")
    
    # Interpretation
    print("\n" + "="*60)
    print("INTERPRETATION")
    print("="*60)
    
    if rmsd_heavy and rmsd_heavy < 1.0:
        print("✓ Heavy chain: Very similar conformation (RMSD < 1 Å)")
        print("  → Chemical promiscuity, not structural adaptation")
    elif rmsd_heavy and rmsd_heavy < 2.0:
        print("✓ Heavy chain: Similar conformation (RMSD < 2 Å)")
        print("  → Minor induced-fit, mainly chemical promiscuity")
    elif rmsd_heavy and rmsd_heavy < 3.0:
        print("✓ Heavy chain: Moderate conformational change (RMSD < 3 Å)")
        print("  → Some induced-fit mechanism")
    else:
        print("✓ Heavy chain: Significant conformational change")
        print("  → Substantial induced-fit mechanism")
    
    # Binding orientation analysis (with aligned structures)
    print("\n" + "="*60)
    print("BINDING ORIENTATION ANALYSIS")
    print("="*60)
    
    p53_vector = calculate_binding_orientation(p53_coords['A'], p53_coords['C'], 
                                                ab_conserved_A, p53_interface)
    calr_vector = calculate_binding_orientation(calr_coords_aligned['A'], calr_coords['C'],
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
    
    # Update summary CSV
    print("\n" + "="*60)
    print("UPDATING SUMMARY")
    print("="*60)
    
    summary = {
        'Metric': ['Heavy_RMSD_Corrected', 'Light_RMSD_Corrected', 'Combined_RMSD_Corrected', 
                   'Binding_Angle_Corrected', 'Superposition_Residues'],
        'Value': [rmsd_heavy, rmsd_light, rmsd_all, angle if p53_vector is not None else None, 
                  len(common_residues)],
        'Unit': ['Å', 'Å', 'Å', '°', 'residues']
    }
    
    with open('/Users/aloksinha/Desktop/binding_mode_comparison_corrected.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Metric', 'Value', 'Unit'])
        for i in range(len(summary['Metric'])):
            writer.writerow([summary['Metric'][i], summary['Value'][i], summary['Unit'][i]])
    
    print("\n✓ Generated binding_mode_comparison_corrected.csv")
    
    # Final interpretation
    print("\n" + "="*60)
    print("FINAL POLYREACTIVITY MECHANISM")
    print("="*60)
    
    if rmsd_all and rmsd_all < 1.0:
        print("**MECHANISM: Chemical Promiscuity**")
        print("- Same paratope conformation for both antigens")
        print("- Versatile chemistry (aromatic/hydrophobic patch) enables diverse binding")
        print("- No significant structural adaptation required")
    elif rmsd_all and rmsd_all < 2.0:
        print("**MECHANISM: Chemical Promiscuity + Minor Induced-Fit**")
        print("- Similar paratope conformation with small adjustments")
        print("- Chemical versatility is primary driver")
        print("- Some conformational fine-tuning for each antigen")
    else:
        print("**MECHANISM: Induced-Fit Flexibility**")
        print("- Significant conformational adaptation for each antigen")
        print("- High flexibility (B-factors) enables structural promiscuity")
        print("- Chemical versatility complements structural adaptation")

if __name__ == '__main__':
    main()

