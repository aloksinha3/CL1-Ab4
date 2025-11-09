#!/usr/bin/env python3
"""Extract antibody sequences from PDB and prepare for ANARCI numbering."""

def extract_sequence_from_pdb(pdb_file, chain_id):
    """Extract amino acid sequence from PDB ATOM records for a specific chain."""
    residues = {}
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                chain = line[21].strip()
                if chain == chain_id:
                    resnum = int(line[22:26].strip())
                    resname = line[17:20].strip()
                    if resnum not in residues:
                        residues[resnum] = resname
    
    # Convert 3-letter to 1-letter code
    aa_code = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
        'H1S': 'H'  # Special case for HIS variant
    }
    
    sequence = ''.join(aa_code.get(residues[i], 'X') for i in sorted(residues.keys()))
    return sequence

def main():
    pdb_file = '/Users/aloksinha/Desktop/CL1-Ab4_p53.pdb'
    
    # Extract heavy chain (A) sequence
    heavy_seq = extract_sequence_from_pdb(pdb_file, 'A')
    print(f"Heavy chain (A) sequence ({len(heavy_seq)} aa):")
    print(heavy_seq)
    
    with open('/Users/aloksinha/Desktop/heavy.fasta', 'w') as f:
        f.write('>CL1-Ab4_Heavy_Chain_A\n')
        # Write in 60-character lines
        for i in range(0, len(heavy_seq), 60):
            f.write(heavy_seq[i:i+60] + '\n')
    
    print("\n✓ Saved to heavy.fasta")
    
    # Extract light chain (B) sequence
    light_seq = extract_sequence_from_pdb(pdb_file, 'B')
    print(f"\nLight chain (B) sequence ({len(light_seq)} aa):")
    print(light_seq)
    
    with open('/Users/aloksinha/Desktop/light.fasta', 'w') as f:
        f.write('>CL1-Ab4_Light_Chain_B\n')
        for i in range(0, len(light_seq), 60):
            f.write(light_seq[i:i+60] + '\n')
    
    print("✓ Saved to light.fasta")
    
    print("\n" + "="*60)
    print("ANARCI NUMBERING INSTRUCTIONS")
    print("="*60)
    
    print("\n1. Install ANARCI (if not already installed):")
    print("   pip install anarci")
    print("   # or")
    print("   conda install -c bioconda anarci")
    
    print("\n2. Run ANARCI numbering:")
    print("   cd /Users/aloksinha/Desktop")
    print("   ANARCI -i heavy.fasta -s kabat -o heavy_kabat.txt")
    print("   ANARCI -i light.fasta -s kabat -o light_kabat.txt")
    
    print("\n3. Alternative: IMGT numbering (recommended for polyreactivity analysis)")
    print("   ANARCI -i heavy.fasta -s imgt -o heavy_imgt.txt")
    print("   ANARCI -i light.fasta -s imgt -o light_imgt.txt")
    
    print("\n4. Or use online tool:")
    print("   - Go to: http://opig.stats.ox.ac.uk/webapps/newsabdab/sabpred/anarci/")
    print("   - Upload heavy.fasta and light.fasta")
    print("   - Select Kabat or IMGT scheme")
    
    print("\n" + "="*60)
    print("CDR IDENTIFICATION")
    print("="*60)
    
    print("\nBased on conserved interface residues, expected CDR locations:")
    print("\nHeavy chain conserved residues:")
    print("  A31-33 → likely HCDR1")
    print("  A50-59 → likely HCDR2")
    print("  A99-110 → likely HCDR3 (long, flexible)")
    
    print("\nLight chain conserved residues:")
    print("  B27-38 → likely LCDR1")
    print("  B55-56 → likely LCDR2")
    print("  B96-102 → likely LCDR3")
    
    print("\nHCDR3 length: ~12 residues (A99-110)")
    print("  → Long HCDR3 correlates with polyreactivity ✓")

if __name__ == '__main__':
    main()

