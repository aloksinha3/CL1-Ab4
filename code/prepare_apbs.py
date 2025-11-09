#!/usr/bin/env python3
"""Prepare antibody structure for APBS electrostatic analysis."""

def extract_antibody(input_pdb, output_pdb):
    """Extract antibody chains (A and B) from complex PDB."""
    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        for line in infile:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chain = line[21].strip()
                if chain in ['A', 'B']:
                    outfile.write(line)
            elif line.startswith('TER'):
                chain = line[21].strip() if len(line) > 21 else ''
                if chain in ['A', 'B']:
                    outfile.write(line)
            elif line.startswith(('HEADER', 'CRYST1', 'END')):
                outfile.write(line)
    print(f"✓ Extracted antibody chains A+B to {output_pdb}")

def create_apbs_input(pqr_file, output_file='apbs.in'):
    """Create APBS input file for electrostatic calculation."""
    apbs_config = f"""read
    mol pqr {pqr_file}
end

elec name antibody
    mg-auto
    dime 129 129 129
    cglen 150 150 150
    fglen 100 100 100
    cgcent mol 1
    fgcent mol 1
    mol 1
    npbe
    bcfl sdh
    ion charge 1 conc 0.150 radius 2.0
    ion charge -1 conc 0.150 radius 1.8
    pdie 2.0
    sdie 78.54
    srfm smol
    chgm spl2
    sdens 10.0
    srad 1.4
    swin 0.3
    temp 298.15
    calcenergy total
    calcforce no
    write pot dx antibody
end

quit
"""
    with open(output_file, 'w') as f:
        f.write(apbs_config)
    print(f"✓ Created APBS input file: {output_file}")
    print("\nTo run APBS electrostatics:")
    print(f"1. pdb2pqr --ff=AMBER --with-ph=7.0 {pqr_file.replace('.pqr', '.pdb')} {pqr_file}")
    print(f"2. apbs {output_file}")
    print("3. Load .dx file in PyMOL or Chimera for visualization")

def main():
    # Extract antibody from both complexes
    extract_antibody('/Users/aloksinha/Desktop/CL1-Ab4_p53.pdb', 
                     '/Users/aloksinha/Desktop/CL1-Ab4_antibody_from_p53.pdb')
    
    extract_antibody('/Users/aloksinha/Desktop/CL1-Ab4_CALR.pdb',
                     '/Users/aloksinha/Desktop/CL1-Ab4_antibody_from_calr.pdb')
    
    # Create APBS input file
    create_apbs_input('CL1-Ab4_antibody_from_p53.pqr', 
                      '/Users/aloksinha/Desktop/apbs_antibody.in')
    
    print("\n" + "="*60)
    print("APBS ELECTROSTATICS WORKFLOW")
    print("="*60)
    print("\nOption 1: Use PyMOL APBS plugin (easiest)")
    print("  - Load CL1-Ab4_antibody_from_p53.pdb in PyMOL")
    print("  - Plugin → APBS Electrostatics")
    print("  - Click 'Set grid' and 'Run APBS'")
    
    print("\nOption 2: Command line (requires pdb2pqr and apbs installed)")
    print("  cd /Users/aloksinha/Desktop")
    print("  pdb2pqr --ff=AMBER --with-ph=7.0 CL1-Ab4_antibody_from_p53.pdb CL1-Ab4_antibody_from_p53.pqr")
    print("  apbs apbs_antibody.in")
    
    print("\nOption 3: Use web server")
    print("  - Go to: https://server.poissonboltzmann.org/")
    print("  - Upload: CL1-Ab4_antibody_from_p53.pdb")
    print("  - Download results and visualize")
    
    print("\n" + "="*60)
    print("HYDROPHOBIC PATCH VISUALIZATION (PyMOL)")
    print("="*60)
    print("\nRun in PyMOL:")
    print("  @/Users/aloksinha/Desktop/electrostatics_hydrophobic.pml")
    print("\nOr manually:")
    print("  load CL1-Ab4_antibody_from_p53.pdb")
    print("  select hydro, resn PHE+TYR+TRP+LEU+ILE+VAL+ALA+MET")
    print("  select paratope_hydro, hydro and (resi 31+32+33+50+52+55+57+59+99+102+103+104+105+106+107+109)")
    print("  show surface, paratope_hydro")
    print("  color yellow, paratope_hydro")

if __name__ == '__main__':
    main()

