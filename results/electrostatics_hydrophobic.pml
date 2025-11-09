# PyMOL script for electrostatic potential and hydrophobic patch analysis
# For CL1-Ab4 antibody paratope

# Load structures
load /Users/aloksinha/Desktop/CL1-Ab4_p53.pdb, p53_complex
load /Users/aloksinha/Desktop/CL1-Ab4_CALR.pdb, calr_complex

# Extract antibody only (chains A and B)
create antibody_p53, p53_complex and chain A+B
create antibody_calr, calr_complex and chain A+B

# Save antibody-only PDB for APBS
save /Users/aloksinha/Desktop/CL1-Ab4_antibody_only.pdb, antibody_p53

# Define paratope residues (conserved interface residues)
select paratope, chain A and resi 31+32+33+50+52+55+57+59+99+102+103+104+105+106+107+109+2+27+28+30+58+98+100+101+108+110 or chain B and resi 2+27+31+98+99+100+102+1+26+28+29+30+32+33+34+35+36+38+55+56+96+97+101

# Hydrophobic residue analysis
# Define hydrophobic residues (F, Y, W, L, I, V, A, M)
select hydrophobic, resn PHE+TYR+TRP+LEU+ILE+VAL+ALA+MET

# Hydrophobic residues in paratope
select paratope_hydro, paratope and hydrophobic

# Aromatic residues in paratope  
select aromatic, paratope and resn PHE+TYR+TRP+HIS

# Basic residues in paratope
select basic, paratope and resn ARG+LYS

# Acidic residues in paratope
select acidic, paratope and resn ASP+GLU

# Visualize
hide everything
show cartoon, antibody_p53
color gray80, antibody_p53

# Color by chemistry
show sticks, paratope
color yellow, paratope_hydro
color orange, aromatic
color blue, basic
color red, acidic

# Surface representation
show surface, antibody_p53
set transparency, 0.4

# APBS electrostatics (requires APBS plugin)
# Run these commands manually or with APBS plugin:
# 1. Select antibody structure
# 2. Tools > APBS Electrostatics
# 3. Configure settings and run

# For manual APBS:
# pdb2pqr --ff=PARSE --with-ph=7.0 CL1-Ab4_antibody_only.pdb CL1-Ab4_antibody.pqr
# apbs apbs.in (need to create apbs.in config file)

# Export hydrophobic patch data
select hydro_surface, antibody_p53 and paratope_hydro and (solvent expand 1.4)

# Print residue counts
print "=== Paratope Composition ==="
count_atoms paratope
count_atoms paratope_hydro  
count_atoms aromatic
count_atoms basic
count_atoms acidic

# Save session
save /Users/aloksinha/Desktop/electrostatics_session.pse

print "Analysis complete. For APBS electrostatics:"
print "1. Use APBS Tools plugin in PyMOL"
print "2. Or run: pdb2pqr --ff=PARSE --with-ph=7.0 CL1-Ab4_antibody_only.pdb CL1-Ab4_antibody.pqr"
print "3. Then: apbs apbs.in"

