# PyMOL script to visualize B-factors
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
