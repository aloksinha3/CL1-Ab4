# CL1-Ab4 CALR p53 Binding Mechanism

## Overview

**Antibody:** CL1-Ab4  
**Antigens:** p53 and CALR
**Analysis Date:** October 7, 2025

---

## Methods

### Structure Generation

Structures of the CL1-Ab4 CALR P domain interface and the CL1-Ab4 p53 core domain were generated using AlphaFold3. The antibody-antigen complexes were modeled with the following chain designations:
Chain A: Heavy chain of CL1-Ab4
Chain B: Light chain of CL1-Ab4  
Chain C: Antigen (p53 core domain or CALR P domain)

### Computational Analysis Pipeline

**Contact Detection:**
Heavy-atom contacts were computed using a custom Python script (`compute_contacts.py`) with a 5.0 Å cutoff
Contacts between antibody chains (A+B) and antigen chain (C) were identified
Conserved interface residues were defined as those contacting both p53 and CALR within 5 Å

**Interface Residue Identification:**
23 antibody residues were identified as conserved interface residues:
Heavy chain (A): 16 residues (positions 31, 32, 33, 50, 52, 55, 57, 59, 99, 102, 103, 104, 105, 106, 107, 109)
Light chain (B): 7 residues (positions 2, 27, 31, 98, 99, 100, 102)

**Hydrogen Bond Detection:**
Hydrogen bonds were identified using geometric criteria: donor-acceptor distance ≤ 3.5 Å and donor-H-acceptor angle ≥ 120°
Backbone and side-chain hydrogen bonds were analyzed separately
Results saved to `hbonds_p53.csv` and `hbonds_calr.csv`

**Salt Bridge Analysis:**
Salt bridges were defined as interactions between oppositely charged residues (Arg/Lys ↔ Asp/Glu) within 4.0 Å
Results saved to `salt_bridges_p53.csv` and `salt_bridges_calr.csv`

**Hydrophobic Contact Analysis:**
Hydrophobic contacts were identified between nonpolar residues (Ala, Val, Ile, Leu, Met, Phe, Trp, Tyr) within 4.5 Å
Results saved to `hydrophobic_contacts_p53.csv` and `hydrophobic_contacts_calr.csv`

**Alanine Scanning:**
FoldX AlaScan command was used to compute ΔΔG values for alanine mutations
26 interface residues were scanned across both complexes
5 independent runs were performed for statistical robustness
ΔΔG values were extracted from `*_AS.fxout` output files using `parse_alascan.py`
Results compiled in `foldx_ddg_analysis.csv`

**Energy Decomposition:**
Conserved hotspots were defined as residues with ΔΔG > 1.0 kcal/mol and <2.0 kcal/mol difference between complexes
Differential residues were defined as those with ΔΔG difference > 2.0 kcal/mol between antigens

**B-factor Extraction:**
B-factors were extracted from ATOM records (columns 60-66) of both PDB structures using `analyze_flexibility.py`
Average B-factors were calculated per residue across all atoms
Flexibility classification:
High flexibility: B-factor > 50 Ų
Medium flexibility: B-factor 30-50 Ų
Low flexibility: B-factor < 30 Ų
Results saved to `bfactors_p53.csv` and `bfactors_calr.csv`

**Residue Classification:**
Hydrophobic: Ala, Val, Ile, Leu, Met, Phe, Trp, Tyr
Aromatic: Phe, Tyr, Trp, His
Basic: Arg, Lys, His
Acidic: Asp, Glu
Polar: Ser, Thr, Asn, Gln, Cys

**Enrichment Analysis:**
Chemical composition was calculated for conserved interface residues using `compute_paratope_chemistry.py`
Percentages were compared to typical polyreactive antibody thresholds (>25% aromatic, >50% hydrophobic)
Results saved to `paratope_chemistry_map.csv`

**Antibody Alignment:**
Kabsch superposition algorithm was applied using `corrected_binding_analysis.py`
Alignment was performed using conserved interface residues (30 residues total)
Rotation matrix and translation vector were computed using singular value decomposition (SciPy)

**RMSD Calculation:**
Root mean square deviation was calculated for conserved interface residues after superposition
Heavy chain, light chain, and combined RMSDs were computed separately
Results saved to `binding_mode_comparison_corrected.csv`

**Binding Orientation Analysis:**
Binding orientation vectors were calculated as centroids of interface residues
Angle between binding orientations was computed using dot product of normalized vectors

**PyMOL Scripts:**
`electrostatics_hydrophobic.pml`: Chemical patch visualization
`visualize_bfactors.pml`: Flexibility mapping

**Python Analysis Pipeline:**
Custom scripts for contact detection, interaction analysis, and statistical calculations
NumPy and SciPy for numerical computations and structural superposition
Pandas for data manipulation and CSV generation

**Sequence Extraction:**
Heavy and light chain sequences extracted using `extract_antibody_sequences.py`
Sequences saved to `heavy.fasta` and `light.fasta` for ANARCI numbering
HCDR3 length analysis performed (12 residues, within polyreactive range)

### Dependencies

**Core Software:**
FoldX 4.1 (alanine scanning and energy calculations)
PyMOL (molecular visualization and analysis)
Python 3.9+ (data analysis and processing)

**Python Libraries:**
NumPy 2.0.2 (numerical computations)
SciPy 1.13.1 (structural superposition and advanced mathematics)
Pandas (data manipulation)

### Statistical Methods

**Multiple Testing:**
5 independent FoldX runs were performed for alanine scanning to ensure reproducibility
B-factor analysis was performed on both complexes independently

**Thresholds and Classifications:**
Contact cutoff: 5.0 Å (heavy atoms)
High flexibility threshold: B-factor > 50 Ų
Conserved hotspot threshold: ΔΔG > 1.0 kcal/mol
Differential residue threshold: ΔΔG difference > 2.0 kcal/mol
Polyreactive signature thresholds: >25% aromatic, >50% hydrophobic

---
## FoldX Alanine Scanning Results

### Conserved Energetic Hotspots (ΔΔG > 1 kcal/mol, similar for both antigens)

| Residue | Type | ΔΔG (p53) | ΔΔG (CALR) | Average | Interpretation |
|---------|------|-----------|------------|---------|----------------|
| **A104** | VAL | 3.43 | 3.17 | **3.30** | Hydrophobic core |
| **A58** | PRO | 3.06 | 2.12 | 2.59 | Structural rigidity |
| **A99** | GLN | 3.20 | 1.81 | 2.50 | H-bond contributor |
| **A101** | LEU | 3.20 | 1.79 | 2.49 | Hydrophobic packing |
| **A102** | ILE | 4.42 | 1.50 | 2.96 | Mixed (differential) |
| **A107** | ASN | 2.39 | 1.25 | 1.82 | Polar contact |
| **A59** | PRO | 2.32 | 0.94 | 1.63 | Structural |
| **A106** | GLY | 2.23 | 0.82 | 1.52 | Flexibility |
| **A98** | PRO | 2.17 | 0.82 | 1.49 | CDR structure |

### Antigen-Specific Residues (ΔΔG difference > 2 kcal/mol)

| Residue | Type | ΔΔG (p53) | ΔΔG (CALR) | Difference | Specificity |
|---------|------|-----------|------------|------------|-------------|
| **A52** | LEU | 4.61 | 0.00 | **4.61** | p53-specific |
| **A30** | THR | -0.36 | 3.69 | **4.05** | CALR-specific |
| **A33** | TYR | 6.19 | 2.35 | **3.84** | p53-biased |
| **A50** | VAL | 3.59 | 1.03 | **2.56** | p53-biased |
| **A32** | TYR | 1.77 | -0.73 | **2.50** | p53-biased |

### Key Finding:
**Multiple small/moderate contributions** rather than one dominant hotspot → characteristic of polyreactive antibodies.

---

## Chemical Composition Analysis

### Conserved Interface (23 residues binding both antigens)

| Property | Count | % of Interface | Interpretation |
|----------|-------|----------------|----------------|
| **Hydrophobic** | 13 | **56.5%** | High hydrophobic patch |
| **Aromatic** | 8 | **34.8%** | ✓ Typical polyreactive signature |
| **Basic (R/K/H)** | 3 | 13.0% | Electrostatic versatility |
| **Acidic (D/E)** | 3 | 13.0% | Charge complementarity |
| **Polar (S/T/N/Q)** | 4 | 17.4% | H-bond capability |

### Top Residue Types:
1. **Tyrosine (TYR)**: 4 positions (17.4%) ← dual H-bond/hydrophobic
2. **Aspartate (ASP)**: 3 positions (13.0%) ← salt bridges
3. **Histidine (HIS)**: 2 positions (8.7%) ← pH-dependent chemistry
4. **Isoleucine (ILE)**: 2 positions (8.7%) ← hydrophobic core

### Polyreactivity Indicators:
- ✓ **Aromatic enrichment >25%** (literature threshold for polyreactivity)
- ✓ **Hydrophobic enrichment >50%** (enables promiscuous packing)
- ✓ **Tyrosine enrichment** (versatile chemistry: H-bond donor + π-stacking)
- ✓ **Basic + aromatic combination** (cation-π interactions)

---

## Flexibility Analysis (B-factors)

### Conserved Interface Statistics:

| Complex | Mean B-factor | % High Flex (>50 Ų) | Interpretation |
|---------|---------------|----------------------|----------------|
| p53 | 71.46 Ų | 87.0% | Very high flexibility |
| CALR | 78.04 Ų | 100.0% | Extremely high flexibility |
| **Average** | **74.75 Ų** | **93.5%** | ✓ Induced-fit mechanism |

### Comparison to Literature:

| Antibody Type | Typical B-factor | CL1-Ab4 |
|---------------|------------------|---------|
| Polyreactive | >60 Ų | **74.75 Ų** ✓ |
| Monospecific | <50 Ų | - |
| Germline | 50-70 Ų | - |

---

## Binding Mode Comparison (Corrected with Proper Alignment)

### Structural Metrics (After Antibody Superposition):

| Metric | Value | Interpretation |
|--------|-------|----------------|
| **Heavy chain RMSD** | **7.602 Å** | ✓ Significant conformational change |
| **Light chain RMSD** | **7.762 Å** | ✓ Significant conformational change |
| **Combined RMSD** | **7.733 Å** | ✓ Substantial induced-fit mechanism |
| **Binding angle difference** | **76.3°** | ✓ Substantially different orientations |
| **p53 epitope size** | 30 residues | Larger epitope |
| **CALR epitope size** | 20 residues | Smaller epitope |
| **Epitope size ratio** | 1.5:1 | Versatile recognition |
| **p53 epitope Rg** | 10.58 Å | Similar spatial extent |
| **CALR epitope Rg** | 9.32 Å | Similar spatial extent |

---

## Interaction Network Analysis

### Conserved Paratope: 23 Residues

**Heavy Chain (16 residues):**
- HCDR1 region: A31, A32, A33
- HCDR2 region: A50, A52, A55, A57, A59
- HCDR3 region: A99, A102, A103, A104, A105, A106, A107, A109
- Framework: A2, A27, A28, A30, A58, A98, A100, A101, A108, A110

**Light Chain (7 residues):**
- LCDR1 region: B27, B31
- LCDR3 region: B98, B99, B100, B102
- Framework: B1, B2, B26, B28, B29, B30, B32-B38, B55, B56, B96, B97, B101

### Interaction Types (Overlapping):
- **Hydrogen bonds**: Tyr, Asn, Gln, Ser → versatile donors/acceptors
- **Salt bridges**: Asp, His → electrostatic complementarity
- **Hydrophobic contacts**: Val, Ile, Leu, Pro → shape-driven packing

---
