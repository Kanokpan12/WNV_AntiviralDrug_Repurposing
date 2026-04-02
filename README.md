# WNV_AntiviralDrug_Repurposing

Molecular docking–based analysis to identify potential drug candidates targeting the NS5 protein of West Nile Virus (WNV).

## Overview

This project demonstrates a computational workflow for identifying potential antiviral compounds targeting the NS5 protein of the West Nile Virus (WNV) Texas strain. The workflow integrates AI-based protein structure prediction, structural validation, druggable pocket identification, and molecular docking to support antiviral drug repurposing efforts.

The AlphaFold2-predicted NS5 structure showed high confidence (mean pLDDT score = 94.5) and strong agreement with experimentally determined structures. Structural alignment confirmed conservation of key catalytic domains:

MTase domain aligned with PDB: 2OY0 (RMSD = 0.489 Å over 255 Cα atoms)
RdRp domain aligned with PDB: 2HCN (RMSD = 1.044 Å over 353 Cα atoms)

Druggable pocket analysis identified a high-confidence binding region within the MTase domain (residues 80–95) with a druggability score of 0.8.

Molecular docking results identified **Tecovirimat** as a promising candidate compound with favorable predicted binding to NS5, suggesting potential as a repurposed antiviral lead against the neuroinvasive WNV Texas strain.

## Key Features
1. AI-driven structural modeling
* AlphaFold2 prediction of WNV NS5 protein structure
2. Structure validation
* RMSD comparison against experimentally resolved MTase and RdRp domains
3. Binding pocket identification
* Drug-target site detection using pocket scoring tools
4. Docking-based compound prioritization
* Screening antiviral compounds from ChEMBL
5. Reproducible computational workflow
* End-to-end pipeline from structure prediction to candidate ranking

## Workflow Steps
1. WNV_TX Protein Structure:
   * Use AlphaFold2 to predict the WNV_TX NS5 protein structure.
   * Compare predicted structures with experimental structures (RMSD measurement) to validate the accuracy.
2. Binding Pocket Identification:
   * Identify druggable sites suitable for small molecule binding
   * Use DOGSiteScorer or similar tools for pocket scoring.
3. Ligand & Protein Preparation:
   * Convert protein and ligand files to docking-ready formats (.pdb → .pdbqt).
4. Molecular Docking:
   * Perform docking using Smina (AutoDock Vina fork) to predict binding conformations.
5. Binding Affinity Estimation & Ranking:
   * Predict binding free energies (K<sub>D</sub>) and rank compounds.
6. Candidate NS5 Inhibitors:
   * Select top compounds for further experimental or computational validation.
  
## Tools & data sources
* Smina (AutoDock Vina fork) for docking simulations
* ChEMBL database for antiviral compounds
* Protein/ligand preparation tools
* Python/R scripts for post-docking analysis and visualization
