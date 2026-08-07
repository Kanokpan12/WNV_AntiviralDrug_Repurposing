# West Nile virus_Antiviral_Drug_Repurposing

Molecular docking–based analysis to identify potential drug candidates targeting the NS5 protein of West Nile Virus (WNV).

## Overview

This project demonstrates a computational workflow for identifying potential antiviral compounds targeting the NS5 protein of the West Nile Virus (WNV) Texas strain. The workflow integrates AI-based protein structure prediction, structural validation, druggable pocket identification, and molecular docking to support antiviral drug repurposing efforts.

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
7. MD simulation
   Run molecular dynamics simulations on top candidates to assess binding stability.

## Requirements
* Python ≥ 3.9 (RDKit, pandas, numpy): data preparation and visualization
* R ≥ 4.0: visualization
* UCSF chimeraX: structure visualization
* Compound library: ChEMBL
* Binding affinity: AutoDock/Smina
* MD simulation: GROMACS
* Structure prediction: AlphaFold2 (or ColabFold)
* Pocket detection: DoGSiteScorer (via ProteinsPlus web server or local install)

## Result 
* The AlphaFold2-predicted NS5 structure showed high confidence (mean pLDDT score = 94.5) and strong agreement with experimentally determined structures. Structural alignment confirmed conservation of key catalytic domains:
 <img width="301" height="290" alt="image" src="https://github.com/user-attachments/assets/0723f68b-d756-4677-b517-2e54855b58df" />


* MTase domain aligned with PDB: 2OY0 (RMSD = 0.489 Å over 255 Cα atoms) 
* RdRp domain aligned with PDB: 2HCN (RMSD = 1.044 Å over 353 Cα atoms)

* Druggable pocket analysis identified a high-confidence binding region within the MTase domain (residues 80–95) with a druggability score of 0.8.

* Molecular docking results identified **Tecovirimat** as a promising candidate compound with favorable predicted binding to NS5, suggesting potential as a repurposed antiviral lead against the neuroinvasive WNV Texas strain.
<img width="288" height="286" alt="image" src="https://github.com/user-attachments/assets/5a3357d3-9df2-47c1-8396-ee43a82f27a8" />

* Protein RMSD. TEC binding improved the structural stability of NS5 by reducing the average backbone RMSD by approximately **20%** (0.51 → 0.40 Å). While the apo NS5 system did not fully converge during the 10 ns simulation, the TEC-bound complex reached a stable equilibrium, indicating that ligand binding actively stabilizes the protein structure.

* Residue Flexibility (RMSF).** TEC substantially decreased residue mobility, reducing the mean RMSF by more than **50%** (3.92 → 1.69 Å). The greatest reduction was observed in **residues 80–107**, suggesting that this binding-loop region becomes significantly more rigid upon ligand binding.

* Ligand RMSD. The bound TEC ligand remained highly stable throughout the simulation, exhibiting an average RMSD of only **0.09 Å** with virtually no displacement from its initial docking pose. The absence of significant drift over 10 ns indicates a persistent and well-maintained binding orientation.

* Hydrogen Bond Analysis. The TEC–NS5 complex consistently maintained **1–2 hydrogen bonds** throughout the simulation. A slight increase in hydrogen-bond fluctuations after approximately **6 ns** coincided with a minor RMSD change, suggesting a small conformational adjustment while preserving overall binding stability.
<img width="647" height="510" alt="image" src="https://github.com/user-attachments/assets/d423bb6d-4cc0-4931-9978-8f85fee3b30f" />

* Free Energy Landscape (FEL). The apo NS5 protein sampled a broad and rugged free energy landscape, indicating high conformational flexibility and dynamic structural fluctuations. In contrast, the TEC-bound NS5 complex collapsed into a single, well-defined energy minimum, demonstrating that TEC restricts the conformational space accessible to the protein. This reduction in conformational entropy suggests that TEC effectively locks NS5 into a stable, low-energy state.
<img width="609" height="231" alt="image" src="https://github.com/user-attachments/assets/aba80e89-2833-4c8d-8407-f323ca7954ae" />
<img width="601" height="251" alt="image" src="https://github.com/user-attachments/assets/4943ab45-404d-45da-8344-8156544d1e50" />






