# WNV_AntiviralDrug_Repurposing

Molecular docking–based analysis to identify potential drug candidates targeting the NS5 protein of West Nile Virus (WNV).

## Overview

This project demonstrates a computational workflow for identifying potential antiviral compounds targeting the NS5 protein of the West Nile Virus (WNV) TX strain. 
The workflow integrates protein structure analysis, ligand preparation, molecular docking, and binding affinity ranking, showcasing reproducible and analytical skills applicable in computational biology and cheminformatics pipelines.

## Key Features
* Structured Workflow: End-to-end pipeline from protein and ligand preparation to docking and ranking.
* Data Handling & Transformation: Efficient processing of protein structures and compound libraries.
* Analytical Decision-Making: Binding affinity predictions guide selection of top candidate compounds.
* Reproducibility: Clear steps and scripts allow full replication of analysis.

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
