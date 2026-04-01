import subprocess
import glob
import os
import pandas as pd

# === Paths ===
smina_path = "/Volumes/Seagate Por/DrugDesign/smina/smina.osx.12"
receptor_path = "/Volumes/Seagate Por/DrugDesign/resultns5pocket/ns5dq176637clean_P_1_res.pdb"
ligand_folder = "/Volumes/Seagate Por/DrugDesign/pdbqt_output"
output_folder = "/Volumes/Seagate Por/DrugDesign/pdbqt_output/docked_results"
os.makedirs(output_folder, exist_ok=True)

# === Docking box parameters ===
center_x, center_y, center_z = 11.44, 13.70, -19.36
size_x, size_y, size_z = 17.85, 17.97, 25.48
exhaustiveness = 8

# === Collect docking results ===
results = []

# === Loop over ligands ===
for ligand_file in glob.glob(os.path.join(ligand_folder, "*.pdb")):
    ligand_name = os.path.basename(ligand_file).replace(".pdb","")
    
    # Prepare ligand path
    prepared_ligand = os.path.join(ligand_folder, f"{ligand_name}_prepared.pdb")
    
    # Prepare 3D + hydrogens using OpenBabel
    subprocess.run([
        "obabel", ligand_file, "-O", prepared_ligand, "--gen3d", "-h"
    ], check=True)
    
    # Output docked file
    docked_file = os.path.join(output_folder, f"{ligand_name}_docked.pdbqt")
    
    # Run Smina docking
    subprocess.run([
        smina_path,
        "-r", receptor_path,
        "-l", prepared_ligand,
        "--center_x", str(center_x),
        "--center_y", str(center_y),
        "--center_z", str(center_z),
        "--size_x", str(size_x),
        "--size_y", str(size_y),
        "--size_z", str(size_z),
        "--exhaustiveness", str(exhaustiveness),
        "--out", docked_file
    ], check=True)
    
    # Parse best binding affinity from docked file
    affinity = None
    with open(docked_file, "r") as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT"):
                affinity = float(line.split()[3])
                break
    if affinity is not None:
        results.append((ligand_name, affinity))

# === Save ranked results ===
df = pd.DataFrame(results, columns=["Ligand", "Affinity_kcal/mol"])
df.sort_values("Affinity_kcal/mol", inplace=True)  # lowest = best
df.to_csv(os.path.join(output_folder, "docking_results_ranked.csv"), index=False)

print("Docking finished! Ranked results saved to docking_results_ranked.csv")
