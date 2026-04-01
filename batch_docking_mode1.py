#!/usr/bin/env python3

import os
import glob
import subprocess
import csv

# === USER SETTINGS ===
receptor = r"/Volumes/Seagate Por/DrugDesign/resultns5pocket/ns5dq176637clean_P_1_res.pdbqt"
ligand_folder = r"/Volumes/Seagate Por/DrugDesign/pdbqt_output"   # 32 ligands here
output_folder = r"/Volumes/Seagate Por/DrugDesign/docked_result"
smina_bin = r"/Volumes/Seagate Por/DrugDesign/smina/smina.osx.12"

# Grid box
center_x, center_y, center_z = 11.44, 13.70, -19.36
size_x, size_y, size_z = 17.85, 17.97, 25.48
exhaustiveness = 8
# ====================

# Create output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Collect results
results = []

# List ligand files
ligand_files = glob.glob(os.path.join(ligand_folder, "*.pdbqt"))
if not ligand_files:
    raise SystemExit(f"❌ No ligand pdbqt files found in {ligand_folder}")

for ligand_file in ligand_files:
    ligand_name = os.path.basename(ligand_file).replace(".pdbqt", "")
    out_file = os.path.join(output_folder, ligand_name + "_docked.pdbqt")

    print(f"\nDocking {ligand_name} ...")

    # Run Smina docking
    cmd = [
        smina_bin,
        "-r", receptor,
        "-l", ligand_file,
        "--center_x", str(center_x),
        "--center_y", str(center_y),
        "--center_z", str(center_z),
        "--size_x", str(size_x),
        "--size_y", str(size_y),
        "--size_z", str(size_z),
        "--exhaustiveness", str(exhaustiveness),
        "--out", out_file
    ]
    subprocess.run(cmd)

    # Parse mode 1 REMARK VINA RESULT
    mode1_affinity = None
    rmsd_lb = None
    rmsd_ub = None

    with open(out_file, "r", errors="ignore") as f:
        for line in f:
            if "REMARK VINA RESULT:" in line:
                parts = line.strip().split()
                if len(parts) >= 4:
                    # First REMARK VINA RESULT = mode 1
                    mode1_affinity = float(parts[3])
                    rmsd_lb = float(parts[4])
                    rmsd_ub = float(parts[5])
                    break  # only mode 1

    results.append([ligand_name, mode1_affinity, rmsd_lb, rmsd_ub, out_file])
    print(f"✅ {ligand_name}: affinity={mode1_affinity} kcal/mol, RMSD_lb={rmsd_lb}, RMSD_ub={rmsd_ub}")

# Write CSV summary
csv_file = os.path.join(output_folder, "docking_mode1_summary.csv")
with open(csv_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Ligand", "Mode1_Affinity_kcal/mol", "RMSD_lb", "RMSD_ub", "Docked_File"])
    writer.writerows(results)

print("\n🎉 All dockings finished!")
print("CSV saved at:", csv_file)
