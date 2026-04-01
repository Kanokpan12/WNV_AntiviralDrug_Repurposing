#!/usr/bin/env python3
"""
smiles_to_pdbqt_full.py
Full pipeline:
 - read CSV
 - sanitize names
 - convert SMILES -> 3D RDKit mol -> PDB
 - convert PDB -> PDBQT using OpenBabel (obabel)
 - log results to CSV
 - zip output folder
"""

import os
import re
import csv
import shutil
import subprocess
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# ----------------------
# SETTINGS - EDIT AS NEEDED
# ----------------------
csv_file = "/Volumes/Seagate Por/DrugDesign/Antiviral2_Simplified.csv"
smiles_column = "Smiles"
name_column = "Name"
fallback_id_column = "ChEMBL_ID"   # fallback if Name missing/empty
output_folder = "/Volumes/Seagate Por/DrugDesign/pdbqt_output"
log_csv = os.path.join(output_folder, "conversion_log.csv")
zip_output = "/Volumes/Seagate Por/DrugDesign/pdbqt_output.zip"
# ----------------------

os.makedirs(output_folder, exist_ok=True)

def safe_filename(s, maxlen=120):
    """Make a filesystem-safe filename from string s."""
    if not isinstance(s, str) or s.strip()=="":
        return ""
    s = s.strip()
    # replace spaces with underscore
    s = s.replace(" ", "_")
    # remove characters not alphanumeric, underscore, hyphen, dot
    s = re.sub(r"[^A-Za-z0-9_\-\.]", "", s)
    # shorten
    if len(s) > maxlen:
        s = s[:maxlen]
    return s

def unique_name(base, existing):
    """Return unique filename based on base not in existing set."""
    candidate = base
    i = 1
    while candidate in existing:
        candidate = f"{base}_{i}"
        i += 1
    existing.add(candidate)
    return candidate

# Read CSV into dataframe
df = pd.read_csv(csv_file, dtype=str)  # read as strings
df = df.fillna("")  # replace NaN with empty strings

# Prepare log rows
log_rows = []
existing_names = set()

for idx, row in df.iterrows():
    original_name = row.get(name_column, "").strip()
    fallback_id = row.get(fallback_id_column, "").strip()
    smi = row.get(smiles_column, "").strip()

    # Build a base name: prefer Name, else fallback id, else index
    base = original_name or fallback_id or f"mol_{idx}"
    base_safe = safe_filename(base)
    if base_safe == "":
        base_safe = f"mol_{idx}"
    file_base = unique_name(base_safe, existing_names)

    log_entry = {
        "row_index": idx,
        "requested_name": original_name,
        "used_name": file_base,
        "ChEMBL_ID": fallback_id,
        "smiles": smi,
        "status": "",
        "message": ""
    }

    if smi == "":
        log_entry["status"] = "FAILED"
        log_entry["message"] = "Empty SMILES"
        log_rows.append(log_entry)
        print(f"[{idx}] {file_base}: Empty SMILES — skipped")
        continue

    # parse SMILES
    mol = None
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise ValueError("RDKit returned None for MolFromSmiles")
    except Exception as e:
        log_entry["status"] = "FAILED"
        log_entry["message"] = f"Bad SMILES / RDKit failed: {e}"
        log_rows.append(log_entry)
        print(f"[{idx}] {file_base}: Bad SMILES — {e}")
        continue

    try:
        mol = Chem.AddHs(mol)
        # Embed - try a few times if embedding fails
        emb_ok = False
        for attempt in range(3):
            res = AllChem.EmbedMolecule(mol, randomSeed=42+attempt)
            if res == 0:
                emb_ok = True
                break
        if not emb_ok:
            raise RuntimeError("Embedding failed after 3 attempts")
        # Optimize geometry
        uff_ok = AllChem.UFFOptimizeMolecule(mol)
        # even if optimization returns non-zero, proceed but log warning
    except Exception as e:
        log_entry["status"] = "FAILED"
        log_entry["message"] = f"3D generation failed: {e}"
        log_rows.append(log_entry)
        print(f"[{idx}] {file_base}: 3D generation failed — {e}")
        continue

    # Save PDB
    pdb_path = os.path.join(output_folder, f"{file_base}.pdb")
    try:
        Chem.MolToPDBFile(mol, pdb_path)
    except Exception as e:
        log_entry["status"] = "FAILED"
        log_entry["message"] = f"Failed to write PDB: {e}"
        log_rows.append(log_entry)
        print(f"[{idx}] {file_base}: Failed to write PDB — {e}")
        continue

    # Convert to PDBQT using obabel (OpenBabel)
    pdbqt_path = os.path.join(output_folder, f"{file_base}.pdbqt")
    try:
        # We already have 3D coords, so no --gen3d necessary; include Gasteiger charges
        cmd = ["obabel", "-ipdb", pdb_path, "-opdbqt", "-O", pdbqt_path, "--partialcharge", "gasteiger"]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            # capture stderr/stdout for debugging
            raise RuntimeError(f"obabel failed (rc={proc.returncode}): {proc.stderr.strip() or proc.stdout.strip()}")
    except Exception as e:
        log_entry["status"] = "FAILED"
        log_entry["message"] = f"OpenBabel failed: {e}"
        log_rows.append(log_entry)
        print(f"[{idx}] {file_base}: OpenBabel conversion failed — {e}")
        continue

    # Success
    log_entry["status"] = "SUCCESS"
    log_entry["message"] = "Converted to PDBQT"
    log_rows.append(log_entry)
    print(f"[{idx}] {file_base}: SUCCESS -> {pdbqt_path}")

# Save log CSV
with open(log_csv, "w", newline="", encoding="utf-8") as fh:
    fieldnames = ["row_index", "requested_name", "used_name", "ChEMBL_ID", "smiles", "status", "message"]
    writer = csv.DictWriter(fh, fieldnames=fieldnames)
    writer.writeheader()
    for r in log_rows:
        writer.writerow(r)

# Zip all pdbqt files
try:
    # shutil.make_archive takes a base name without extension
    base_zip_name = os.path.splitext(zip_output)[0]
    # remove existing zip if exists
    if os.path.exists(zip_output):
        os.remove(zip_output)
    shutil.make_archive(base_zip_name, 'zip', output_folder)
    print(f"Zipped output to {zip_output}")
except Exception as e:
    print(f"Failed to make zip: {e}")

print("Pipeline finished. See log:", log_csv)
