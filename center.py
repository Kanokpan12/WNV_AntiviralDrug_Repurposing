from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

# Load your pocket PDB
pdb_file = "ns5dq176637clean_P_1_res.pdb"
mol = Chem.MolFromPDBFile(pdb_file, removeHs=False, sanitize=False)

# Get coordinates of all atoms
conf = mol.GetConformer()
coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])

# Center of pocket
center = coords.mean(axis=0)
center_x, center_y, center_z = center
print(f"Center coordinates: {center_x:.2f}, {center_y:.2f}, {center_z:.2f}")

# Pocket dimensions (max-min) + padding
padding = 2.0  # Å
size = coords.max(axis=0) - coords.min(axis=0) + padding
size_x, size_y, size_z = size
print(f"Box size: {size_x:.2f}, {size_y:.2f}, {size_z:.2f}")
