import os
import subprocess
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation

DATA_DIR = "data"
RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

# --- 1. Prepare Ligand ---
print("Preparing Pyrimethamine ligand...")
pyrimethamine_smiles = "C1=CC(=CC=C1C2=NC(=NC(=N2)N)N)Cl"
mol = Chem.MolFromSmiles(pyrimethamine_smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, randomSeed=42)
preparator = MoleculePreparation()
preparator.prepare(mol)
ligand_pdbqt = os.path.join(DATA_DIR, "pyrimethamine.pdbqt")
preparator.write_pdbqt_file(ligand_pdbqt)
print("Ligand ready.")

# --- 2. Clean Receptors & Find Center ---
def clean_and_center(filepath):
    if not os.path.exists(filepath):
        return None
    with open(filepath, 'r') as f:
        lines = f.readlines()
    clean_lines = [line for line in lines if line.startswith(('ATOM', 'HETATM'))]
    with open(filepath, 'w') as f:
        f.writelines(clean_lines)
    
    coords = []
    for line in clean_lines:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        coords.append([x, y, z])
    
    if not coords:
        return None
    coords = np.array(coords)
    center_x = float(np.mean(coords[:, 0]))
    center_y = float(np.mean(coords[:, 1]))
    center_z = float(np.mean(coords[:, 2]))
    return center_x, center_y, center_z

# We calculate center based on the whole protein (Wild-Type)
center_coords = clean_and_center(os.path.join(DATA_DIR, "wild_type_receptor.pdbqt"))
if center_coords:
    center_x, center_y, center_z = center_coords
    print(f"Calculated protein center: X={center_x:.2f}, Y={center_y:.2f}, Z={center_z:.2f}")
else:
    center_x, center_y, center_z = 0, 0, 0

# --- 3. Run Docking (Blind Docking) ---
# Increase box size to cover the whole protein
size_x, size_y, size_z = 60, 60, 60 

for struct in ["wild_type", "mutant"]:
    print(f"\n--- Docking {struct} (Blind Docking) ---")
    receptor = os.path.join(DATA_DIR, f"{struct}_receptor.pdbqt") 
    out_file = os.path.join(RESULTS_DIR, f"{struct}_docked.pdbqt")
    log_file = os.path.join(RESULTS_DIR, f"{struct}_log.txt")

    command = [
        "vina.exe",
        "--receptor", receptor,
        "--ligand", ligand_pdbqt,
        "--center_x", str(center_x),
        "--center_y", str(center_y),
        "--center_z", str(center_z),
        "--size_x", str(size_x),
        "--size_y", str(size_y),
        "--size_z", str(size_z),
        "--out", out_file,
        "--exhaustiveness", "8"
    ]

    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        with open(log_file, "w", encoding="utf-8") as f:
            f.write(result.stdout)
            f.write(result.stderr)
        
        lines = result.stderr.split('\n')
        for line in lines:
            if line.strip().startswith("1 "):
                print(f">>> BEST AFFINITY FOR {struct.upper()}: {line.split()[1]} kcal/mol <<<")
                break
    except subprocess.CalledProcessError as e:
        print(f"Vina error for {struct}: {e.stderr}")

print("\nAll done!")