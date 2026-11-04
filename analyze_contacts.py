import os
from Bio.PDB import PDBParser, NeighborSearch

DATA_DIR = "data"
RESULTS_DIR = "results"

def get_best_pose_ligand_coords(pdbqt_file):
    """Extracts coordinates of the first (best) pose from a Vina output PDBQT."""
    coords = []
    with open(pdbqt_file, 'r') as f:
        for line in f:
            # The best pose is usually the first MODEL or just the first block of ATOM/HETATM
            if line.startswith("ENDMDL"):
                break # Stop after the first model
            if line.startswith(("ATOM", "HETATM")):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append((x, y, z))
    return coords

def find_contacts(protein_pdb, ligand_coords, cutoff=4.0):
    """Finds protein residues within a cutoff distance of the ligand."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", protein_pdb)
    
    # Get all atoms from the protein
    atoms = [atom for atom in structure.get_atoms()]
    
    # Create a NeighborSearch object
    ns = NeighborSearch(atoms)
    
    # Find all protein atoms within the cutoff of any ligand atom
    nearby_residues = set()
    for lx, ly, lz in ligand_coords:
        # Find neighbors within cutoff distance
        neighbors = ns.search((lx, ly, lz), cutoff)
        for atom in neighbors:
            residue = atom.get_parent()
            # Get residue name and number (e.g., ALA 45)
            resname = residue.get_resname()
            resnum = residue.get_id()[1]
            nearby_residues.add(f"{resname}{resnum}")
            
    return sorted(list(nearby_residues))

# Main Analysis
print("Analyzing binding site interactions (cutoff 4.0 Å)...\n")

# Analyze Wild-Type
wt_ligand = os.path.join(RESULTS_DIR, "wild_type_docked.pdbqt")
wt_protein = os.path.join(DATA_DIR, "wild_type.pdb")
if os.path.exists(wt_ligand) and os.path.exists(wt_protein):
    coords = get_best_pose_ligand_coords(wt_ligand)
    contacts = find_contacts(wt_protein, coords)
    print(f"Wild-Type Contacts ({len(contacts)} residues):")
    print(", ".join(contacts))
else:
    print("Wild-Type files not found.")

print("\n" + "-"*40 + "\n")

# Analyze Mutant
mut_ligand = os.path.join(RESULTS_DIR, "mutant_docked.pdbqt")
mut_protein = os.path.join(DATA_DIR, "mutant.pdb")
if os.path.exists(mut_ligand) and os.path.exists(mut_protein):
    coords = get_best_pose_ligand_coords(mut_ligand)
    contacts = find_contacts(mut_protein, coords)
    print(f"Mutant Contacts ({len(contacts)} residues):")
    print(", ".join(contacts))
else:
    print("Mutant files not found.")