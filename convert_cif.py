import gemmi
import os

DATA_DIR = "data"

print("Converting CIF to PDB using gemmi...")

# Convert Wild-Type
st_wt = gemmi.read_structure(os.path.join(DATA_DIR, "wild_type.cif"))
st_wt.setup_entities()
st_wt.write_pdb(os.path.join(DATA_DIR, "wild_type.pdb"))
print("Saved wild_type.pdb")

# Convert Mutant
st_mut = gemmi.read_structure(os.path.join(DATA_DIR, "mutant.cif"))
st_mut.setup_entities()
st_mut.write_pdb(os.path.join(DATA_DIR, "mutant.pdb"))
print("Saved mutant.pdb")

print("Done!")