from openbabel import openbabel as ob
from openbabel import pybel as pb

from rdkit import Chem
from rdkit.Chem import AllChem


def generate_conformers_rdkit(xyz_file, num_conformers=10, seed=42):
    """
    Generates conformers for a molecule using pybel to read XYZ and RDKit for conformer generation.

    Args:
        xyz_file: Path to the XYZ file.
        num_conformers: Number of conformers to generate.
        seed: Random seed for conformer generation.
    """

    # Read the molecule from the XYZ file using pybel
    mol = next(pb.readfile("xyz", xyz_file))
    mol.make3D()
    
    # Convert the pybel molecule to a mol block string
    initial_mol = mol.write("mol")

    # Create an RDKit Mol object from the mol block string
    mol = Chem.MolFromMolBlock(initial_mol, removeHs=False)
    if mol is None:
        raise ValueError("Could not create RDKit molecule from mol block")

    # Embed multiple conformers
    AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, randomSeed=seed)

    # Write conformers to separate XYZ files
    for i in range(mol.GetNumConformers()):
        Chem.MolToXYZFile(mol, f"conformer_{i+1}.xyz", confId=i)

    print ("{} conformers have been printed".format(len(mol.GetNumConformers())))


if __name__ == "__main__":
    # Example to generate conformers
    xyz_file = "molecule.xyz"  # Replace with your XYZ file path
    generate_conformers_rdkit(xyz_file, num_conformers=20)
