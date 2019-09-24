from rdkit import Chem
from rdkit.Chem import inchi

identifier_mol_mapping = {"smiles": Chem.MolFromSmiles, "inchi": inchi.MolFromInchi}

mol_identifier_mapping = {
    "smiles": Chem.MolToSmiles,
    "inchi": inchi.MolToInchi,
    "inchi_key": inchi.MolToInchiKey,
}
