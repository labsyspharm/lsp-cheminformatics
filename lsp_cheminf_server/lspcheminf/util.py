from rdkit import Chem
from rdkit.Chem import inchi

identifier_mol_mapping = {
    "smiles": Chem.MolFromSmiles,
    "inchi": inchi.MolFromInchi,
    "smarts": Chem.MolFromSmarts,
}

mol_identifier_mapping = {
    "smiles": Chem.MolToSmiles,
    "inchi": inchi.MolToInchi,
    "inchi_key": inchi.MolToInchiKey,
}


def convert_compound_request(compounds):
    names = compounds.get("names", list(range(len(compounds["compounds"]))))
    return {
        n: identifier_mol_mapping[compounds["identifier"]](m)
        for n, m in zip(names, compounds["compounds"])
    }
