from typing import Mapping, List, Tuple, Any

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
    "smarts": Chem.MolToSmarts,
}


def convert_compound_request(
    compounds: Mapping
) -> Tuple[Mapping[str, Chem.Mol], List[Tuple[str, str]]]:
    names = map(str, compounds.get("names", list(range(len(compounds["compounds"])))))
    mapper = identifier_mol_mapping[compounds["identifier"]]
    skipped = []
    converted = {}
    for n, m in zip(names, compounds["compounds"]):
        try:
            mol = mapper(m)
        except Exception:
            skipped.append((n, m))
            continue
        if mol is not None:
            converted[n] = mol
        else:
            skipped.append((n, m))
    return converted, skipped
