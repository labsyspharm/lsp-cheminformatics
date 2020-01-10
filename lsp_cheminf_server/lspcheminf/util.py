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
    "smarts": Chem.MolToSmarts
}


def convert_compound_request(
    compounds: Mapping
) -> Tuple[Mapping[Any, Chem.Mol], List[Tuple[Any, str]]]:
    names = compounds.get("names", list(range(len(compounds["compounds"]))))
    skipped = []
    converted = {}
    for n, m in zip(names, compounds["compounds"]):
        try:
            mol = identifier_mol_mapping[compounds["identifier"]](m)
        except Exception:
            skipped.append((n, m))
            continue
        converted[n] = mol
    return converted, skipped
