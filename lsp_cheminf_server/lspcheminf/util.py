import warnings
from typing import Mapping, List, Tuple, Any, Union

from rdkit import Chem
from rdkit.Chem import inchi

try:
    import chemfp

    chemfp_available = True
    import chemfp.search
    import chemfp.arena
    import chemfp.rdkit_types
except ImportError:
    chemfp_available = False

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

Molmap = Mapping[str, Chem.Mol]
MolmapArena = Union[Molmap, chemfp.arena.FingerprintArena]


def convert_compound_request(
    compounds: Mapping
) -> Tuple[Molmap, List[Tuple[str, str]]]:
    names = map(str, compounds.get("names", list(range(len(compounds["compounds"])))))
    mapper = identifier_mol_mapping[compounds["identifier"]]
    skipped = []
    converted = {}
    for n, m in zip(names, compounds["compounds"]):
        try:
            with warnings.catch_warnings():
                mol = mapper(m)
        except Exception:
            skipped.append((n, m))
            continue
        if mol is not None:
            converted[n] = mol
        else:
            skipped.append((n, m))
    return converted, skipped
