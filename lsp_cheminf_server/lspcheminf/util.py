import warnings
from typing import Mapping, List, Tuple, Any, Union, Callable

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

identifier_mol_mapping: Mapping[str, Callable] = {
    "smiles": Chem.MolFromSmiles,
    "inchi": inchi.MolFromInchi,
    "smarts": Chem.MolFromSmarts,
}

mol_identifier_mapping: Mapping[str, Callable] = {
    "smiles": Chem.MolToSmiles,
    "inchi": inchi.MolToInchi,
    "inchi_key": inchi.MolToInchiKey,
    "smarts": Chem.MolToSmarts,
}

Molmap = Mapping[str, Chem.Mol]
MolmapArena = Union[Molmap, chemfp.arena.FingerprintArena]


def convert_compound_request(
    compounds: Mapping, field: str = "compounds"
) -> Tuple[Union[Molmap, Mapping[str, str]], List[str]]:
    names = map(str, compounds.get("names", list(range(len(compounds["compounds"])))))
    mapper = identifier_mol_mapping[compounds["identifier"]]
    skipped = []
    converted = {}
    if field == "fingerprints":
        return {n: fp for n, fp in zip(names, compounds["fingerprints"])}, []
    for n, m in zip(names, compounds["compounds"]):
        try:
            with warnings.catch_warnings():
                mol = mapper(m)
        except Exception:
            skipped.append(n)
            continue
        if mol is not None:
            converted[n] = mol
        else:
            skipped.append(n)
    return converted, skipped
