from functools import partial
from typing import Mapping, List, Any, Tuple

from rdkit.Chem import Mol
from rdkit import DataStructs
from rdkit.Chem import AllChem

fingerprint_functions = {
    "topological": AllChem.RDKFingerprint,
    "morgan": partial(AllChem.GetMorganFingerprintAsBitVect, radius=2),
}


def calculate_fingerprint(
    mol: Mol, fingerprint_type: str = "morgan", fingerprint_args: Mapping[str, Any] = {}
) -> DataStructs.ExplicitBitVect:
    if not fingerprint_type in fingerprint_functions:
        raise ValueError(
            "`fingerprint_type` must be one of ", list(fingerprint_functions.keys())
        )
    return fingerprint_functions[fingerprint_type](mol, **fingerprint_args)


def calculate_similarity(
    query: Mol,
    targets: Mapping[str, Mol],
    fingerprint_type: str = "morgan",
    fingerprint_args: Mapping[str, Any] = {},
) -> Mapping[str, float]:
    query_fp = calculate_fingerprint(query, fingerprint_type, fingerprint_args)
    return {
        n: DataStructs.FingerprintSimilarity(
            query_fp, calculate_fingerprint(t, fingerprint_type, fingerprint_args)
        )
        for n, t in targets.items()
    }


def find_substructure_matches(
    query: Mol, targets: Mapping[str, Mol], substructure_args: Mapping[str, Any] = {}
) -> Mapping[str, Tuple[Tuple[int]]]:
    res = {}
    for n, t in targets.items():
        match = t.GetSubstructMatches(query, **substructure_args)
        if len(match) > 0:
            res[n] = match
    return res
