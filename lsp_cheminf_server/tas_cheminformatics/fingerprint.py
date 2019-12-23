from functools import partial

from rdkit import DataStructs
from rdkit.Chem import AllChem

fingerprint_functions = {
    "topological": AllChem.RDKFingerprint,
    "morgan": partial(AllChem.GetMorganFingerprintAsBitVect, radius=2),
}


def calculate_fingerprint(mol, fingerprint_type="morgan", fingerprint_args={}):
    if not fingerprint_type in fingerprint_functions:
        raise ValueError(
            "`fingerprint_type` must be one of ", list(fingerprint_functions.keys())
        )
    return fingerprint_functions[fingerprint_type](mol, **fingerprint_args)


def calculate_similarity(
    query, targets, fingerprint_type="morgan", fingerprint_args={}
):
    query_fp = calculate_fingerprint(query, fingerprint_type, fingerprint_args)
    return {
        n: DataStructs.FingerprintSimilarity(
            query_fp, calculate_fingerprint(t, fingerprint_type, fingerprint_args)
        )
        for n, t in targets.items()
    }
