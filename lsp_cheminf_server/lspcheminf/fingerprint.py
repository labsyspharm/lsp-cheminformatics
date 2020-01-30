from functools import partial
from typing import Mapping, List, Any, Tuple, Optional

from rdkit.Chem import Mol
from rdkit import DataStructs
from rdkit.Chem import AllChem

try:
    import chemfp

    chemfp_available = True
    import chemfp.search
    import chemfp.arena
    import chemfp.rdkit_types
except ImportError:
    chemfp_available = False

fingerprint_functions = {
    "topological": AllChem.RDKFingerprint,
    "morgan": partial(AllChem.GetMorganFingerprintAsBitVect, radius=2),
}

if chemfp_available:

    def chemfp_fp_func(fun, **kwargs):
        def fp_func(arg_dict):
            for k, v in kwargs.items():
                arg_dict.setdefault(k, v)
            return fun(arg_dict)

        return fp_func

    chemfp_fingerprint_functions = {
        "topological": chemfp_fp_func(
            chemfp.rdkit_types.RDKitFingerprintType_v2, fpSize=2048
        ),
        "morgan": chemfp_fp_func(
            chemfp.rdkit_types.RDKitMorganFingerprintType_v1, fpSize=2048, radius=2
        ),
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


def make_fingerprint_arena(
    mols: Mapping[str, Mol],
    fingerprint_type: str = "morgan",
    fingerprint_args: Mapping[str, Any] = {},
) -> chemfp.arena.FingerprintArena:
    fp_maker = chemfp_fingerprint_functions[fingerprint_type](
        fingerprint_args
    ).make_fingerprinter()
    fp_generator = ((str(n), fp_maker(m)) for n, m in mols.items())
    fp = fp_maker(next(iter(mols.values())))
    arena = chemfp.load_fingerprints(
        fp_generator, metadata=chemfp.Metadata(num_bits=len(fp)*8)
    )
    return arena


def find_similarity_matches(
    query: chemfp.arena.FingerprintArena,
    target: Optional[chemfp.arena.FingerprintArena],
    threshold: float = 0.95,
) -> Mapping[str, Mapping[str, float]]:
    if target is not None:
        match_res = chemfp.search.threshold_tanimoto_search_arena(
            query, target, threshold=threshold
        )
    else:
        match_res = chemfp.search.threshold_tanimoto_search_symmetric(
            query, threshold=threshold, include_lower_triangle=False
        )
    return {
        q_id: {t: s for t, s in targets}
        for q_id, targets in zip(query.ids, match_res.iter_ids_and_scores())
        if len(targets) > 0
    }
