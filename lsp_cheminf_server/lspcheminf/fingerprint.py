from functools import partial
from math import isclose
from typing import Mapping, List, Any, Tuple, Optional

from rdkit.Chem import Mol
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Descriptors import MolWt

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
            arg_dict = arg_dict.copy()
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
        fp_generator, metadata=chemfp.Metadata(num_bits=len(fp) * 8)
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


def compound_identity(
    query: Mapping[str, Mol], target: Optional[Mapping[str, Mol]]
) -> Mapping[str, List[str]]:
    target_set = set((target if target is not None else query).keys())
    match_sets = {q: target_set.copy() for q in query.keys()}
    for fp_type in ["morgan", "topological"]:
        query_arena = make_fingerprint_arena(query, fingerprint_type=fp_type)
        target_arena = (
            make_fingerprint_arena(target, fingerprint_type=fp_type)
            if target is not None
            else None
        )
        matches = find_similarity_matches(query_arena, target_arena, threshold=1)
        for q in match_sets.keys():
            match_sets[q] &= set(x for x in matches.get(q, {}).keys())
    query_weights = {k: MolWt(m) for k, m in query.items()}
    target_weights = (
        {k: MolWt(m) for k, m in target.items()}
        if target is not None
        else query_weights
    )
    for q, ts in match_sets.items():
        weight_matches = set()
        for t in ts:
            if isclose(query_weights[q], target_weights[t], rel_tol=0.001):
                weight_matches.add(t)
        match_sets[q] &= weight_matches
    return {k: list(v) for k, v in match_sets.items()}
