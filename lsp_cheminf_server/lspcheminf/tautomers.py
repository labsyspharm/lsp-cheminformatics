import io
from typing import Mapping, List, Any, Tuple, Optional, Callable
from functools import partial
from concurrent.futures import TimeoutError

from pebble import concurrent, ProcessExpired
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Draw, inchi

from molvs import standardize as mol_standardize
from molvs import tautomer

from chembl_structure_pipeline import (
    standardize_mol as standardize_chembl,
    get_parent_mol,
)

molvs_standardizer = mol_standardize.Standardizer(prefer_organic=True)


def standardize_molvs(mol: Chem.Mol):
    st = molvs_standardizer.standardize(mol)
    st = molvs_standardizer.charge_parent(st, skip_standardize=True)
    st = molvs_standardizer.isotope_parent(st, skip_standardize=True)
    st = molvs_standardizer.standardize(st)
    return st


def standardize_chembl_parent(mol: Chem.Mol):
    return get_parent_mol(mol, check_exclusion=False)[0]


STANDARDIZERS: Mapping[str, Callable] = {
    "molvs": standardize_molvs,
    "chembl": partial(standardize_chembl, check_exclusion=False),
    "chembl-parent": standardize_chembl_parent,
}


def make_tautomers(mol, max_tautomers=10):
    enum = tautomer.TautomerEnumerator(max_tautomers=max_tautomers)
    return enum(mol)


def canonicalize_compound(
    mol: Chem.Mol, canonicalizer_fun: Callable, standardizer_fun: Optional[Callable],
):
    if standardizer_fun is not None:
        try:
            mol = standardizer_fun(mol)
        except Exception as e:
            pass
    try:
        can = canonicalizer_fun(mol)
    except Exception as e:
        raise RuntimeError(f"Could not canonicalize\n{e}")
    if standardizer_fun is not None:
        try:
            can = standardizer_fun(can)
        except Exception as e:
            pass
    return can


def canonicalize(
    compounds: Mapping[str, Chem.Mol],
    standardize: bool = False,
    standardizer: str = "chembl",
    progress_callback: Optional[Callable] = None,
    timeout: Optional[int] = None,
) -> Tuple[Mapping[str, Chem.Mol], List[str]]:
    @concurrent.process(timeout=timeout)
    def process_compound(*args, **kwargs):
        return canonicalize_compound(*args, **kwargs)

    canonicalizer_fun = tautomer.TautomerCanonicalizer()
    res = {}
    skipped = list()
    standardizer_fun = STANDARDIZERS[standardizer]
    # Suppress pesky warning messages
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.ERROR)
    for i, (k, mol) in enumerate(compounds.items()):
        if i % 100 == 0 and progress_callback is not None:
            progress_callback(i)
        future = process_compound(
            mol, canonicalizer_fun=canonicalizer_fun, standardizer_fun=standardizer_fun
        )
        try:
            can = future.result()
        except TimeoutError as error:
            print(f"Processing `{k}` took longer than {timeout}s. Skipping.")
            skipped.append(k)
        except Exception as error:
            print(f"Error canonicalizing {k}. Skipping.\n{error}")
            skipped.append(k)
        res[k] = can
    return (res, skipped)


def tautomerize(mol, max_tautomers=10):
    tauts = make_tautomers(mol, max_tautomers=max_tautomers)
    return tauts
