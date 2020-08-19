import io
from typing import Mapping, List, Any, Tuple, Optional, Callable
from functools import partial

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Draw, inchi

from molvs import standardize as mol_standardize
from molvs import tautomer

from chembl_structure_pipeline import standardize_mol as standardize_chembl

molvs_standardizer = mol_standardize.Standardizer(prefer_organic=True)


def standardize_molvs(mol: Chem.Mol):
    st = molvs_standardizer.standardize(mol)
    st = molvs_standardizer.charge_parent(st, skip_standardize=True)
    st = molvs_standardizer.isotope_parent(st, skip_standardize=True)
    st = molvs_standardizer.standardize(st)
    return st


STANDARDIZERS: Mapping[str, Callable] = {
    "molvs": standardize_molvs,
    "chembl": partial(standardize_chembl, check_exclusion=False),
}


def make_tautomers(mol, max_tautomers=10):
    enum = tautomer.TautomerEnumerator(max_tautomers=max_tautomers)
    return enum(mol)


def canonicalize(
    compounds: Mapping[str, Chem.Mol], standardize: bool = False, standardizer="chembl"
) -> Tuple[Mapping[str, Chem.Mol], List[str]]:
    canonicalizer = tautomer.TautomerCanonicalizer()
    res = {}
    skipped = list()
    standardizer_fun = STANDARDIZERS[standardizer]
    # Suppress pesky warning messages
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.ERROR)
    for k, mol in compounds.items():
        if standardize:
            try:
                mol = standardizer_fun(mol)
            except Exception as e:
                print(f"Can't standardize {k}, using input unstandardized\n{e}")
        try:
            can = canonicalizer(mol)
        except Exception as e:
            print(f"Skipping {k}: Could not canonicalize\n{e}")
            skipped.append(k)
            continue
        if standardize:
            try:
                can = standardizer_fun(can)
            except Exception as e:
                print(
                    f"Can't standardize canonical tautomer for {k}, using unstandardized version\n{e}"
                )
        res[k] = can
    return (res, skipped)


def tautomerize(mol, max_tautomers=10):
    tauts = make_tautomers(mol, max_tautomers=max_tautomers)
    return tauts
