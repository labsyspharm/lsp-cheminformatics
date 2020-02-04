import io
from typing import Mapping, List, Any, Tuple, Optional

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Draw, inchi

from molvs import standardize as mol_standardize
from molvs import tautomer


def make_tautomers(mol, max_tautomers=10):
    enum = tautomer.TautomerEnumerator(max_tautomers=max_tautomers)
    return enum(mol)


def canonicalize(
    compounds: Mapping[str, Chem.Mol], standardize: bool = False
) -> Tuple[Mapping[str, Chem.Mol], List[str]]:
    canonicalizer = tautomer.TautomerCanonicalizer()
    if standardize:
        standardizer = mol_standardize.Standardizer(prefer_organic=True)
    res = {}
    skipped = list()
    # Suppress pesky warning messages
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.ERROR)
    for k, mol in compounds.items():
        if standardize:
            try:
                st = standardizer.standardize(mol)
                st = standardizer.charge_parent(st, skip_standardize=True)
                st = standardizer.isotope_parent(st, skip_standardize=True)
                st = standardizer.standardize(st)
                mol = st
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
                can = standardizer.standardize(can)
            except Exception as e:
                print(
                    f"Can't standardize canonical tautomer for {k}, using unstandardized version\n{e}"
                )
        res[k] = can
    return (res, skipped)


def tautomerize(mol, max_tautomers=10):
    tauts = make_tautomers(mol, max_tautomers=max_tautomers)
    return tauts
