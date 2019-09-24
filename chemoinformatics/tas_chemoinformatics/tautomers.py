import io

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Draw, inchi

from molvs import standardize as mol_standardize
from molvs import tautomer

from .util import identifier_mol_mapping


def make_tautomers(mol, max_tautomers=10):
    enum = tautomer.TautomerEnumerator(max_tautomers=max_tautomers)
    return enum(mol)


def canonicalize(compounds, id_used, standardize=False):
    canonicalizer = tautomer.TautomerCanonicalizer()
    mol_mapping = identifier_mol_mapping[id_used]
    if standardize:
        standardizer = mol_standardize.Standardizer(prefer_organic=True)
    res = list()
    skipped = list()
    # Suppress pesky warning messages
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.ERROR)
    for i, ms in enumerate(compounds):
        if i % 1000 == 0:
            print(f"Completed {i}")
        mol = mol_mapping(ms)
        if mol is None:
            print(f"Skipping {ms}: Could not parse compound string")
            skipped.append(ms)
            continue
        if standardize:
            try:
                st = standardizer.standardize(mol)
                st = standardizer.charge_parent(st, skip_standardize=True)
                st = standardizer.isotope_parent(st, skip_standardize=True)
                st = standardizer.standardize(st)
                mol = st
            except Exception as e:
                print(f"Can't standardize {ms}, using input unstandardized\n{e}")
        try:
            can = canonicalizer(mol)
        except Exception as e:
            print(f"Skipping {ms}: Could not canonicalize\n{e}")
            skipped.append(ms)
            continue
        if standardize:
            try:
                can = standardizer.standardize(can)
            except Exception as e:
                print(
                    f"Can't standardize canonical tautomer for {ms}, using unstandardized version\n{e}"
                )
        can_smiles = Chem.MolToSmiles(can)
        can_inchi = inchi.MolToInchi(can)
        res.append((ms, can_smiles, can_inchi))
    return (res, skipped)


def tautomerize(compound, id_used, max_tautomers=10):
    mol = identifier_mol_mapping[id_used](compound)
    tauts = make_tautomers(mol, max_tautomers=max_tautomers)
    return [
        {
            "smiles": Chem.MolToSmiles(t),
            "inchi": inchi.MolToInchi(t),
            "inchi_key": inchi.MolToInchiKey(t),
        }
        for t in tauts
    ]


def draw_molecules(compounds, id_used, names=None):
    if not isinstance(compounds, list):
        compounds = [compounds]
    mols = [identifier_mol_mapping[id_used](c) for c in compounds]
    for m in mols:
        AllChem.Compute2DCoords(m)
    img = Draw.MolsToGridImage(
        mols, molsPerRow=4, subImgSize=(200, 200), legends=names, useSVG=True
    )
    return img
