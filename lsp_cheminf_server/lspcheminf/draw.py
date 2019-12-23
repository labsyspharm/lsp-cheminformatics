from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Draw, inchi

from molvs import standardize as mol_standardize
from molvs import tautomer


def draw_molecule_grid(compounds, names=None):
    for m in compounds:
        AllChem.Compute2DCoords(m)
    img = Draw.MolsToGridImage(
        compounds, molsPerRow=4, subImgSize=(200, 200), legends=names, useSVG=True
    )
    return img
