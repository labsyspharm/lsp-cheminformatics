from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Draw, inchi

from molvs import standardize as mol_standardize
from molvs import tautomer


def draw_molecule_grid(compounds, names=None, common_core=None, draw_args={}):
    if common_core is not None:
        AllChem.Compute2DCoords(common_core)
    for m in compounds:
        if common_core is not None:
            AllChem.GenerateDepictionMatching2DStructure(
                m, common_core, acceptFailure=True
            )
        else:
            AllChem.Compute2DCoords(m)
    if names is not None:
        names = [str(n) for n in names]
    draw_args.setdefault("molsPerRow", 4)
    draw_args.setdefault("subImgSize", (200, 200))
    img = Draw.MolsToGridImage(compounds, legends=names, useSVG=True, **draw_args)
    return img
