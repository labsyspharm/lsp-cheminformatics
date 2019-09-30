from flask import request
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Draw, inchi

from tas_chemoinformatics import app
from .schemas import TautomerizeSchema, CanonicalizeSchema, TautomerizeResultSchema, CanonicalizeResultSchema
from .tautomers import canonicalize, tautomerize
from .util import identifier_mol_mapping, mol_identifier_mapping

# Run using:
# FLASK_APP=__init__.py flask run -p 5000
# OR to support multiple parallel requests to speed things up:
# gunicorn --workers=4 -b 127.0.0.1:5000 -t 600 tas_chemoinformatics


@app.route("/query/convert", methods=["POST"])
def convert_ids():
    format_in = request.json["in"]
    if format_in not in identifier_mol_mapping:
        raise ValueError(
            "input must be one of: ", repr(list(identifier_mol_mapping.keys))
        )
    format_out = request.json["out"]
    if format_out not in mol_identifier_mapping:
        raise ValueError(
            "output must be one of: ", repr(list(mol_identifier_mapping.keys))
        )
    ids_in = request.json["value"]
    if not isinstance(ids_in, list):
        ids_in = [ids_in]
    ids_in = set(ids_in)
    # print(f"Requested conversion of {format_in} to {format_out}")
    mols_out = dict()
    for m in ids_in:
        try:
            mol = identifier_mol_mapping[format_in](m)
            mol_out = mol_identifier_mapping[format_out](mol)
            mols_out[m] = mol_out
        except Exception as e:
            # print(f"Can't convert {m}: {str(e)}")
            mols_out[m] = None
    return mols_out
