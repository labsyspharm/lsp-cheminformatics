from flask import request
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Draw, inchi

from tas_chemoinformatics import app
from .tautomers import canonicalize, tautomerize
from .util import identifier_mol_mapping, mol_identifier_mapping

# Run using:
# FLASK_APP=id_mapping/tautomer_server.py flask run -p 5000
# OR to support multiple parallel requests to speed things up:
# gunicorn --workers=4 -b 127.0.0.1:5000 -t 600 tas_chemoinformatics.tautomer_server:app


@app.route("/query/tautomers", methods=["POST"])
def tautomerize_route():
    ids_used = [k for k in identifier_mol_mapping.keys() if k in request.form]
    if not len(ids_used) == 1:
        raise ValueError(
            "Request needs to contain exactly one of these identifies: ",
            repr(list(identifier_mol_mapping.keys())),
        )
    id_used = ids_used[0]
    mol_input = request.form[id_used]
    print(f"Requested tautomers for {id_used} {mol_input}")
    mol = identifier_mol_mapping[id_used](mol_input)
    mol_smiles = Chem.MolToSmiles(mol)
    mol_inchi = inchi.MolToInchi(mol)
    mol_inchi_key = inchi.MolToInchiKey(mol)
    tauts = tautomerize(mol_input, id_used, request.form.get("max_tautomers", 10))
    print(f"Found tautomers: {len(tauts)}")
    return {
        "request": {
            "smiles": mol_smiles,
            "inchi": mol_inchi,
            "inchi_key": mol_inchi_key,
        },
        "tautomers": tauts,
    }


@app.route("/query/canonicalize", methods=["POST"])
def canonicalize_route():
    ids_used = [k for k in identifier_mol_mapping.keys() if k in request.json]
    if not len(ids_used) == 1:
        raise ValueError(
            "Request needs to contain exactly one of these identifies: ",
            repr(list(identifier_mol_mapping.keys())),
        )
    id_used = ids_used[0]
    mol_input = request.json[id_used]
    standardize = bool(request.json.get("standardize", False))
    if not isinstance(mol_input, list):
        mol_input = [mol_input]
    print(f"Requested canonical tautomer for {id_used} input")
    res, skipped = canonicalize(mol_input, id_used, standardize=standardize)
    res = list(zip(*res))
    return {
        "canonicalized": {"query": res[0], "smiles": res[1], "inchi": res[2]},
        "skipped": skipped,
    }


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
    print(f"Requested conversion of {format_in} to {format_out}")
    mols_out = dict()
    for m in ids_in:
        try:
            mol = identifier_mol_mapping[format_in](m)
            mol_out = mol_identifier_mapping[format_out](mol)
            mols_out[m] = mol_out
        except Exception as e:
            print(f"Can't convert {m}: {str(e)}")
            mols_out[m] = None
    return mols_out
