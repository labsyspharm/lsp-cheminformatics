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


@app.route("/query/tautomers", methods=["POST"])
def tautomerize_route():
    """Enumerate tautomers for the given compounds.
    ---
    post:
      summary: Enumerate tautomers for the given compounds.
      requestBody:
        required: true
        content:
          application/json:
            schema: TautomerizeSchema
      responses:
        '200':
          content:
            application/json:
              schema: TautomerizeResultSchema
    """
    data = TautomerizeSchema().load(request.json)
    id_used = data["compounds"]["identifier"]
    # print(f"Requested tautomers for {id_used} {mol_input}")
    out_tauts = {}
    for in_id in data["compounds"]["compounds"]:
        try:
            tauts = tautomerize(
                identifier_mol_mapping[id_used](in_id),
                data.get("max_tautomers", 10)
            )
        except Exception as e:
            print(e)
            continue
        print("found", len(tauts))
        out_tauts[in_id] = [mol_identifier_mapping[id_used](t) for t in tauts]
    out = {"tautomers": out_tauts}
    print("totl", len(out_tauts))
    # print(f"Found tautomers: {len(tauts)}")
    TautomerizeResultSchema().validate(out)
    return out


@app.route("/query/canonicalize", methods=["POST"])
def canonicalize_route():
    """Find canonical representation of compounds.
    ---
    post:
      summary: Find canonical representation of compounds.
      description: >-
        Find canonical tautomer for all compounds in the `COMPOUND_FILE` input.
        This file should be a csv file where one column contains compounds as InChI or SMILES.
        The canonical compounds will be written as InChI and SMILES to the `OUTPUT_FILE`.
      requestBody:
        required: true
        content:
          application/json:
            schema: CanonicalizeSchema
      responses:
        '200':
          content:
            application/json:
              schema: CanonicalizeResultSchema
    """
    data = CanonicalizeSchema().load(request.json)
    id_used = data["compounds"]["identifier"]
    mol_input = data["compounds"]["compounds"]
    print(f"Requested canonical tautomer for {id_used} input")
    res, skipped = canonicalize(mol_input, id_used, standardize=data["standardize"])
    out = {
        "canonical": {
            k: {"identifier": id_used, "compounds": mol_identifier_mapping[id_used](v)} for k, v in res.items()
        },
        "skipped": skipped,
    }
    CanonicalizeResultSchema().validate(out)
    return out


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
