from flask import request
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Draw, inchi

from tas_cheminformatics import app
from .schemas import (
    TautomerizeSchema,
    CanonicalizeSchema,
    TautomerizeResultSchema,
    CanonicalizeResultSchema,
)
from .tautomers import canonicalize, tautomerize
from .util import identifier_mol_mapping, mol_identifier_mapping


@app.route("/tautomers/enumerate", methods=["POST"])
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
                identifier_mol_mapping[id_used](in_id), data.get("max_tautomers", 10)
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


@app.route("/tautomers/canonicalize", methods=["POST"])
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
    # print(f"Requested canonical tautomer for {id_used} input")
    res, skipped = canonicalize(mol_input, id_used, standardize=data["standardize"])
    out = {
        "canonical": {
            k: {"identifier": id_used, "compounds": mol_identifier_mapping[id_used](v)}
            for k, v in res.items()
        },
        "skipped": skipped,
    }
    CanonicalizeResultSchema().validate(out)
    return out
