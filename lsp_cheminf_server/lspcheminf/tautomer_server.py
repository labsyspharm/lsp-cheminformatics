from flask import request
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Draw, inchi

from lspcheminf import app
from .schemas import (
    TautomerizeSchema,
    CanonicalizeSchema,
    TautomerizeResultSchema,
    CanonicalizeResultSchema,
)
from .tautomers import canonicalize, tautomerize
from .util import (
    identifier_mol_mapping,
    mol_identifier_mapping,
    convert_compound_request,
)


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
    compounds, skipped = convert_compound_request(data["compounds"])
    id_used = data["compounds"]["identifier"]
    out_tauts = {}
    for n, m in compounds.items():
        try:
            tauts = tautomerize(m, data.get("max_tautomers", 10))
        except Exception as e:
            print(e)
            skipped.append(n)
            continue
        print("found", len(tauts))
        out_tauts[n] = {
            "compounds": list(set(mol_identifier_mapping[id_used](t) for t in tauts)),
            "identifier": id_used,
        }
    out = {"tautomers": out_tauts, "skipped": skipped}
    TautomerizeResultSchema().validate(out)
    return out


@app.route("/tautomers/canonicalize", methods=["POST"])
def canonicalize_route():
    """Find canonical representation of compounds.
    ---
    post:
      summary: Find canonical representation of compounds.
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
    compounds, skipped = convert_compound_request(data["compounds"])
    res, skipped = canonicalize(compounds, standardize=data["standardize"])
    mol_to_inchi = mol_identifier_mapping["inchi"]
    out = {
        "canonical": {
            "identifier": "inchi",
            "compounds": {k: mol_to_inchi(mol) for k, mol in res.items()},
        },
        "skipped": skipped,
    }
    CanonicalizeResultSchema().validate(out)
    return out
