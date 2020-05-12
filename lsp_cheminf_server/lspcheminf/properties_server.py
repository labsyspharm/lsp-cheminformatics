from flask import request
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Descriptors import MolWt

from lspcheminf import app
from .schemas import (
    ConvertIDSchema,
    ConvertIdResultSchema,
    CalculateMassSchema,
    CalculateMassResultSchema,
    IsOrganicSchema,
    IsOrganicResultSchema,
)
from .util import convert_compound_request, mol_identifier_mapping


@app.route("/properties/convert", methods=["POST"])
def convert_ids():
    """Convert compound identifiers.
    ---
    post:
      summary: Convert compound identifiers.
      requestBody:
        required: true
        content:
          application/json:
            schema: ConvertIDSchema
      responses:
        '200':
          content:
            application/json:
              schema: ConvertIdResultSchema
    """
    data = ConvertIDSchema().load(request.json)
    compounds, skipped = convert_compound_request(data["compounds"])
    mol_out_mapping = mol_identifier_mapping[data["target_identifier"]]
    compounds_out = {
        "compounds": [],
        "names": [],
        "identifier": data["target_identifier"],
    }
    for n, m in compounds.items():
        try:
            compounds_out["compounds"].append(mol_out_mapping(m))
            compounds_out["names"].append(n)
        except Exception as e:
            skipped.append(mol_identifier_mapping[data["compounds"]["identifier"]](m))
    out = {"compounds": compounds_out, "skipped": skipped}
    ConvertIdResultSchema().validate(out)
    return out


@app.route("/properties/mass", methods=["POST"])
def calculate_mass_route():
    """Calculate compound molecular mass.
    ---
    post:
      summary: Calculate compound molecular mass.
      requestBody:
        required: true
        content:
          application/json:
            schema: CalculateMassSchema
      responses:
        '200':
          content:
            application/json:
              schema: CalculateMassResultSchema
    """
    data = CalculateMassSchema().load(request.json)
    compounds, skipped = convert_compound_request(data["compounds"])
    mass_out = {}
    for n, m in compounds.items():
        try:
            mass_out[n] = MolWt(m)
        except Exception as e:
            skipped.append(n)
    out = {"mass": mass_out, "skipped": skipped}
    CalculateMassResultSchema().validate(out)
    return out


@app.route("/properties/organic", methods=["POST"])
def is_organic_route():
    """Check if compound is organic. This is true if it contains at least one carbon.
    ---
    post:
      summary: Check if compound is organic.
      requestBody:
        required: true
        content:
          application/json:
            schema: IsOrganicSchema
      responses:
        '200':
          content:
            application/json:
              schema: IsOrganicResultSchema
    """
    data = CalculateMassSchema().load(request.json)
    compounds, skipped = convert_compound_request(data["compounds"])
    organic_out = {}
    for n, m in compounds.items():
        try:
            organic_out[n] = any(a.GetSymbol() == "C" for a in m.GetAtoms())
        except Exception as e:
            skipped.append(n)
    out = {"organic": organic_out, "skipped": skipped}
    CalculateMassResultSchema().validate(out)
    return out
