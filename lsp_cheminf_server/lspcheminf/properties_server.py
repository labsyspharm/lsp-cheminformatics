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
)
from .util import identifier_mol_mapping, mol_identifier_mapping


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
    # print(f"Requested conversion of {format_in} to {format_out}")
    mol_in_mapping = identifier_mol_mapping[data["compounds"]["identifier"]]
    mol_out_mapping = mol_identifier_mapping[data["target_identifier"]]
    skipped = []
    compounds_out = {"compounds": [], "identifier": data["target_identifier"]}
    for m in data["compounds"]["compounds"]:
        try:
            mol = mol_in_mapping(m)
            compounds_out["compounds"].append(mol_out_mapping(mol))
        except Exception as e:
            skipped.append(m)
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
    # print(f"Requested conversion of {format_in} to {format_out}")
    mol_in_mapping = identifier_mol_mapping[data["compounds"]["identifier"]]
    skipped = []
    mass_out = {}
    for m in data["compounds"]["compounds"]:
        try:
            mol = mol_in_mapping(m)
            mass_out[m] = MolWt(mol)
        except Exception as e:
            skipped.append(m)
    out = {"mass": mass_out, "skipped": skipped}
    CalculateMassResultSchema().validate(out)
    return out
