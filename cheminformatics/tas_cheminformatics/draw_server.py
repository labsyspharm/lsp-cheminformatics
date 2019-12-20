from flask import request
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Draw, inchi

from tas_cheminformatics import app
from .schemas import DrawGridSchema, DrawGridResultSchema
from .draw import draw_molecule_grid
from .util import identifier_mol_mapping, mol_identifier_mapping

# Run using:
# FLASK_APP=__init__.py flask run -p 5000
# OR to support multiple parallel requests to speed things up:
# gunicorn --workers=4 -b 127.0.0.1:5000 -t 600 tas_cheminformatics


@app.route("/draw/grid", methods=["POST"])
def draw_molecule_grid_route():
    """Draw grid of molecules.
    ---
    post:
      summary: Draw grid of molecules.
      requestBody:
        required: true
        content:
          application/json:
            schema: DrawGridSchema
      responses:
        '200':
          content:
            application/json:
              schema: DrawGridResultSchema
    """
    data = DrawGridSchema().load(request.json)
    # print(f"Requested conversion of {format_in} to {format_out}")
    names = data["compounds"].get("names", list(range(len(data["compounds"]["compounds"]))))
    mol_in_mapping = identifier_mol_mapping[data["compounds"]["identifier"]]
    skipped = []
    mols = {}
    for n, m in zip(names, data["compounds"]["compounds"]):
        try:
            mol = mol_in_mapping(m)
            mols[n] = mol
        except Exception as e:
            skipped.append(m)
    svg = draw_molecule_grid(list(mols.values()), names=list(mols.keys()))
    out = {"svg": svg, "skipped": skipped}
    DrawGridResultSchema().validate(out)
    return out
