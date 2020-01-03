from flask import request
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Draw, inchi

from lspcheminf import app
from .schemas import DrawGridSchema, DrawGridResultSchema
from .draw import draw_molecule_grid
from .util import convert_compound_request


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
    molecules, skipped = convert_compound_request(data["compounds"])
    svg = draw_molecule_grid(list(molecules.values()), names=list(molecules.keys()))
    out = {"svg": svg, "skipped": skipped}
    DrawGridResultSchema().validate(out)
    return out
