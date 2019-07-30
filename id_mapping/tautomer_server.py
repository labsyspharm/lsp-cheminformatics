from flask import Flask, request
from rdkit import Chem
from rdkit.Chem import inchi
from molvs import tautomer

# Run by setting
# export FLASK_APP=id_mapping/tautomer_server.py
# flask run

def make_tautomers(mol, max_tautomers = 10):
    enum = tautomer.TautomerEnumerator(max_tautomers = max_tautomers)
    return enum(mol)
    
identifier_mol_mapping = {
    "smiles": Chem.MolFromSmiles,
    "inchi": inchi.MolFromInchi,
}

mol_identifier_mapping = {
    "smiles": Chem.MolToSmiles,
    "inchi": inchi.MolToInchi,
    "inchi_key": inchi.MolToInchiKey,
}

app = Flask(__name__)

@app.route("/query/tautomers", methods=["POST"])
def process_tautomers():
    ids_used = [k for k in identifier_mol_mapping.keys() if k in request.form]
    if not len(ids_used) == 1:
        raise ValueError(
            "Request needs to contain exactly one of these identifies: ",
            repr(list(identifier_mol_mapping.keys()))
        )
    id_used = ids_used[0]
    mol_input = request.form[id_used]
    print(f"Requested tautomers for {id_used} {mol_input}")
    mol = identifier_mol_mapping[id_used](mol_input)
    mol_smiles = Chem.MolToSmiles(mol)
    mol_inchi = inchi.MolToInchi(mol)
    mol_inchi_key = inchi.MolToInchiKey(mol)
    max_tautomers = request.form.get("max_tautomers", 10)
    tauts = make_tautomers(mol, max_tautomers = max_tautomers)
    print(f"Found tautomers: {len(tauts)}")
    return {
        "request": {
            "smiles": mol_smiles,
            "inchi": mol_inchi,
            "inchi_key": mol_inchi_key,
        },
        "tautomers": [
            {
                "smiles": Chem.MolToSmiles(t),
                "inchi": inchi.MolToInchi(t),
                "inchi_key": inchi.MolToInchiKey(t),
            }
            for t in tauts
        ],
    }
    
@app.route("/query/convert", methods=["POST"])
def convert_ids():
    format_in = request.form["in"]
    if format_in not in identifier_mol_mapping:
        raise ValueError("input must be one of: ", repr(list(identifier_mol_mapping.keys)))
    format_out = request.form["out"]
    if format_out not in mol_identifier_mapping:
        raise ValueError("output must be one of: ", repr(list(mol_identifier_mapping.keys)))
    id_in = request.form["value"]
    print(f"Requested conversion of {format_in} {id_in} to {format_out}")
    mol = identifier_mol_mapping[format_in](id_in)
    mol_out = mol_identifier_mapping[format_out](mol)
    return {
        "request": {
            "in": format_in,
            "out": format_out,
            "value": id_in,
        },
        "converted": {
            format_out: mol_out,
        },
    }
