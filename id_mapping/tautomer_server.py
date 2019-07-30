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
    mol = identifier_mol_mapping[id_used](mol_input)
    mol_smiles = Chem.MolToSmiles(mol)
    mol_inchi = inchi.MolToInchi(mol)
    mol_inchi_key = inchi.MolToInchiKey(mol)
    print(f"Requested tautomers for {id_used} {mol_input}")
    max_tautomers = request.form.get("tautomers", 10)
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
