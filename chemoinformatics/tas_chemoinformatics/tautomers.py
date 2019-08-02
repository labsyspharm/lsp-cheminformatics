from flask import Flask, request
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi
from molvs import tautomer

# Run using:
# FLASK_APP=id_mapping/tautomer_server.py flask run -p 5000
# OR to support multiple parallel requests to speed things up:
# PYTHONPATH=id_mapping gunicorn --workers=4 -b 127.0.0.1:5000 -t 600 tautomer_server:app


def make_tautomers(mol, max_tautomers=10):
    enum = tautomer.TautomerEnumerator(max_tautomers=max_tautomers)
    return enum(mol)


identifier_mol_mapping = {"smiles": Chem.MolFromSmiles, "inchi": inchi.MolFromInchi}

mol_identifier_mapping = {
    "smiles": Chem.MolToSmiles,
    "inchi": inchi.MolToInchi,
    "inchi_key": inchi.MolToInchiKey,
}


def canonicalize(compounds, id_used):
    canonicalizer = tautomer.TautomerCanonicalizer()
    mol_mapping = identifier_mol_mapping[id_used]
    res = list()
    skipped = list()
    # Suppress pesky warning messages
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.ERROR)
    for i, ms in enumerate(compounds):
        if i % 1000 == 0:
            print(f"Completed {i}")
        mol = mol_mapping(ms)
        if mol is None:
            print(f"Skipping {ms}: Could not parse compound string")
            skipped.append(ms)
            continue
        try:
            can = canonicalizer(mol)
        except Exception as e:
            print(f"Skipping {ms}: Could not canonicalize\n{e}")
            skipped.append(ms)
            continue
        can_smiles = Chem.MolToSmiles(mol)
        can_inchi = inchi.MolToInchi(mol)
        res.append((ms, can_smiles, can_inchi))
    return (res, skipped)


def tautomerize(compound, id_used, max_tautomers=10):
    mol = identifier_mol_mapping[id_used](compound)
    max_tautomers = request.form.get("max_tautomers", 10)
    tauts = make_tautomers(mol, max_tautomers=max_tautomers)
    return [
        {
            "smiles": Chem.MolToSmiles(t),
            "inchi": inchi.MolToInchi(t),
            "inchi_key": inchi.MolToInchiKey(t),
        }
        for t in tauts
    ]


app = Flask(__name__)


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
    if not isinstance(mol_input, list):
        mol_input = [mol_input]
    print(f"Requested canonical tautomer for {id_used} input")
    res, skipped = canonicalize(mol_input, id_used)
    return {
        "canonicalized": {"query": res[0], "smiles": res[1], "inchi": res[2]},
        "skipped": skipped,
    }


@app.route("/query/convert", methods=["POST"])
def convert_ids():
    format_in = request.form["in"]
    if format_in not in identifier_mol_mapping:
        raise ValueError(
            "input must be one of: ", repr(list(identifier_mol_mapping.keys))
        )
    format_out = request.form["out"]
    if format_out not in mol_identifier_mapping:
        raise ValueError(
            "output must be one of: ", repr(list(mol_identifier_mapping.keys))
        )
    id_in = request.form["value"]
    print(f"Requested conversion of {format_in} {id_in} to {format_out}")
    mol = identifier_mol_mapping[format_in](id_in)
    mol_out = mol_identifier_mapping[format_out](mol)
    return {
        "request": {"in": format_in, "out": format_out, "value": id_in},
        "converted": {format_out: mol_out},
    }
