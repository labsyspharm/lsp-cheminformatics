import argparse
import pickle
import sys
import gzip
import tempfile
import rdkit
import pandas as pd
from flask import Flask, request, send_file
from rdkit import Chem
from rdkit.Chem import inchi
from rdkit.Chem.Fingerprints import FingerprintMols
from chemfp.commandline import rdkit2fps

identifier_mol_mapping = {"smiles": Chem.MolFromSmiles, "inchi": inchi.MolFromInchi}

app = Flask(__name__)

@app.route("/fingerprints/fingerprint_db")
def fingerprint_arena():
    cmpd_json = request.json["compounds"]
    cmpd_encoding = request.json["request"]["encoding"]
    input_df = pd.read_json(cmpd_json)
    print(input_df.head())
    mol_func = identifier_mol_mapping[cmpd_encoding]
    temp_sdf = tempfile.NamedTemporaryFile(mode="wb", suffix=".sdf.gz")
    temp_sdf_gz = gzip.GzipFile(fileobj=temp_sdf)
    sdf_writer = Chem.SDWriter(temp_sdf_gz)
    for cmpd in input_df.iterrows():
      m = mol_func(cmpd["compound"])
      m.SetProp("name", cmpd["name"])
      sdf_writer.write(m)
    print("Wrote molecules to sdf")
    temp_fps = tempfile.mkstemp(suffix=".fps")
    rdkit2fps.main(
      [
        temp_sdf.name,
        "--id-tag", "name",
        "-o", temp_fps[1],
      ]
    )
    return send_file(
      temp_fps[1]
    )
