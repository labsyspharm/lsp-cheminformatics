from __future__ import print_function
import argparse
import pickle
import os
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

@app.route("/fingerprints/fingerprint_db", methods=["POST"])
def fingerprint_arena():
    cmpd_json = request.json["compounds"]
    cmpd_encoding = request.json["request"]["encoding"]
    input_df = pd.DataFrame(cmpd_json)
    print(input_df.head())
    mol_func = identifier_mol_mapping[cmpd_encoding]
    temp_sdf = tempfile.NamedTemporaryFile(mode="wb", suffix=".sdf.gz", delete=False)
    temp_sdf_gz = gzip.GzipFile(fileobj=temp_sdf)
    sdf_writer = Chem.SDWriter(temp_sdf_gz)
    for cmpd in input_df.itertuples():
      m = mol_func(str(cmpd.compound))
      m.SetProp("name", str(cmpd.name))
      sdf_writer.write(m)
    sdf_writer.close()
    temp_sdf_gz.close()
    temp_sdf.close()
    print("Wrote molecules to sdf", temp_sdf.name)
    temp_fps = tempfile.mkstemp(suffix=".fps")
    rdkit2fps.main(
      [
        "-o", temp_fps[1],
        "--id-tag", "name",
        temp_sdf.name,
      ]
    )
    os.remove(temp_sdf.name)
    return send_file(
      temp_fps[1]
    )
