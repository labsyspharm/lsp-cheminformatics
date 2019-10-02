from __future__ import print_function

import argparse
import base64
import gzip
import os
import pickle
import sys
import tempfile
from tempfile import NamedTemporaryFile

try:
    from itertools import izip
except ImportError:
    izip = zip

import pandas as pd
from flask import Flask, request, send_file
from rdkit import Chem
from rdkit.Chem import inchi

from tas_chemoinformatics import app
from tas_chemoinformatics.util import identifier_mol_mapping

# from tas_chemoinformatics.schemas import FingerprintDBSchema, FingerprintDBResultSchema

try:
    from chemfp.commandline import rdkit2fps
    from chemfp.commandline import simsearch
except ImportError:
    pass

# Run using:
# FLASK_APP=id_mapping/fingerprint_server.py flask run -p 8000
# OR to support multiple parallel requests to speed things up:
# PYTHONPATH=id_mapping gunicorn --workers=4 -b 127.0.0.1:8000 -t 600 fingerprint_server:app

fp_types = {"topological": ["--RDK"], "morgan": ["--morgan"]}


@app.route("/fingerprints/fingerprint_db", methods=["POST"])
def fingerprint_db():
    """Create database of fingerprints from given compounds.
    ---
    post:
      summary: Create database of fingerprints from given compounds.
      description: Fingerprint database is returned as base64 encoded database file in fps format.
      requestBody:
        required: true
        content:
          application/json:
            schema: FingerprintDBSchema
      responses:
        '200':
          content:
            application/json:
              schema: FingerprintDBResultSchema
    """
    data = request.json
    cmpd_data = data["compounds"]
    id_used = cmpd_data["identifier"]
    fp_type = data.get("fingerprint_type", "topological")
    fp_args = data.get("fingerprint_args", [])
    if not isinstance(fp_args, list):
        fp_args = [fp_args]
    cmpd_names = cmpd_data.get("names", list(range(len(cmpd_data["compounds"]))))
    print("Processing {} compound fingerprints".format(len(cmpd_data["compounds"])))
    if len(cmpd_data["compounds"]) > 100000:
        raise ValueError(
            "No more than 100,000 compounds per request for performance reasons"
        )
    mol_func = identifier_mol_mapping[id_used]
    temp_sdf = NamedTemporaryFile(mode="wb", suffix=".sdf.gz", delete=False)
    temp_sdf_gz = gzip.GzipFile(fileobj=temp_sdf)
    sdf_writer = Chem.SDWriter(temp_sdf_gz)
    skipped = list()
    for i, (n, c) in enumerate(izip(cmpd_names, cmpd_data["compounds"])):
        if i % 1000 == 0:
            print("{} done".format(i))
        m = mol_func(str(c), sanitize=True)
        if m is None:
            skipped.append(c)
            continue
        m.SetProp("name", str(n))
        sdf_writer.write(m)
    sdf_writer.close()
    temp_sdf_gz.close()
    temp_sdf.close()
    print("Wrote molecules to sdf", temp_sdf.name)
    temp_fps = tempfile.mkstemp(suffix=".fps")
    rdkit_cmd = ["-o", temp_fps[1], "--id-tag", "name"]
    rdkit_cmd.extend(fp_types[fp_type])
    rdkit_cmd.extend(fp_args)
    rdkit2fps.main(rdkit_cmd + [temp_sdf.name])
    os.remove(temp_sdf.name)
    with open(temp_fps[1], "rb") as f:
        fps_b64 = base64.b64encode(f.read())
    os.remove(temp_fps[1])
    out = {"fingerprint_db": fps_b64, "skipped": skipped}
    # FingerprintDBResultSchema().validate(out)
    return out


def pairwise(iterable):
    "s -> (s0, s1), (s2, s3), (s4, s5), ..."
    a = iter(iterable)
    return izip(a, a)


def parse_sim_result(sim_result_file):
    res = list()
    with open(sim_result_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            l = line.split()
            if len(l) <= 2:
                continue
            query = l[1]
            for match, score in pairwise(l[2:]):
                res.append((query, match, score))
    return zip(*res)


def fps_newline(data):
    if data.endswith("\n"):
        return data
    return data + "\n"


@app.route("/fingerprints/simsearch", methods=["POST"])
def fingerprint_search():
    """Scan database of fingerprints for matches.
    ---
    post:
      summary: Scan database of fingerprints for matches.
      requestBody:
        required: true
        content:
          application/json:
            schema: FingerprintScanSchema
      responses:
        '200':
          content:
            application/json:
              schema: FingerprintScanResultSchema
    """
    data = request.json
    fingerprint_db = data["fingerprint_db"]
    fingerprint_query = data.get("fingerprint_query", None)
    threshold = data.get("threshold", 0.95)
    print("Starting processing")
    with NamedTemporaryFile(suffix=".fps") as temp_db, NamedTemporaryFile(
        suffix=".txt"
    ) as temp_out, NamedTemporaryFile(suffix=".fps") as temp_query:
        print(
            "temp_db {}\ntemp_out {}\ntemp_query {}".format(
                temp_db.name, temp_out.name, temp_query.name
            )
        )
        temp_db.write(fps_newline(base64.b64decode(fingerprint_db)))
        args = []
        if fingerprint_query:
            temp_query.write(fps_newline(base64.b64decode(fingerprint_query)))
            args.extend(["-q", temp_query.name])
        else:
            args.append("--NxN")
        args.extend(["--memory", "-o", temp_out.name, "-t", threshold, temp_db.name])
        args = map(str, args)
        print(args)
        # Make sure files are completely written to disk
        for f in [temp_db, temp_out, temp_query]:
            f.flush()
            os.fsync(f.fileno())
        simsearch.main(args)
        sim_results = parse_sim_result(temp_out.name)
    return {"query": sim_results[0], "match": sim_results[1], "score": sim_results[2]}
