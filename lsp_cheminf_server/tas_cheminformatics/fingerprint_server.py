import base64
import gzip
import os
import pickle
import sys
import tempfile
import pandas as pd
from flask import Flask, request, send_file
from rdkit.Chem import AllChem
from rdkit.Chem import inchi

from tas_cheminformatics import app
from tas_cheminformatics.util import convert_compound_request
from tas_cheminformatics.fingerprint import calculate_similarity
from tas_cheminformatics.schemas import SimilaritySchema, SimilarityResultSchema


@app.route("/fingerprints/similarity", methods=["GET", "POST"])
def similarity_route():
    """Calculate chemical similarity between compounds using chemical fingerprinting.
    ---
    post:
      summary: Calculate chemical similarity
      description: Calculate similarity between query and target compounds using chemical fingerprinting.
      requestBody:
        required: true
        content:
          application/json:
            schema: SimilaritySchema
      responses:
        '200':
          content:
            application/json:
              schema: SimilarityResultSchema
    """
    data = SimilaritySchema().load(request.json)
    query = data["query"]
    target = data["target"]
    if len(query["compounds"]) * len(target["compounds"]) > 100000:
        raise ValueError(
            "Must request fewer than 100,000 similarities per request for performance reasons."
        )
    query_mols = convert_compound_request(query)
    target_mols = convert_compound_request(target)
    similarities = {
        qn: calculate_similarity(
            q,
            target_mols,
            fingerprint_type=data["fingerprint_type"],
            # **data.get(data["fingerprint_args"], {}),
        )
        for qn, q in query_mols.items()
    }
    out = {"query": [], "target": [], "score": []}
    for k, v in similarities.items():
        out["query"].extend([k] * len(v))
        out["target"].extend(v.keys())
        out["score"].extend(v.values())
    SimilarityResultSchema().validate(out)
    return out
