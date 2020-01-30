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

from lspcheminf import app
from lspcheminf.util import convert_compound_request
from lspcheminf.fingerprint import (
    calculate_similarity,
    find_substructure_matches,
    find_similarity_matches,
    make_fingerprint_arena,
)
from lspcheminf.schemas import (
    SimilaritySchema,
    SimilarityResultSchema,
    SimilarityThresholdSchema,
    SubstructureSchema,
    SubstructureResultSchema,
)


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
    query_mols, query_skipped = convert_compound_request(query)
    target_mols, target_skipped = convert_compound_request(target)
    similarities = {
        qn: calculate_similarity(
            q,
            target_mols,
            fingerprint_type=data["fingerprint_type"],
            fingerprint_args=data["fingerprint_args"],
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


@app.route("/fingerprints/substructure", methods=["GET", "POST"])
def substructure_route():
    """Find substructures from the queries that also occur in the targets.
    ---
    post:
      summary: Find substructure matches
      description: Find targets whose substructure matches with a query.
      requestBody:
        required: true
        content:
          application/json:
            schema: SubstructureSchema
      responses:
        '200':
          content:
            application/json:
              schema: SubstructureResultSchema
    """
    data = SubstructureSchema().load(request.json)
    query = data["query"]
    target = data["target"]
    if len(query["compounds"]) * len(target["compounds"]) > 100000:
        raise ValueError(
            "Must request fewer than 100,000 substructures per request for performance reasons."
        )
    query_mols, query_skipped = convert_compound_request(query)
    target_mols, target_skipped = convert_compound_request(target)
    substructures = {
        qn: find_substructure_matches(
            q, target_mols, substructure_args=data["substructure_args"]
        )
        for qn, q in query_mols.items()
    }
    out = {"query": [], "target": [], "match": []}
    for k, v in substructures.items():
        out["query"].extend([k] * len(v))
        out["target"].extend(v.keys())
        out["match"].extend(v.values())
    SubstructureResultSchema().validate(out)
    return out


@app.route("/fingerprints/similarity_threshold", methods=["GET", "POST"])
def similarity_matches_route():
    """Find chemical similarity matches between compounds up to a certain threshold of similarity.
    ---
    post:
      summary: Find similar chemicals
      description: Find chemically similar compounds up to a certain threshold of similarity.
        This is much faster than /fingerprints/similarity because it uses the chemfp package
        and only returns matches that pass the threshold.
      requestBody:
        required: true
        content:
          application/json:
            schema: SimilarityThresholdSchema
      responses:
        '200':
          content:
            application/json:
              schema: SimilarityResultSchema
    """
    data = SimilarityThresholdSchema().load(request.json)
    query = data["query"]
    target = data["target"]
    query_mols, query_skipped = convert_compound_request(query)
    query_arena = make_fingerprint_arena(
        query_mols,
        fingerprint_type=data["fingerprint_type"],
        fingerprint_args=data["fingerprint_args"],
    )
    target_arena = None
    if target is not None:
        target_mols, target_skipped = convert_compound_request(target)
        target_arena = make_fingerprint_arena(
            target_mols,
            fingerprint_type=data["fingerprint_type"],
            fingerprint_args=data["fingerprint_args"],
        )
    matches = find_similarity_matches(
        query_arena, target_arena, threshold=data["threshold"]
    )
    out = {"query": [], "target": [], "score": []}
    for k, v in matches.items():
        out["query"].extend([k] * len(v))
        out["target"].extend(v.keys())
        out["score"].extend(v.values())
    SimilarityResultSchema().validate(out)
    return out
