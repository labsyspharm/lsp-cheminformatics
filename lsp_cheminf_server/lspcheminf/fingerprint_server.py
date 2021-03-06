import base64
import gzip
import os
import pickle
import sys
import tempfile
import pandas as pd
from flask import Flask, request, send_file
from rdkit.Chem import AllChem, inchi, rdFMCS
from rdkit.Chem.Scaffolds import MurckoScaffold

from lspcheminf import app
from lspcheminf.util import (
    convert_compound_request,
    mol_identifier_mapping,
    identifier_mol_mapping,
)
from lspcheminf.fingerprint import (
    calculate_similarity,
    find_substructure_matches,
    find_similarity_matches,
    make_fingerprint_arena,
    compound_identity,
    calculate_fingerprints,
)
from lspcheminf.schemas import (
    SimilaritySchema,
    SimilarityResultSchema,
    SimilarityThresholdSchema,
    SubstructureSchema,
    SubstructureResultSchema,
    CompoundIdentitySchema,
    CompoundIdentityResultSchema,
    CalculateFingerprintsSchema,
    CalculateFingerprintsResultSchema,
    MCSSchema,
    MCSResultSchema,
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
    target = data.get("target", None)
    query_mols, query_skipped = convert_compound_request(
        query, field="fingerprints" if "fingerprints" in query else "compounds"
    )
    query_arena = make_fingerprint_arena(
        query_mols,
        fingerprint_type=data["fingerprint_type"],
        fingerprint_args=data["fingerprint_args"],
    )
    del query_mols
    target_arena = None
    if target is not None:
        target_mols, target_skipped = convert_compound_request(
            target, field="fingerprints" if "fingerprints" in target else "compounds"
        )
        target_arena = make_fingerprint_arena(
            target_mols,
            fingerprint_type=data["fingerprint_type"],
            fingerprint_args=data["fingerprint_args"],
        )
        del target_mols
    matches = find_similarity_matches(
        query_arena,
        target_arena,
        threshold=data["threshold"],
        n_threads=data["n_threads"],
    )
    out = {"query": [], "target": [], "score": []}
    for k, v in matches.items():
        out["query"].extend([k] * len(v))
        out["target"].extend(v.keys())
        out["score"].extend(v.values())
    SimilarityResultSchema().validate(out)
    return out


@app.route("/fingerprints/compound_identity", methods=["GET", "POST"])
def compound_identity_route():
    """Establish identity between compounds using both morgan and topological fingerprints and molecular weight.
    ---
    post:
      summary: Find identical compounds
      description: Compare fingerprints between query and target compounds and
        find identical compounds using both morgan and topological fingerprints as well as
        molecular weight.
      requestBody:
        required: true
        content:
          application/json:
            schema: CompoundIdentitySchema
      responses:
        '200':
          content:
            application/json:
              schema: CompoundIdentityResultSchema
    """
    data = CompoundIdentitySchema().load(request.json)
    query = data["query"]
    target = data.get("target", None)
    query_mols, query_skipped = convert_compound_request(query)
    target_mols = None
    if target is not None:
        target_mols, target_skipped = convert_compound_request(target)
    matches = compound_identity(query_mols, target_mols)
    out = {"query": [], "target": []}
    for k, v in matches.items():
        out["query"].extend([k] * len(v))
        out["target"].extend(v)
    CompoundIdentityResultSchema().validate(out)
    return out


@app.route("/fingerprints/calculate", methods=["GET", "POST"])
def calculate_fingerprint_route():
    """Calculate chemical fingerprints.
    ---
    post:
      summary: Calculate chemical fingerprints
      description: Calculate chemical fingerprints for the supplied compound identifiers.
      requestBody:
        required: true
        content:
          application/json:
            schema: CalculateFingerprintsSchema
      responses:
        '200':
          content:
            application/json:
              schema: CalculateFingerprintsResultSchema
    """
    data = CalculateFingerprintsSchema().load(request.json)
    query = data["query"]
    if len(query["compounds"]) > 100000:
        raise ValueError(
            "Must request fewer than 100,000 compounds per request for performance reasons."
        )
    query_mols, query_skipped = convert_compound_request(query)
    fingerprints = calculate_fingerprints(
        query_mols,
        fingerprint_type=data["fingerprint_type"],
        fingerprint_args=data["fingerprint_args"],
    )
    out = {
        "names": list(fingerprints.keys()),
        "fingerprints": list(fingerprints.values()),
    }
    CalculateFingerprintsResultSchema().validate(out)
    return out


@app.route("/fingerprints/maximum_common_substructure", methods=["GET", "POST"])
def maximum_common_substructure_route():
    """Calculate maximum common substructure (MCS) of compounds.
    ---
    post:
      summary: Calculate MCS
      description: Calculate maximum common substructure (MCS) of compounds.
      requestBody:
        required: true
        content:
          application/json:
            schema: MCSSchema
      responses:
        '200':
          content:
            application/json:
              schema: MCSResultSchema
    """
    data = MCSSchema().load(request.json)
    query = data["query"]
    query_mols, query_skipped = convert_compound_request(query)
    substructure = rdFMCS.FindMCS(list(query_mols.values()))
    substructure_mol = identifier_mol_mapping["smarts"](substructure.smartsString)
    out = {
        "substructure": {
            "compounds": mol_identifier_mapping[query["identifier"]](substructure_mol),
            "identifier": query["identifier"],
        },
        "skipped": query_skipped,
    }
    MCSResultSchema().validate(out)
    return out


@app.route("/fingerprints/murcko_scaffold", methods=["GET", "POST"])
def murcko_scaffold_route():
    """Calculate Murcko Scaffold of compounds.
    ---
    post:
      summary: Calculate Murcko scaffold
      description: Calculate Murcko scaffold of compounds.
      requestBody:
        required: true
        content:
          application/json:
            schema: MCSSchema
      responses:
        '200':
          content:
            application/json:
              schema: MCSResultSchema
    """
    data = MCSSchema().load(request.json)
    query = data["query"]
    query_mols, query_skipped = convert_compound_request(query)
    scaffolds = {
        n: MurckoScaffold.GetScaffoldForMol(mol) for n, mol in query_mols.items()
    }
    out = {
        "scaffolds": {
            "compounds": [
                mol_identifier_mapping[query["identifier"]](mol)
                for mol in scaffolds.values()
            ],
            "names": list(scaffolds.keys()),
            "identifier": query["identifier"],
        },
        "skipped": query_skipped,
    }
    MCSResultSchema().validate(out)
    return out
