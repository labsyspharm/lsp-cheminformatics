import json
import os
import tempfile

import pytest

import lspcheminf
from test_common import test_compounds, client


@pytest.mark.parametrize(
    "fingerprint_type,fingerprint_args",
    [("morgan", {}), ("morgan", {"radius": 3}), ("topological", {})],
)
def test_similarity(client, fingerprint_type, fingerprint_args):
    res = client.post(
        "/fingerprints/similarity",
        json={
            "query": {"identifier": "inchi", "compounds": test_compounds["inchi"][:2]},
            "target": {
                "identifier": "smiles",
                "compounds": test_compounds["smiles"][1:4],
            },
            "fingerprint_type": fingerprint_type,
            "fingerprint_args": fingerprint_args,
        },
    )
    assert res.status_code == 200
    res_json = res.get_json()
    assert all(len(x) == 6 for x in res_json.values())


@pytest.mark.parametrize(
    "substructure_args", [{}, {"useChirality": True, "maxMatches": 3}]
)
def test_substructure(client, substructure_args):
    res = client.post(
        "/fingerprints/substructure",
        json={
            "query": {"identifier": "inchi", "compounds": test_compounds["inchi"][:2]},
            "target": {
                "identifier": "smiles",
                "compounds": test_compounds["smiles"][1:4],
            },
            "substructure_args": substructure_args,
        },
    )
    assert res.status_code == 200
    res_json = res.get_json()
    assert all(len(x) == 1 for x in res_json.values())
    assert len(res_json["match"][0][0]) == 41


@pytest.mark.parametrize(
    "fingerprint_type,fingerprint_args",
    [("morgan", {}), ("morgan", {"radius": 3}), ("topological", {})],
)
def test_similarity_threshold(client, fingerprint_type, fingerprint_args):
    res = client.post(
        "/fingerprints/similarity_threshold",
        json={
            "query": {"identifier": "inchi", "compounds": test_compounds["inchi"][:2]},
            "target": {
                "identifier": "smiles",
                "compounds": test_compounds["smiles"][1:4],
            },
            "fingerprint_type": fingerprint_type,
            "fingerprint_args": fingerprint_args,
            "threshold": 0.8,
        },
    )
    assert res.status_code == 200
    res_json = res.get_json()
    assert all(len(x) == 1 for x in res_json.values())
    assert res_json["query"][0] == "1"
    assert res_json["target"][0] == "0"


def test_similarity_threshold_fp(client):
    res = client.post(
        "/fingerprints/similarity_threshold",
        json={
            "query": {"identifier": "inchi", "compounds": ["a", "b", "c"], "fingerprints": ["880DF", "881DF", "870E0"]},
            "threshold": 0.8,
        },
    )
    assert res.status_code == 200
    res_json = res.get_json()
    assert all(len(x) == 1 for x in res_json.values())
    assert res_json["query"][0] == "0"
    assert res_json["target"][0] == "1"


def test_compound_identity(client):
    res = client.post(
        "/fingerprints/compound_identity",
        json={
            "query": {"identifier": "inchi", "compounds": test_compounds["inchi"][:2]},
            "target": {
                "identifier": "smiles",
                "compounds": test_compounds["smiles"][1:4],
            },
        },
    )
    assert res.status_code == 200
    res_json = res.get_json()
    assert all(len(x) == 1 for x in res_json.values())
    assert res_json["query"][0] == "1"
    assert res_json["target"][0] == "0"


@pytest.mark.parametrize(
    "fingerprint_type,fingerprint_args",
    [("morgan", {}), ("morgan", {"radius": 3}), ("topological", {})],
)
def test_calculate_fingerprints(client, fingerprint_type, fingerprint_args):
    res = client.post(
        "/fingerprints/calculate",
        json={
            "query": {"identifier": "inchi", "compounds": test_compounds["inchi"][:2]},
            "fingerprint_type": fingerprint_type,
            "fingerprint_args": fingerprint_args,
        },
    )
    assert res.status_code == 200
    res_json = res.get_json()
    assert all(len(x) == 2 for x in res_json.values())
