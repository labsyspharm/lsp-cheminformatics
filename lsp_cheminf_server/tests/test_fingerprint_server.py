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
            "query": {
                "identifier": "inchi",
                "compounds": ["a", "b", "c"],
                "fingerprints": ["880DF", "881DF", "870E0"],
            },
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


def test_maximum_common_substructure(client):
    res = client.post(
        "/fingerprints/maximum_common_substructure",
        json={
            "query": {
                "identifier": "smiles",
                "compounds": [
                    "O=C(NCc1cc(OC)c(O)cc1)CCCC/C=C/C(C)C",
                    "CC(C)CCCCCC(=O)NCC1=CC(=C(C=C1)O)OC",
                    "c1(C=O)cc(OC)c(O)cc1",
                ],
            },
        },
    )
    assert res.status_code == 200
    res_json = res.get_json()
    assert res_json["substructure"]["compounds"] == "COC1:C:C(C):C:C:C:1O"


def test_murcko_scaffold(client):
    res = client.post(
        "/fingerprints/murcko_scaffold",
        json={
            "query": {
                "identifier": "smiles",
                "compounds": [
                    "O=C(NCc1cc(OC)c(O)cc1)CCCC/C=C/C(C)C",
                    "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5",
                ],
            },
        },
    )
    assert res.status_code == 200
    res_json = res.get_json()
    assert res_json["scaffolds"]["compounds"][0] == "c1ccccc1"
