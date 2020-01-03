import json
import os
import tempfile

import pytest

import lspcheminf
from test_common import test_compounds, client


def test_mass(client):
    res = client.post(
        "/properties/mass",
        json={
            "compounds": {"identifier": "inchi", "compounds": test_compounds["inchi"]}
        },
    )
    assert res.status_code == 200
    res_json = res.get_json()
    assert len(res_json["mass"]) == 4
    assert all(x > 200 for x in res_json["mass"].values())


def test_id_conversion(client):
    res = client.post(
        "/properties/convert",
        json={
            "compounds": {"identifier": "inchi", "compounds": test_compounds["inchi"]},
            "target_identifier": "smiles",
        },
    )
    assert res.status_code == 200
    res_json = res.get_json()
    assert len(res_json["compounds"]["compounds"]) == 4
    assert all(len(x) > 10 for x in res_json["compounds"]["compounds"])
