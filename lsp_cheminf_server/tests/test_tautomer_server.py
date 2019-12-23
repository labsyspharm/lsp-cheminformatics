import json
import os
import tempfile

import pytest

import lspcheminf

from test_common import test_compounds, client


def test_tautomerize(client):
    res = client.post(
        "/tautomers/enumerate",
        json={
            "compounds": {
                "identifier": "inchi",
                "compounds": test_compounds["inchi"][:2],
            },
            "max_tautomers": 5,
        },
    )
    assert res.status_code == 200
    res_json = res.get_json()
    assert len(res_json["tautomers"]) == 2
    assert tuple(len(v) for v in res_json["tautomers"].values()) == (1, 6)


def test_canonicalize(client):
    res = client.post(
        "/tautomers/canonicalize",
        json={
            "compounds": {
                "identifier": "inchi",
                "compounds": test_compounds["inchi"][:2],
            },
            "standardize": True,
        },
    )
    assert res.status_code == 200
    res_json = res.get_json()
    assert len(res_json["canonical"]) == 2
    assert len(res_json["skipped"]) == 0
