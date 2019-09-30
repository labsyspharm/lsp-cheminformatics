import json
import os
import tempfile

import pytest

import tas_chemoinformatics

@pytest.fixture
def client():
    tas_chemoinformatics.app.testing = True
    with tas_chemoinformatics.app.test_client() as client:
        yield client

def test_tautomerize(client):
    res = client.post(
        "/tautomers/enumerate",
        json={
            "compounds": {
                "identifier": "inchi",
                "compounds": [
                    "InChI=1S/C23H31N5O4/c1-15(2)14-28-21(29)19-20(25(3)23(28)30)24-22-26(10-6-11-27(19)22)12-9-16-7-8-17(31-4)18(13-16)32-5/h7-8,13,15H,6,9-12,14H2,1-5H3",
                    "InChI=1S/C31H30F3N3O3S/c1-19(2)26-13-12-25(16-20(26)3)37(17-21-4-6-23(7-5-21)29(40)35-15-14-28(38)39)30-36-27(18-41-30)22-8-10-24(11-9-22)31(32,33)34/h4-13,16,18-19H,14-15,17H2,1-3H3,(H,35,40)(H,38,39)",
                ],
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
                "compounds": [
                    "InChI=1S/C23H31N5O4/c1-15(2)14-28-21(29)19-20(25(3)23(28)30)24-22-26(10-6-11-27(19)22)12-9-16-7-8-17(31-4)18(13-16)32-5/h7-8,13,15H,6,9-12,14H2,1-5H3",
                    "InChI=1S/C31H30F3N3O3S/c1-19(2)26-13-12-25(16-20(26)3)37(17-21-4-6-23(7-5-21)29(40)35-15-14-28(38)39)30-36-27(18-41-30)22-8-10-24(11-9-22)31(32,33)34/h4-13,16,18-19H,14-15,17H2,1-3H3,(H,35,40)(H,38,39)",
                ],
            },
            "standardize": True
        },
    )
    assert res.status_code == 200
    res_json = res.get_json()
    assert len(res_json["canonical"]) == 2
    assert len(res_json["skipped"]) == 0
