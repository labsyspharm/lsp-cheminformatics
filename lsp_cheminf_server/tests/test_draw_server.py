import json
import os
import tempfile

import pytest

import lspcheminf
from test_common import test_compounds, client


def test_draw_grid(client):
    res = client.post(
        "/draw/grid",
        json={
            "compounds": {"identifier": "inchi", "compounds": test_compounds["inchi"]}
        },
    )
    assert res.status_code == 200
    res_json = res.get_json()
    assert len(res_json["svg"]) > 50000
