import json
import marshmallow
from apispec import APISpec
from apispec.ext.marshmallow import MarshmallowPlugin
from apispec_webframeworks.flask import FlaskPlugin

from lspcheminf import app
import lspcheminf.schemas as schemas

spec = APISpec(
    title="TAS Chemoinformatics tools",
    version="0.4.0",
    openapi_version="3.0.2",
    plugins=[FlaskPlugin(), MarshmallowPlugin()],
)

for i in dir(schemas):
    c = getattr(schemas, i)
    try:
        if issubclass(c, marshmallow.Schema):
            spec.definition(i, c)
    except Exception:
        pass

with app.test_request_context():
    for p in app.view_functions.values():
        spec.path(view=p)

if __name__ == "__main__":
    with open("./lsp_cheminf_server/lspcheminf/static/apidoc.json", "w") as f:
        json.dump(spec.to_dict(), f)
    with open("./docs/api/apidoc.json", "w") as f:
        json.dump(spec.to_dict(), f)
