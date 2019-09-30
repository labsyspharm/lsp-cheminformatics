import json
import marshmallow
from apispec import APISpec
from apispec.ext.marshmallow import MarshmallowPlugin
from apispec_webframeworks.flask import FlaskPlugin

from tas_chemoinformatics import app
import tas_chemoinformatics.schemas as schemas

spec = APISpec(
    title="TAS Chemoinformatics tools",
    version="0.2",
    openapi_version="3.0.2",
    plugins=[FlaskPlugin(), MarshmallowPlugin()]
)

for i in dir(schemas):
    c = getattr(schemas, i)
    try:
        if isinstance(c, marshmallow.Schema):
            spec.definition(i, c)
    except Exception:
        pass

with app.test_request_context():
    for p in app.view_functions.values():
        spec.path(view=p)

if __name__ == "__main__":
    with open("./tas_chemoinformatics/static/apidoc.json", "w") as f:
        json.dump(spec.to_dict(), f)
