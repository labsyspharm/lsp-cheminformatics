from marshmallow import Schema, fields, validate

from .util import identifier_mol_mapping


class CompoundsSchema(Schema):
    compounds = fields.List(fields.String, required=True)
    identifier = fields.String(
        validate=validate.OneOf(identifier_mol_mapping.keys()), required=True
    )

class TautomerizeSchema(Schema):
    compounds = fields.Nested(CompoundsSchema, required=True)
    max_tautomers = fields.Integer(default=10)

class TautomerizeResultSchema(Schema):
    request = fields.Nested(TautomerizeSchema)
    tautomers = fields.Mapping(keys=fields.String, values=fields.Nested(CompoundsSchema), required=True)

class CanonicalizeSchema(Schema):
    compounds = fields.Nested(CompoundsSchema, required=True)
    standardize = fields.Boolean(default=True)

class CanonicalizeResultSchema(Schema):
    request = fields.Nested(CanonicalizeSchema)
    canonical = fields.Mapping(keys=fields.String, values=fields.Nested(CompoundsSchema), required=True)
    skipped = fields.List(fields.String)
