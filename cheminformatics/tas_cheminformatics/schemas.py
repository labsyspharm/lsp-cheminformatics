from marshmallow import Schema, fields, validate

from .util import identifier_mol_mapping
from .fingerprint import fingerprint_functions

identifier_field = fields.String(
    validate=validate.OneOf(identifier_mol_mapping.keys()),
    required=True,
    description="The type of compound identifiers used",
    example="inchi",
)

skipped_field = fields.List(
    fields.String,
    description="A list of compound identifiers that couldn't be processed",
)


class CompoundsSchema(Schema):
    compounds = fields.List(
        fields.String,
        required=True,
        description="A list of compound identifiers",
        example="['InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)', 'InChI=1S/C8H11NO2/c9-4-3-6-1-2-7(10)8(11)5-6/h1-2,5,10-11H,3-4,9H2']",
    )
    names = fields.List(
        fields.String,
        description="An optional list of compound names",
        example="['Aspirin', 'Dopamine']",
    )
    identifier = identifier_field


class TautomerizeSchema(Schema):
    compounds = fields.Nested(CompoundsSchema, required=True)
    max_tautomers = fields.Integer(
        missing=10, description="Maximum number of tautomers to generate"
    )


class TautomerizeResultSchema(Schema):
    tautomers = fields.Mapping(
        keys=fields.String, values=fields.Nested(CompoundsSchema), required=True
    )


class CanonicalizeSchema(Schema):
    compounds = fields.Nested(CompoundsSchema, required=True)
    standardize = fields.Boolean(
        missing=True,
        description="Standardize input molecule before canonicalization. "
        "See https://molvs.readthedocs.io/en/latest/guide/standardize.html",
    )


class CanonicalizeResultSchema(Schema):
    canonical = fields.Mapping(
        keys=fields.String, values=fields.Nested(CompoundsSchema), required=True
    )
    skipped = skipped_field


class SimilaritySchema(Schema):
    query = fields.Nested(
        CompoundsSchema,
        description="Query compounds to be compared with the target compounds",
        required=True,
    )
    target = fields.Nested(
        CompoundsSchema,
        description="Target compounds to be compared with the query compounds",
        required=True,
    )
    fingerprint_type = fields.String(
        validate=validate.OneOf(fingerprint_functions.keys()),
        missing="topological",
        description="The type of fingerprinting algorithm to be used",
    )
    fingerprint_args = fields.Mapping(
        keys=fields.String,
        values=fields.Field,
        missing={},
        description="Optional additional arguments passed to RDKit fingerprinting functions. "
        "See https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html#rdkit.Chem.rdmolops.RDKFingerprint "
        "and https://www.rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html#rdkit.Chem.rdMolDescriptors.GetMorganFingerprint",
        example='{"minPath": 2, "useHs": false}',
    )


class SimilarityResultSchema(Schema):
    query = fields.List(fields.String, required=True)
    target = fields.List(fields.String, required=True)
    score = fields.List(fields.Float, required=True)


class ConvertIDSchema(Schema):
    compounds = fields.Nested(CompoundsSchema, required=True)
    target_identifier = identifier_field


class ConvertIdResultSchema(Schema):
    compounds = fields.Nested(CompoundsSchema, required=True)
    skipped = skipped_field


class DrawGridSchema(Schema):
    compounds = fields.Nested(CompoundsSchema, required=True)


class DrawGridResultSchema(Schema):
    svg = fields.String()
    skipped = skipped_field


class CalculateMassSchema(Schema):
    compounds = fields.Nested(CompoundsSchema, required=True)


class CalculateMassResultSchema(Schema):
    mass = fields.Mapping(keys=fields.String, values=fields.Float, required=True)
    skipped = skipped_field
