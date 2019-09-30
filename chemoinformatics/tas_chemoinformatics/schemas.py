from marshmallow import Schema, fields, validate

from .util import identifier_mol_mapping
from .fingerprint_server import fp_types


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
    identifier = fields.String(
        validate=validate.OneOf(identifier_mol_mapping.keys()),
        required=True,
        description="The type of compound identifiers used",
        example="inchi",
    )


class TautomerizeSchema(Schema):
    compounds = fields.Nested(CompoundsSchema, required=True)
    max_tautomers = fields.Integer(
        missing=10, description="Maximum number of tautomers to generate"
    )


class TautomerizeResultSchema(Schema):
    request = fields.Nested(TautomerizeSchema)
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
    request = fields.Nested(CanonicalizeSchema)
    canonical = fields.Mapping(
        keys=fields.String, values=fields.Nested(CompoundsSchema), required=True
    )
    skipped = fields.List(fields.String)


class FingerprintDBSchema(Schema):
    compounds = fields.Nested(CompoundsSchema, required=True)
    fingerprint_type = fields.String(
        validate=validate.OneOf(fp_types.keys()),
        missing="topological",
        description="The type of fingerprinting algorithm to be used",
    )
    fingerprint_args = fields.List(
        fields.String,
        missing=[],
        description="Optional additional arguments passed to chemfp. "
        "See https://chemfp.readthedocs.io/en/chemfp-1.5/tool-help.html#rdkit2fps-command-line-options",
        example="['--useChirality', '0']",
    )


class FingerprintDBResultSchema(Schema):
    request = fields.Nested(FingerprintDBSchema)
    fingerprint_db = fields.String(
        required=True,
        description="Fingerprint database file in .fps format as base64 encoded string",
    )
    skipped = fields.List(fields.String)


class FingerprintScanSchema(Schema):
    fingerprint_db = fields.String(
        required=True,
        description="Reference database that is searched for fingerprints matching the query database. "
        "Fingerprint database file in .fps format as base64 encoded string.",
    )
    threshold = fields.Float(
        missing=0.95,
        description="Threshold for reporting matches calculated by Tanimoto similarity",
    )
    query = fields.String(
        description="Query database with fingerprints to be matched to the reference database. "
        "If ommitted, an  all-by-all search using  all  reference fingerprints is performed. "
        "Fingerprint database file in .fps format as base64 encoded string."
    )


class FingerprintScanResultSchema(Schema):
    request = fields.Nested(FingerprintScanSchema)
    query = fields.List(fields.String, required=True)
    match = fields.List(fields.String, required=True)
    score = fields.List(fields.Float, required=True)
