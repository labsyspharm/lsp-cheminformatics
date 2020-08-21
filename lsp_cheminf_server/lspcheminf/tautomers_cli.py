from __future__ import absolute_import, print_function

import click
import pandas as pd
from progress.bar import Bar

from . import tautomers
from .util import convert_compound_request, mol_identifier_mapping


@click.group()
def cli():
    pass


@cli.command(
    short_help="Find canonical tautomer for input compounds.",
    help="Find canonical tautomer for all compounds in the `COMPOUND_FILE` input. "
    "This file should be a csv file where one column contains compounds as InChI or SMILES."
    "The canonical compounds will be written as InChI and SMILES to the `OUTPUT_FILE`.",
)
@click.argument(
    "compound_file",
    type=click.File("r"),
    # help="CSV file with column containing compound string",
)
@click.argument(
    "output_file",
    type=click.File("w"),
    # help="Output file where found canonical tautomers are written",
)
@click.option(
    "--compound-col",
    default="compound",
    show_default=True,
    help="Name of column in compound_file in which compound strings are stored",
)
@click.option(
    "--compound-encoding",
    default="inchi",
    show_default=True,
    type=click.Choice(["inchi", "smiles"]),
    help="Encoding format of compound strings",
)
@click.option(
    "--standardize/--no-standardize",
    default=True,
    show_default=True,
    help="Standardize input molecule before canonicalization: https://molvs.readthedocs.io/en/latest/guide/standardize.html",
)
@click.option(
    "--standardizer",
    default="chembl",
    show_default=True,
    type=click.Choice(tautomers.STANDARDIZERS.keys()),
    help="Standardizer used. Either molvs, chembl or chembl-parent (https://github.com/chembl/ChEMBL_Structure_Pipeline)",
)
def canonicalize(
    compound_file,
    output_file,
    compound_col,
    compound_encoding,
    standardize,
    standardizer,
):
    input_df = pd.read_csv(compound_file)
    click.echo("Read input")
    compounds = convert_compound_request(
        {"compounds": input_df[compound_col], "identifier": compound_encoding}
    )
    progress = Bar(max=len(input_df.index))
    canonical, skipped = tautomers.canonicalize(
        compounds[0],
        standardize=standardize,
        standardizer=standardizer,
        progress_callback=progress.goto,
    )
    progress.finish()
    click.echo("Finished canonicalization")
    mol_to_inchi = mol_identifier_mapping["inchi"]
    out_df = pd.DataFrame(
        {
            "row": list(canonical.keys()),
            "inchi": list(mol_to_inchi(x) for x in canonical.values()),
        }
    )
    out_df = out_df.append(pd.DataFrame({"row": skipped}), sort=True)
    out_df.to_csv(output_file, index=False)
