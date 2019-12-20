from __future__ import absolute_import, print_function

import click
import pandas as pd

from . import tautomers


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
def canonicalize(
    compound_file, output_file, compound_col, compound_encoding, standardize
):
    input_df = pd.read_csv(compound_file)
    click.echo("Read input")
    canonical, skipped = tautomers.canonicalize(
        input_df[compound_col], compound_encoding, standardize=standardize
    )
    click.echo("Finished canonicalization")
    out_df = pd.DataFrame(canonical, columns=["query", "smiles", "inchi"])
    out_df = out_df.append(pd.DataFrame({"query": skipped}), sort=True)
    out_df.to_csv(output_file, index=False)
