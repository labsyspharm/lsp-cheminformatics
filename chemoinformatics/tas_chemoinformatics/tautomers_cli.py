from __future__ import print_function, absolute_import
from . import tautomers
import click
import pandas as pd


@click.group()
def cli():
    pass


@cli.command()
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
def canonicalize(compound_file, output_file, compound_col, compound_encoding):
    input_df = pd.read_csv(compound_file)
    click.echo("Read input")
    canonical, skipped = tautomers.canonicalize(
        input_df[compound_col], compound_encoding
    )
    click.echo("Finished canonicalization")
    out_df = pd.DataFrame(canonical, columns=["query", "smiles", "inchi"])
    out_df = out_df.append(pd.DataFrame({"query": skipped}), sort=True)
    out_df.to_csv(output_file, index=False)
