import argparse
import pickle
import pathlib
import sys
import tempfile
import rdkit
import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi
from rdkit.Chem.Fingerprints import FingerprintMols

identifier_mol_mapping = {"smiles": Chem.MolFromSmiles, "inchi": inchi.MolFromInchi}


def main(argvs):
    parser = argparse.ArgumentParser(
        "Write molecular fingerprints to pickle file",
        description="Calculates topological fingerprints using rdkit and default values. "
        "See https://www.rdkit.org/docs/GettingStartedInPython.html#topological-fingerprints",
    )
    parser.add_argument(
        "input_file",
        type=argparse.FileType(mode="r"),
        help="Input csv file of compounds including a chemical string column "
        "(Inchi, SMILES) and a compound identifier column.",
    )
    parser.add_argument(
        "output_file",
        type=argparse.FileType(mode="wb"),
        help="Output pickle file containg fingerprints.",
    )
    parser.add_argument(
        "--compound-col",
        default="inchi",
        choices=["inchi", "smiles"],
        help="Column containing chemical identifier in input file.",
    )
    parser.add_argument(
        "--name-col",
        default="chembl_id",
        help="Column containing compound name or database identifier.",
    )
    args = parser.parse_args(argvs)
    input_df = pd.read_csv(args.input_file.name)
    print(input_df.head())
    mol_func = identifier_mol_mapping[args.compound_col]
    mols = [mol_func(m) for m in input_df[args.compound_col]]
    print(f"Generated {len(mols)} molecules.")
    fps = FingerprintMols.FingerprintsFromMols(zip(input_df[args.name_col], mols))
    pickle.dump(fps, args.output_file)


if __name__ == "__main__":
    main(sys.argv[1:])
