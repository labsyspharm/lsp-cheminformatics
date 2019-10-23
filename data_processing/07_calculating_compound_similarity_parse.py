#!/usr/bin/env python3
import logging
import sys

import click

# log = logging.getLogger(__name__)

def pairwise(iterable):
    "s -> (s0, s1), (s2, s3), (s4, s5), ..."
    a = iter(iterable)
    return zip(a, a)


@click.command()
@click.argument("input", type=click.File("r"))
@click.argument("output", type=click.File("w"))
def parse_sim_result_cli(input, output):
    parse_sim_result(input, output)

def parse_sim_result(sim_result_it, out_handle):
    # with open(sim_result_file, "r") as f, open(out_file, "w") as o:
    for i, line in enumerate(sim_result_it):
        if i % 1000 == 0:
            print("Processed", i, file=sys.stderr, flush=True)
        if line.startswith("#"):
            continue
        l = line.split()
        if len(l) <= 2:
            continue
        query = l[1]
        for match, score in pairwise(l[2:]):
            o_ids = (query, match) if query < match else (match, query)
            print(o_ids[0], o_ids[1], score, sep="\t", file=out_handle)
    print("Done", file=sys.stderr, flush=True)


if __name__ == "__main__":
    parse_sim_result_cli()
