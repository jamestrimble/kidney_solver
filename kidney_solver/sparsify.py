"""
A command-line tool to sparsify an instance.
"""

import argparse
import random
import sys

from . import kidney_digraph
from . import kidney_ndds

def write_edges(n, edges):
    """Writes an instance in .input or .ndds format

    Args:
        n: Number of donor-patient pairs (for .input) or NDDs (for .ndds)
        edges: The edges in (src_as_int, tgt_as_int, weight_as_float) format
    """

    print(("{}\t{}".format(n, len(edges))))
    for edge in edges:
        print(("{}\t{}\t{}".format(edge[0], edge[1], edge[2])))

    print("-1\t-1\t-1")

if __name__=="__main__":
    parser = argparse.ArgumentParser(
            "Sparsify a kidney-exchange instance, by keeping each edge with probability p.")
    parser.add_argument("p", type=float,
            help="the probability with which to keep each edge")
            
    args = parser.parse_args()

    if not 0 < args.p < 1:
        raise ValueError("p should be in the interval (0, 1)")

    input_lines = [line for line in sys.stdin]
    n_digraph_edges = int(input_lines[0].split()[1])
    digraph_lines = input_lines[:n_digraph_edges + 2]

    d = kidney_digraph.read_digraph(digraph_lines)

    if len(input_lines) > len(digraph_lines):
        ndd_lines = input_lines[n_digraph_edges + 2:]
        ndds = kidney_ndds.read_ndds(ndd_lines, d)
    else:
        ndds = []

    pair_pair_edges = []
    for e in d.es:
        if random.random() < args.p:
            pair_pair_edges.append((e.src.id, e.tgt.id, e.score))
    write_edges(d.n, pair_pair_edges)

    ndd_pair_edges = []
    for i, ndd in enumerate(ndds):
        for e in ndd.edges:
            if random.random() < args.p:
                ndd_pair_edges.append((i, e.target_v.id, e.score))
    write_edges(len(ndds), ndd_pair_edges)
