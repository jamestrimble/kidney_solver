"""
A program for counting the number of cycles and chains in a kidney-exchange
instance, up to given cycle and chain caps.
"""

import argparse
import sys

from . import kidney_digraph
from . import kidney_ndds

def count_cycles(digraph, max_length):
    """Find cycles of length up to max_length in the digraph.

    Return value: a list of cycles, in which each cycle is represented
    as a list of vertices, with the first vertex not repeated.
    """

    counts_by_size = [0] * (max_length + 1)

    if max_length < 2:
        return counts_by_size

    vtx_used = [False] * len(digraph.vs)

    def cycle(first_vtx, last_vtx, num_vertices):
        if digraph.edge_exists(last_vtx, first_vtx):
            counts_by_size[num_vertices] += 1
        if num_vertices < max_length:
            for e in last_vtx.edges: 
                v = e.tgt
                if (num_vertices + shortest_paths_to_low_vtx[v.id] <= max_length
                            and not vtx_used[v.id]):
                    vtx_used[v.id] = True
                    cycle(first_vtx, v, num_vertices + 1)
                    vtx_used[v.id] = False

    # Adjacency lists for transpose graph
    transp_adj_lists = [[] for v in digraph.vs]
    for edge in digraph.es:
        transp_adj_lists[edge.tgt.id].append(edge.src)

    for v in digraph.vs:
        shortest_paths_to_low_vtx = digraph.calculate_shortest_path_lengths(
                v, max_length - 1,
                lambda u: (w for w in transp_adj_lists[u.id] if w.id > v.id))
        vtx_used[v.id] = True
        cycle(v, v, 1)
        vtx_used[v.id] = False
    return counts_by_size

def count_chains(digraph, ndds, max_chain):
    vtx_used = [False] * len(digraph.vs)
    def find_chains_recurse(last_vtx, length):
        counts_by_size[length] += 1
        if length < max_chain:
            for e in digraph.vs[last_vtx].edges:
                if not vtx_used[e.tgt.id]:
                    vtx_used[e.tgt.id] = True
                    find_chains_recurse(e.tgt.id, length+1)
                    vtx_used[e.tgt.id] = False
    counts_by_size = [0] * (max_chain + 1)

    if max_chain == 0:
        return counts_by_size
    for ndd_idx, ndd in enumerate(ndds):
        for e in ndd.edges:
            vtx_used[e.target_v.id] = True
            find_chains_recurse(e.target_v.id, 1)
            vtx_used[e.target_v.id] = False
    return counts_by_size


def start():
    parser = argparse.ArgumentParser(
            "Count cycles and chains in a kidney-exchange instance")
    parser.add_argument("cycle_cap", type=int,
            help="The maximum permitted cycle length")
    parser.add_argument("chain_cap", type=int,
            help="The maximum permitted number of edges in a chain")
            
    args = parser.parse_args()

    input_lines = [line for line in sys.stdin]
    n_digraph_edges = int(input_lines[0].split()[1])
    digraph_lines = input_lines[:n_digraph_edges + 2]

    d = kidney_digraph.read_digraph(digraph_lines)

    if len(input_lines) > len(digraph_lines):
        ndd_lines = input_lines[n_digraph_edges + 2:]
        altruists = kidney_ndds.read_ndds(ndd_lines, d)
    else:
        altruists = []
        
    cycle_counts_by_size = count_cycles(d, args.cycle_cap)
    chain_counts_by_size = count_chains(d, altruists, args.chain_cap)
    print(("cycles: {}".format(sum(cycle_counts_by_size))))
    print(("chains: {}".format(sum(chain_counts_by_size))))
    print("cycles by length:")
    for i in range(2, args.cycle_cap+1):
        print(("{}\t{}".format(i, cycle_counts_by_size[i])))
    print("chains by length:")
    for i in range(1, args.chain_cap+1):
        print(("{}\t{}".format(i, chain_counts_by_size[i])))

if __name__=="__main__":
    start()
