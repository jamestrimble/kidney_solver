from kidney_solver.kidney_digraph import *
import kidney_solver.kidney_ip as k_ip
import kidney_solver.kidney_ndds as k_ndds
import kidney_solver.kidney_utils as k_utils
import kidney_solver.count_cycles_and_chains as ccc

def read_with_ndds(basename):
    with open(basename + ".input") as f:
        lines = f.readlines()
    d = read_digraph(lines)
    with open(basename + ".ndds") as f:
        lines = f.readlines()
    ndds = k_ndds.read_ndds(lines, d)
    return d, ndds

def test_counter():
    d, ndds = read_with_ndds("test-fixtures/100-random-weights")
    max_length = 4
    cycles = d.find_cycles(max_length)
    cycle_counts = [0] * (max_length + 1)
    for c in cycles:
        cycle_counts[len(c)] += 1
    assert cycle_counts == ccc.count_cycles(d, max_length)

    chains = k_ndds.find_chains(d, ndds, max_length)
    chain_counts = [0] * (max_length + 1)
    for c in chains:
        chain_counts[len(c.vtx_indices)] += 1
    assert chain_counts == ccc.count_chains(d, ndds, max_length)
    
    
