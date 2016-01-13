from kidney_solver.kidney_digraph import *
import kidney_solver.kidney_ip as k_ip
import kidney_solver.kidney_ndds as k_ndds

def read_with_ndds(basename):
    with open(basename + ".input") as f:
        lines = f.readlines()
    d = read_digraph_without_prob(lines)
    with open(basename + ".ndds") as f:
        lines = f.readlines()
    ndds = k_ndds.read_ndds(lines, d)
    return d, ndds

def test_chains_only_instance():
    d, ndds = read_with_ndds("test-fixtures/no-cycles")
    d.create_adj_mat()
    fns = [k_ip.optimise_uuef, k_ip.optimise_hpief_prime,
            k_ip.optimise_picef, k_ip.optimise_ccf]
    for max_chain in [0, 1, 2]:
        for fn in fns:
            opt_result = fn(d, ndds, 3, max_chain, None)
            assert len(opt_result.cycles) == 0
            if fn == k_ip.optimise_uuef or max_chain > 0:
                assert len(opt_result.chains) == 2
            else:
                assert len(opt_result.chains) == 0

def test_one_cycle_instance():
    d, ndds = read_with_ndds("test-fixtures/one-cycle")
    d.create_adj_mat()
    fns = [k_ip.optimise_uuef, k_ip.optimise_hpief_prime,
            k_ip.optimise_picef, k_ip.optimise_ccf]
    for max_chain in [0, 1, 2]:
        for fn in fns:
            opt_result = fn(d, ndds, 3, max_chain, None)
            assert len(opt_result.cycles) == 1
            if fn == k_ip.optimise_uuef or max_chain > 0:
                assert len(opt_result.chains) == 1
            else:
                assert len(opt_result.chains) == 0
