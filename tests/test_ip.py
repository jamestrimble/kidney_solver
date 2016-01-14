from kidney_solver.kidney_digraph import *
import kidney_solver.kidney_ip as k_ip
import kidney_solver.kidney_ndds as k_ndds

def read_with_ndds(basename):
    with open(basename + ".input") as f:
        lines = f.readlines()
    d = read_digraph(lines)
    with open(basename + ".ndds") as f:
        lines = f.readlines()
    ndds = k_ndds.read_ndds(lines, d)
    return d, ndds

def test_chains_only_instance():
    d, ndds = read_with_ndds("test-fixtures/no-cycles")
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

def test_single_cycle_instance():
    d, ndds = read_with_ndds("test-fixtures/one-cycle")
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

def test_weighted_instance():
    """Checks that the capped formulations agree on the optimal
    result for an instance with weighted edges.
    """
    d, ndds = read_with_ndds("test-fixtures/100-random-weights")
    fns = [k_ip.optimise_hpief_prime,
            k_ip.optimise_picef, k_ip.optimise_ccf]
    for max_cycle in [0, 1, 2, 3, 4]:
        for max_chain in [0, 1, 2, 3]:
            opt_result_0 = fns[0](d, ndds, max_cycle, max_chain, None)
            print opt_result_0.total_score
            for fn in fns[1:]:
                opt_result = fn(d, ndds, max_cycle, max_chain, None)
                assert opt_result.total_score == opt_result_0.total_score

def test_relabelled_solve():
    """Checks whether the vertex-ordering heuristic affects the
    result for a weighted instance
    """
    EPS = 0.000001
    d, ndds = read_with_ndds("test-fixtures/100-random-weights")
    for max_cycle in [0, 3]:
        for max_chain in [0, 5]:
            opt_result_0 = k_ip.optimise_picef(d, ndds, max_cycle, max_chain, None)
            print opt_result_0.total_score
            opt_result = k_ip.optimise_relabelled(
                    k_ip.optimise_picef, d, ndds, max_cycle, max_chain, None)
            print "   ", opt_result.total_score
            assert abs(opt_result.total_score - opt_result_0.total_score) < EPS
