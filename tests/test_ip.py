from kidney_solver.kidney_digraph import *
import kidney_solver.kidney_ip as k_ip
import kidney_solver.kidney_ndds as k_ndds
import kidney_solver.kidney_utils as k_utils

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
    fns = [k_ip.optimise_uuef,
            k_ip.optimise_hpief_prime, k_ip.optimise_hpief_2prime,
            k_ip.optimise_hpief_prime_full_red, k_ip.optimise_hpief_2prime_full_red,
            k_ip.optimise_eef, k_ip.optimise_eef_full_red,
            k_ip.optimise_picef, k_ip.optimise_ccf]
    for max_chain in [0, 1, 2]:
        for fn in fns:
            opt_result = fn(k_ip.OptConfig(d, ndds, 3, max_chain, None))
            assert len(opt_result.cycles) == 0
            if fn == k_ip.optimise_uuef or max_chain > 0:
                assert len(opt_result.chains) == 2
            else:
                assert len(opt_result.chains) == 0

def test_single_cycle_instance():
    d, ndds = read_with_ndds("test-fixtures/one-cycle")
    fns = [k_ip.optimise_uuef,
            k_ip.optimise_hpief_prime,
            k_ip.optimise_hpief_2prime,
            k_ip.optimise_hpief_prime_full_red,
            k_ip.optimise_hpief_2prime_full_red,
            k_ip.optimise_eef,
            k_ip.optimise_eef_full_red,
            k_ip.optimise_picef, k_ip.optimise_ccf]
    for max_chain in [0, 1, 2]:
        for fn in fns:
            opt_result = fn(k_ip.OptConfig(d, ndds, 3, max_chain))
            assert len(opt_result.cycles) == 1
            if fn == k_ip.optimise_uuef or max_chain > 0:
                assert len(opt_result.chains) == 1
            else:
                assert len(opt_result.chains) == 0

def test_instance_with_two_2cycles():
    d, ndds = read_with_ndds("test-fixtures/two_2cycles")
    fns = [k_ip.optimise_uuef,
            k_ip.optimise_eef, k_ip.optimise_eef_full_red,
            k_ip.optimise_hpief_prime_full_red, k_ip.optimise_hpief_2prime_full_red,
            k_ip.optimise_hpief_prime, k_ip.optimise_hpief_2prime,
            k_ip.optimise_picef, k_ip.optimise_ccf]
    for fn in fns:
        opt_result = fn(k_ip.OptConfig(d, ndds, 3, 0))
        assert len(opt_result.cycles) == 1
        assert opt_result.total_score == 40

def test_preflib_instance_with_zero_chain_cap():
    d, ndds = read_with_ndds("test-fixtures/MD-00001-00000100")
    fns = [ k_ip.optimise_eef, k_ip.optimise_eef_full_red,
            k_ip.optimise_hpief_prime_full_red, k_ip.optimise_hpief_2prime_full_red,
            k_ip.optimise_hpief_prime, k_ip.optimise_hpief_2prime,
            k_ip.optimise_eef, k_ip.optimise_eef_full_red,
            k_ip.optimise_picef, k_ip.optimise_ccf]
    for fn in fns:
        max_cycle = 4
        max_chain = 0
        opt_result = fn(k_ip.OptConfig(d, ndds, max_cycle, max_chain))
        k_utils.check_validity(opt_result, d, ndds, max_cycle, max_chain)
        print(fn)
        assert opt_result.total_score == 39 

def test_failure_aware():
    # Use instance with a 3-cycle and a 3-chain
    d, ndds = read_with_ndds("test-fixtures/3cycle_and_3chain")
    fns = [k_ip.optimise_picef, k_ip.optimise_ccf]
    for fn in fns:
        for max_cycle, max_chain, expected_score in [
                (3, 2, 3*.125 + .75),
                (3, 3, 3*.125 + .875),
                (0, 3, .875),
                (3, 0, 3*.125)]:
        
            opt_result = fn(k_ip.OptConfig(d, ndds, max_cycle, max_chain, edge_success_prob=.5))
            # 3 edges * .5**3 for the cycle; .5 + .5**2 for the chain
            assert opt_result.total_score == expected_score
            k_utils.check_validity(opt_result, d, ndds, max_cycle, max_chain)

def test_weighted_instance():
    """Checks that the capped formulations agree on the optimal
    result for an instance with weighted edges.
    """
    EPS = 0.000001
    d, ndds = read_with_ndds("test-fixtures/100-random-weights")
    fns = [k_ip.optimise_hpief_prime,
            k_ip.optimise_hpief_2prime,
            k_ip.optimise_hpief_prime_full_red,
            k_ip.optimise_hpief_2prime_full_red,
            k_ip.optimise_picef, k_ip.optimise_ccf,
            k_ip.optimise_eef, k_ip.optimise_eef_full_red]
    for max_cycle in [0, 1, 2, 3, 4]:
        for max_chain in [0, 1, 2, 3]:
            opt_result_0 = fns[0](k_ip.OptConfig(d, ndds, max_cycle, max_chain))
            k_utils.check_validity(opt_result_0, d, ndds, max_cycle, max_chain)
            for fn in fns[1:]:
                opt_result = fn(k_ip.OptConfig(d, ndds, max_cycle, max_chain))
                k_utils.check_validity(opt_result, d, ndds, max_cycle, max_chain)
                print(max_cycle, max_chain, opt_result.total_score, \
                        opt_result.ip_model.obj_val, opt_result_0.total_score, fn)
                assert abs(opt_result.total_score - opt_result_0.total_score) < EPS

def test_weighted_instance_failure_aware():
    """Checks that the CF and PICEF formulations agree on the optimal
    result for an instance with weighted edges, using failure-aware solving.
    """
    EPS = 0.000001
    d, ndds = read_with_ndds("test-fixtures/100-random-weights")
    fns = [k_ip.optimise_picef, k_ip.optimise_ccf]
    for max_cycle in [0, 1, 2, 3, 4]:
        for max_chain in [0, 1, 2, 3]:
            cfg = k_ip.OptConfig(d, ndds, max_cycle, max_chain, edge_success_prob=0.7) 
            opt_result_cf = k_ip.optimise_ccf(cfg)
            opt_result_picef = k_ip.optimise_picef(cfg)
            opt_result_cf_relabelled = k_ip.optimise_relabelled(k_ip.optimise_ccf, cfg)
            opt_result_picef_relabelled = k_ip.optimise_relabelled(k_ip.optimise_picef, cfg)
            assert abs(opt_result_cf.total_score - opt_result_picef.total_score) < EPS
            assert abs(opt_result_cf.total_score - opt_result_cf_relabelled.total_score) < EPS
            assert abs(opt_result_cf.total_score - opt_result_picef_relabelled.total_score) < EPS

def test_relabelled_solve():
    """Checks whether the vertex-ordering heuristic affects the
    result for a weighted instance
    """
    EPS = 0.000001
    d, ndds = read_with_ndds("test-fixtures/100-random-weights")
    for max_cycle in [0, 3]:
        for max_chain in [0, 5]:
            opt_result_0 = k_ip.optimise_picef(k_ip.OptConfig(d, ndds, max_cycle, max_chain))
            print(opt_result_0.total_score)
            opt_result = k_ip.optimise_relabelled(
                    k_ip.optimise_picef, k_ip.OptConfig(d, ndds, max_cycle, max_chain))
            print("   ", opt_result.total_score)
            assert abs(opt_result.total_score - opt_result_0.total_score) < EPS
