"""
IP formulations for kidney exchange, including PICEF
"""

import argparse
import time
import sys

from . import kidney_digraph
from . import kidney_ip
from . import kidney_utils
from . import kidney_ndds

def solve_kep(cfg, formulation, use_relabelled=True):

    formulations = {
        "uef":  ("Uncapped edge formulation", kidney_ip.optimise_uuef),
        "eef": ("EEF", kidney_ip.optimise_eef),
        "eef_full_red": ("EEF with full reduction by cycle generation", kidney_ip.optimise_eef_full_red),
        "hpief_prime": ("HPIEF'", kidney_ip.optimise_hpief_prime),
        "hpief_prime_full_red": ("HPIEF' with full reduction by cycle generation", kidney_ip.optimise_hpief_prime_full_red),
        "hpief_2prime": ("HPIEF''", kidney_ip.optimise_hpief_2prime),
        "hpief_2prime_full_red": ("HPIEF'' with full reduction by cycle generation", kidney_ip.optimise_hpief_2prime_full_red),
        "picef": ("PICEF", kidney_ip.optimise_picef),
        "cf":   ("Cycle formulation",
                  kidney_ip.optimise_ccf)
    }
    
    if formulation in formulations:
        formulation_name, formulation_fun = formulations[formulation]
        if use_relabelled:
            opt_result = kidney_ip.optimise_relabelled(formulation_fun, cfg)
        else:
            opt_result = formulation_fun(cfg)
        kidney_utils.check_validity(opt_result, cfg.digraph, cfg.ndds, cfg.max_cycle, cfg.max_chain)
        opt_result.formulation_name = formulation_name
        return opt_result
    else:
        raise ValueError("Unrecognised IP formulation name")

def start():
    parser = argparse.ArgumentParser("Solve a kidney-exchange instance")
    parser.add_argument("cycle_cap", type=int,
            help="The maximum permitted cycle length")
    parser.add_argument("chain_cap", type=int,
            help="The maximum permitted number of edges in a chain")
    parser.add_argument("formulation",
            help="The IP formulation (uef, eef, eef_full_red, hpief_prime, hpief_2prime, hpief_prime_full_red, hpief_2prime_full_red, picef, cf)")
    parser.add_argument("--use-relabelled", "-r", required=False,
            action="store_true",
            help="Relabel vertices in descending order of in-deg + out-deg")
    parser.add_argument("--eef-alt-constraints", "-e", required=False,
            action="store_true",
            help="Use slightly-modified EEF constraints (ignored for other formulations)")
    parser.add_argument("--timelimit", "-t", required=False, default=None,
            type=float,
            help="IP solver time limit in seconds (default: no time limit)")
    parser.add_argument("--verbose", "-v", required=False,
            action="store_true",
            help="Log Gurobi output to screen and log file")
    parser.add_argument("--edge-success-prob", "-p", required=False,
            type=float, default=1.0,
            help="Edge success probability, for failure-aware matching. " +
                 "This can only be used with PICEF and cycle formulation. (default: 1)")
    parser.add_argument("--lp-file", "-l", required=False, default=None,
            metavar='FILE',
            help="Write the IP model to FILE, then exit.")
    parser.add_argument("--relax", "-x", required=False,
            action='store_true',
            help="Solve the LP relaxation.")
            
    args = parser.parse_args()
    args.formulation = args.formulation.lower()

    input_lines = [line for line in sys.stdin if len(line.strip()) > 0]
    n_digraph_edges = int(input_lines[0].split()[1])
    digraph_lines = input_lines[:n_digraph_edges + 2]

    d = kidney_digraph.read_digraph(digraph_lines)

    if len(input_lines) > len(digraph_lines):
        ndd_lines = input_lines[n_digraph_edges + 2:]
        altruists = kidney_ndds.read_ndds(ndd_lines, d)
    else:
        altruists = []
        
    start_time = time.time()
    cfg = kidney_ip.OptConfig(d, altruists, args.cycle_cap, args.chain_cap, args.verbose,
                              args.timelimit, args.edge_success_prob, args.eef_alt_constraints,
                              args.lp_file, args.relax)
    opt_solution = solve_kep(cfg, args.formulation, args.use_relabelled)
    time_taken = time.time() - start_time
    print(("formulation: {}".format(args.formulation)))
    print(("formulation_name: {}".format(opt_solution.formulation_name)))
    print(("using_relabelled: {}".format(args.use_relabelled)))
    print(("cycle_cap: {}".format(args.cycle_cap)))
    print(("chain_cap: {}".format(args.chain_cap)))
    print(("edge_success_prob: {}".format(args.edge_success_prob)))
    print(("ip_time_limit: {}".format(args.timelimit)))
    print(("ip_vars: {}".format(opt_solution.ip_model.numVars)))
    print(("ip_constrs: {}".format(opt_solution.ip_model.numConstrs)))
    print(("total_time: {}".format(time_taken)))
    print(("ip_solve_time: {}".format(opt_solution.ip_model.runtime)))
    print(("solver_status: {}".format(opt_solution.ip_model.status)))
    print(("total_score: {}".format(opt_solution.total_score)))
    opt_solution.display()

if __name__=="__main__":
    start()
