"""
IP formulations for kidney exchange, including PICEF
"""

import argparse
import time
import sys

import kidney_digraph
import kidney_ip
import kidney_utils
import kidney_ndds

def solve_kep(digraph, ndds, max_cycle, max_chain, formulation, timelimit, use_relabelled=True):

    formulations = {
        "uef":  ("Uncapped edge formulation", kidney_ip.optimise_uuef),
        "hpief_prime": ("HPIEF'", kidney_ip.optimise_hpief_prime),
        "picef": ("PICEF", kidney_ip.optimise_picef),
        "cf":   ("Cycle formulation",
                  kidney_ip.optimise_ccf)
    }
    
    if formulation in formulations:
        formulation_name, formulation_fun = formulations[formulation]
        if use_relabelled:
            opt_result = kidney_ip.optimise_relabelled(formulation_fun,
                    digraph, ndds, max_cycle, max_chain, timelimit)
        else:
            opt_result = formulation_fun(digraph, ndds, max_cycle, max_chain, timelimit)
        kidney_utils.check_validity(opt_result, digraph, ndds, max_cycle, max_chain)
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
            help="The IP formulation (uef, hpief_prime, picef, cf)")
    parser.add_argument("--use-relabelled", "-r", required=False,
            action="store_true",
            help="Relabel vertices in descending order of in-deg + out-deg")
    parser.add_argument("--timelimit", "-t", required=False, default=None,
            type=float,
            help="IP solver time limit in seconds (default: no time limit)")
            
    args = parser.parse_args()
    args.formulation = args.formulation.lower()

    input_lines = [line for line in sys.stdin]
    n_digraph_edges = int(input_lines[0].split()[1])
    digraph_lines = input_lines[:n_digraph_edges + 2]

    d = kidney_digraph.read_digraph_without_prob(digraph_lines)

    if len(input_lines) > len(digraph_lines):
        ndd_lines = input_lines[n_digraph_edges + 2:]
        altruists = kidney_ndds.read_ndds(ndd_lines, d)
    else:
        altruists = []
        
    start_time = time.time()
    opt_solution = solve_kep(d, altruists, args.cycle_cap, args.chain_cap,
                             args.formulation, args.timelimit, args.use_relabelled)
    time_taken = time.time() - start_time
    print "formulation: {}".format(args.formulation)
    print "formulation_name: {}".format(opt_solution.formulation_name)
    print "using_relabelled: {}".format(args.use_relabelled)
    print "cycle_cap: {}".format(args.cycle_cap)
    print "chain_cap: {}".format(args.chain_cap)
    print "ip_time_limit: {}".format(args.timelimit)
    print "ip_vars: {}".format(opt_solution.ip_model.numVars)
    print "ip_constrs: {}".format(opt_solution.ip_model.numConstrs)
    print "total_time: {}".format(time_taken)
    print "ip_solve_time: {}".format(opt_solution.ip_model.runtime)
    print "solver_status: {}".format(opt_solution.ip_model.status)
    print "total_score: {}".format(opt_solution.total_score)
    opt_solution.display()

if __name__=="__main__":
    start()
