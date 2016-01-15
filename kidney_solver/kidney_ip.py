"""Solving the kidney-exchange problem using the Gurobi IP solver."""

from kidney_digraph import *
from kidney_ndds import *
import kidney_utils

from gurobipy import *

###################################################################################################
#                                                                                                 #
#                                  Code used by all formulations                                  #
#                                                                                                 #
###################################################################################################

class OptSolution(object):
    """An optimal solution for a kidney-exchange problem instance.
    
    Data members:
        ip_model: The Gurobi Model object
        cycles: A list of cycles in the optimal solution, each represented
            as a list of vertices
        chains: A list of chains in the optimal solution, each represented
            as a Chain object
        total_score: The total score of the solution
    """

    def __init__(self, ip_model, cycles, chains, digraph, edge_success_prob=1):
        self.ip_model = ip_model
        self.cycles = cycles
        self.chains = chains
        self.digraph = digraph
        self.total_score = (sum(c.score for c in chains) +
                sum(failure_aware_cycle_score(c, digraph, edge_success_prob) for c in cycles))

    def display(self):
        """Print the optimal cycles and chains to standard output."""

        print "cycle_count: {}".format(len(self.cycles))
        print "chain_count: {}".format(len(self.chains))
        print "cycles:"
        # cs is a list of cycles, with each cycle represented as a list of vertex IDs
        cs = [[v.id for v in c] for c in self.cycles]
        # Put the lowest-indexed vertex at the start of each cycle
        for i in range(len(cs)):
            min_index_pos = cs[i].index(min(cs[i]))
            cs[i] = cs[i][min_index_pos:] + cs[i][:min_index_pos]
        # Sort the cycles
        cs.sort()
        for c in cs:
            print "\t".join(str(v_id) for v_id in c)
        print "chains:"
        for c in self.chains:
            print str(c.ndd_index) + "\t" + "\t".join(str(v) for v in c.vtx_indices)

    def relabelled_copy(self, old_to_new_vertices, new_digraph):
        """Create a copy of the solution with vertices relabelled.

        If the solution was found on a relabelled copy of the instance digraph, this
        method can be used to transform the solution back to the original digraph. Each
        Vertex v in the OptSolution on which this method is called is replaced in the
        returned copy by old_to_new_vertices[v.id].
        """

        relabelled_cycles = [[old_to_new_vertices[v.id] for v in c] for c in self.cycles]
        relabelled_chains = [Chain(c.ndd_index,
                                   [old_to_new_vertices[i].id for i in c.vtx_indices],
                                   c.score)
                             for c in self.chains]
        return OptSolution(self.ip_model, relabelled_cycles, relabelled_chains, new_digraph)

def optimise_relabelled(formulation_fun, digraph, ndds, max_cycle, max_chain, timelimit, edge_success_prob=1):
    """Optimise on a relabelled graph such that vertices are sorted in descending
        order of (indegree + outdegree)"""

    in_degs = [0] * len(digraph.vs)
    for e in digraph.es:
        in_degs[e.tgt.id] += 1

    sorted_vertices = sorted(digraph.vs,
                             key=lambda v: len(v.edges) + in_degs[v.id],
                             reverse=True)
    
    relabelled_digraph = digraph.induced_subgraph(sorted_vertices)

    # old_to_new_vtx[i] is the vertex in the new graph corresponding to vertex
    # i in the original digraph
    old_to_new_vtx = [None] * len(digraph.vs)
    for i, v in enumerate(sorted_vertices):
        old_to_new_vtx[v.id] = relabelled_digraph.vs[i]

    relabelled_ndds = create_relabelled_ndds(ndds, old_to_new_vtx)

    opt_result = formulation_fun(relabelled_digraph, relabelled_ndds,
                                 max_cycle, max_chain, timelimit, edge_success_prob)

    return opt_result.relabelled_copy(sorted_vertices, digraph)

def create_ip_model(time_limit):
    """Create a Gurobi Model."""

    m = Model("kidney-mip")
    m.params.outputflag = 0
    m.params.mipGap = 0
    if time_limit is not None:
        m.params.timelimit = time_limit
    return m

###################################################################################################
#                                                                                                 #
#                                       Uncapped formulation                                      #
#                                                                                                 #
###################################################################################################

def add_unlimited_vars_and_constraints(digraph, ndds, m):
    """Add the IP variables and constraints for chains in the uncapped edge formulation. 

    Args:
        digraph: the instance digraph
        ndds: a list of NDDs in the instance
        m: The Gurobi model
    """

    for v in digraph.vs:
        v.grb_vars_in  = []
        v.grb_vars_out = []

    for ndd in ndds:
        ndd_edge_vars = []
        for e in ndd.edges:
            edge_var = m.addVar(vtype=GRB.BINARY)
            e.edge_var = edge_var
            ndd_edge_vars.append(edge_var)
            e.target_v.grb_vars_in.append(edge_var)
        m.update()
        m.addConstr(quicksum(ndd_edge_vars) <= 1)

    # Add pair->pair edge variables
    for e in digraph.es:
        e.grb_vars = []
        edge_var = m.addVar(vtype=GRB.BINARY)
        e.grb_vars.append(edge_var)
        e.src.grb_vars_out.append(edge_var)
        e.tgt.grb_vars_in.append(edge_var)

    m.update()

    for v in digraph.vs:
        if len(v.grb_vars_in) > 1:
            m.addConstr(quicksum(v.grb_vars_in) <= 1)

    # Sum of edges into a vertex must be >= sum of edges out
    for v in digraph.vs:
        m.addConstr(quicksum(v.grb_vars_in) >= quicksum(v.grb_vars_out))

def optimise_uuef(digraph, ndds, max_cycle, max_chain, timelimit, edge_success_prob=1):
    """Optimise using the uncapped edge formulation.

    Args:
        ndds: NDDs in the instance
        max_cycle: the cycle cap
        max_chain: the chain cap
        timelimit: the Gurobi timeout in seconds, or None for no timeout

    Returns:
        an OptSolution object
    """

    if edge_success_prob != 1:
        raise ValueError("This formulation does not support failure-aware matching.")

    m = create_ip_model(timelimit)

    add_unlimited_vars_and_constraints(digraph, ndds, m)

    obj_expr = ( quicksum(e.score * e.edge_var for ndd in ndds for e in ndd.edges) +
                 quicksum(e.score * var for e in digraph.es for var in e.grb_vars) )
   
    m.setObjective(obj_expr, GRB.MAXIMIZE)
    m.optimize()

    # Try all possible cycle start positions
    cycle_start_vv = range(len(digraph.vs))

    cycle_next_vv = {}
    for e in digraph.es:
        for var in e.grb_vars:
            if var.x > 0.1:
                cycle_next_vv[e.src.id] = e.tgt.id

    return OptSolution(ip_model=m,
                       cycles=kidney_utils.selected_edges_to_cycles(digraph, cycle_start_vv, cycle_next_vv),
                       chains=kidney_utils.get_optimal_chains(digraph, ndds),
                       digraph=digraph)
        
###################################################################################################
#                                                                                                 #
#                  Chain vars and constraints (used by HPIEF', HPIEF'' and PICEF)                 #
#                                                                                                 #
###################################################################################################

def add_chain_vars_and_constraints(digraph, ndds, max_chain, m, vtx_to_vars, store_edge_positions=False):
    """Add the IP variables and constraints for chains in PICEF and HPIEF'.

    Args:
        ndds: a list of NDDs in the instance
        max_chain: the chain cap
        m: The Gurobi model
        vtx_to_vars: A list such that for each Vertex v in the Digraph,
            vtx_to_vars[v.id] will contain the Gurobi variables representing
            edges pointing to v.
    """

    if max_chain > 0:
        for v in digraph.vs:
            v.grb_vars_in  = [[] for i in range(max_chain-1)]
            v.grb_vars_out = [[] for i in range(max_chain-1)]

        for ndd in ndds:
            ndd_edge_vars = []
            for e in ndd.edges:
                edge_var = m.addVar(vtype=GRB.BINARY)
                e.edge_var = edge_var
                ndd_edge_vars.append(edge_var)
                vtx_to_vars[e.target_v.id].append(edge_var)
                if max_chain>1: e.target_v.grb_vars_in[0].append(edge_var)
            m.update()
            m.addConstr(quicksum(ndd_edge_vars) <= 1)

        dists_from_ndd = kidney_utils.get_dist_from_nearest_ndd(digraph, ndds)

        # Add pair->pair edge variables, indexed by position in chain
        for e in digraph.es:
            e.grb_vars = []
            if store_edge_positions:
                e.grb_var_positions = []
            for i in range(max_chain-1):
                if dists_from_ndd[e.src.id] <= i+1:
                    edge_var = m.addVar(vtype=GRB.BINARY)
                    e.grb_vars.append(edge_var)
                    if store_edge_positions:
                        e.grb_var_positions.append(i+1)
                    vtx_to_vars[e.tgt.id].append(edge_var)
                    e.src.grb_vars_out[i].append(edge_var)
                    if i < max_chain-2:
                        e.tgt.grb_vars_in[i+1].append(edge_var)

        m.update()

        # At each chain position, sum of edges into a vertex must be >= sum of edges out
        for i in range(max_chain-1):
            for v in digraph.vs:
                m.addConstr(quicksum(v.grb_vars_in[i]) >= quicksum(v.grb_vars_out[i]))

###################################################################################################
#                                                                                                 #
#                                Code shared by HPIEF' and HPIEF''                                #
#                                                                                                 #
###################################################################################################

def add_hpief_prime_vars_partial_red(max_cycle, digraph, m, hpief_2_prime=False):
    vars_and_edges = [] # A list of (gurobi_var, position, edge, low_vertex) tuples

    # max_pos is the maximum edge position for which variables may be created
    max_pos = max_cycle-2 if hpief_2_prime else max_cycle-1
    
    # Index i is in the list edge_vars_in[pos][v][low_v] if and only if
    # vars_and_edges[i] corresponds to an edge at position pos, pointing to vertex
    # v, in low_v's graph copy 
    edge_vars_in = [[[[] for __ in range(digraph.n)] for __ in range(digraph.n)] for __ in range(max_pos + 1)]

    # Index i is in the list edge_vars_out[pos][v][low_v] if and only if
    # vars_and_edges[i] corresponds to an edge at position pos, leaving vertex
    # v, in low_v's graph copy 
    edge_vars_out = [[[[] for __ in range(digraph.n)] for __ in range(digraph.n)] for __ in range(max_pos + 1)]

    for low_vtx in range(digraph.n-2 if hpief_2_prime else digraph.n-1):
        # Length of shortest path from low vertex to each vertex with a higher index
        # Default value is 999999999 (which represents infinity)
        shortest_path_from_lv = digraph.get_shortest_path_from_low_vtx(low_vtx, max_cycle-1) 
        shortest_path_to_lv = digraph.get_shortest_path_to_low_vtx(low_vtx, max_cycle-1) 

        for v1 in digraph.vs[low_vtx+1:]:
            for e in v1.edges:
                if e.tgt.id >=low_vtx:
                    for pos in xrange(1, max_pos + 1):
                        if (shortest_path_from_lv[e.src.id] <= pos and
                                    shortest_path_to_lv[e.tgt.id] < max_cycle - pos):
                            new_var = m.addVar(vtype=GRB.BINARY)
                            vars_and_edges.append((new_var, pos, e, low_vtx))
                            idx = len(vars_and_edges) - 1 # Index of tuple just added
                            edge_vars_in[pos][e.tgt.id][low_vtx].append(idx)
                            edge_vars_out[pos][e.src.id][low_vtx].append(idx)
    m.update()
    return vars_and_edges, edge_vars_in, edge_vars_out

def add_hpief_prime_vars_full_red(max_cycle, digraph, m, hpief_2_prime=False):
    vars_and_edges = [] # A list of (gurobi_var, position, edge, low_vertex) tuples

    # max_pos is the maximum edge position for which variables may be created
    max_pos = max_cycle-2 if hpief_2_prime else max_cycle-1
    
    edge_vars_in = [[[[] for __ in range(digraph.n)] for __ in range(digraph.n)] for __ in range(max_pos + 1)]
    edge_vars_out = [[[[] for __ in range(digraph.n)] for __ in range(digraph.n)] for __ in range(max_pos + 1)]

    edges_seen = set()  # (low_v_id, src_v_id, tgt_v_id, pos) tuples
    for cycle in digraph.generate_cycles(max_cycle):
        for i in range(1, len(cycle)-1):
            edges_seen.add((cycle[0].id, cycle[i].id, cycle[i+1].id, i))
        if not hpief_2_prime or len(cycle) < max_cycle:
            edges_seen.add((cycle[0].id, cycle[-1].id, cycle[0].id, len(cycle)-1))
            
    for low_v, src_v, tgt_v, pos in edges_seen:
        new_var = m.addVar(vtype=GRB.BINARY)
        e = digraph.adj_mat[src_v][tgt_v]
        vars_and_edges.append((new_var, pos, e, low_v))
        idx = len(vars_and_edges) - 1 # Index of tuple just added
        edge_vars_in[pos][tgt_v][low_v].append(idx)
        edge_vars_out[pos][src_v][low_v].append(idx)
    m.update()
    return vars_and_edges, edge_vars_in, edge_vars_out

def add_hpief_prime_vars_and_constraints(max_cycle, digraph, vtx_to_in_edges, m, full_red, hpief_2_prime=False):
    max_pos = max_cycle-2 if hpief_2_prime else max_cycle-1

    if full_red:
        vars_and_edges, edge_vars_in, edge_vars_out = add_hpief_prime_vars_full_red(max_cycle, digraph, m, hpief_2_prime)
    else:
        vars_and_edges, edge_vars_in, edge_vars_out = add_hpief_prime_vars_partial_red(max_cycle, digraph, m, hpief_2_prime)
    
    for grb_var, pos, edge, low_vtx in vars_and_edges:
        vtx_to_in_edges[edge.tgt.id].append(grb_var)
        if pos==1:
            vtx_to_in_edges[edge.src.id].append(grb_var)
        if hpief_2_prime and pos == max_cycle - 2 and edge.tgt.id != low_vtx:
            vtx_to_in_edges[low_vtx].append(grb_var)
        
    # Capacity constraint for vertices
    for l in vtx_to_in_edges:
        if len(l) > 0:
            m.addConstr(quicksum(l) <= 1)
    
    # Cycle flow-conservation constraint for vertices
    for pos in range(1, max_pos):
        for v in range(digraph.n):
            for low_v_id in range(v):
                in_vars  = [vars_and_edges[i][0] for i in edge_vars_in[pos][v][low_v_id]]
                out_vars = [vars_and_edges[i][0] for i in edge_vars_out[pos+1][v][low_v_id]]
                if len(in_vars) > 0 or len(out_vars) > 0:
                    m.addConstr(quicksum(in_vars) == quicksum(out_vars))

    return vars_and_edges

def optimise_hpief_prime(digraph, ndds, max_cycle, max_chain, timelimit, edge_success_prob=1, full_red=False, hpief_2_prime=False):
    """Optimise using the HPIEF' or HPIEF'' formulation.

    The HPIEF' model is based on HPIEF, but does not include cycle-edge variables at position zero.
    HPIEF'' also removes variables corresponding to edges at the last possible position of a cycle. 

    Args:
        ndds: NDDs in the instance
        max_cycle: the cycle cap
        max_chain: the chain cap
        timelimit: the Gurobi timeout in seconds, or None for no timeout
        full_red: True if cycles should be generated in order to reduce number of variables further
        hpief_2_prime: Use HPIEF''? Default: HPIEF'

    Returns:
        an OptSolution object
    """

    if edge_success_prob != 1:
        raise ValueError("This formulation does not support failure-aware matching.")

    if max_cycle < 3:
        hpief_2_prime = False

    m = create_ip_model(timelimit)
    m.params.method = 2
    m.params.presolve = 0

    # For each vertex v, a list of variables corresponding to in-edges to v
    vtx_to_in_edges = [[] for __ in digraph.vs]

    add_chain_vars_and_constraints(digraph, ndds, max_chain, m, vtx_to_in_edges)

    vars_and_edges = add_hpief_prime_vars_and_constraints(max_cycle, digraph, vtx_to_in_edges, m, full_red, hpief_2_prime)

    obj_terms = []
    for var, pos, edge, low_v_id in vars_and_edges:
        score = edge.score
        if pos==1:
            score += digraph.adj_mat[low_v_id][edge.src.id].score
        if hpief_2_prime and pos==max_cycle - 2 and edge.tgt.id != low_v_id:
            score += digraph.adj_mat[edge.tgt.id][low_v_id].score
        obj_terms.append(score * var)

    obj_expr = quicksum(obj_terms)
   
    if max_chain > 0:
        obj_expr += quicksum(e.score * e.edge_var for ndd in ndds for e in ndd.edges) 
        obj_expr += quicksum(e.score * var for e in digraph.es for var in e.grb_vars)
    
    m.setObjective(obj_expr, GRB.MAXIMIZE)
    m.optimize()

    cycle_start_vv = []
    cycle_next_vv = {}
    
    for var, pos, edge, low_v_id in vars_and_edges:
        if var.x > 0.1:
            cycle_next_vv[edge.src.id] = edge.tgt.id
            if pos == 1:
                cycle_start_vv.append(low_v_id)
                cycle_next_vv[low_v_id] = edge.src.id
            if hpief_2_prime and pos == max_cycle - 2 and edge.tgt.id != low_v_id:
                cycle_next_vv[edge.tgt.id] = low_v_id
        
    return OptSolution(ip_model=m,
                       cycles=kidney_utils.selected_edges_to_cycles(digraph, cycle_start_vv, cycle_next_vv),
                       chains=[] if max_chain==0 else kidney_utils.get_optimal_chains(digraph, ndds),
                       digraph=digraph)

###################################################################################################
#                                                                                                 #
#                                               HPIEF'                                            #
#                                                                                                 #
###################################################################################################

def optimise_hpief_prime_full_red(digraph, ndds, max_cycle, max_chain, timelimit, edge_success_prob=1):
    return optimise_hpief_prime(digraph, ndds, max_cycle, max_chain, timelimit, edge_success_prob, True)

###################################################################################################
#                                                                                                 #
#                                             HPIEF''                                             #
#                                                                                                 #
###################################################################################################

def optimise_hpief_2prime(digraph, ndds, max_cycle, max_chain, timelimit, edge_success_prob=1, full_red=False):
    return optimise_hpief_prime(digraph, ndds, max_cycle, max_chain, timelimit,
            edge_success_prob, full_red, hpief_2_prime=True)

def optimise_hpief_2prime_full_red(digraph, ndds, max_cycle, max_chain, timelimit, edge_success_prob=1):
    return optimise_hpief_2prime(digraph, ndds, max_cycle, max_chain, timelimit, edge_success_prob, full_red=True)

###################################################################################################
#                                                                                                 #
#                                              PICEF                                              #
#                                                                                                 #
###################################################################################################

def optimise_picef(digraph, ndds, max_cycle, max_chain, timelimit, edge_success_prob=1):
    """Optimise using the PICEF formulation.

    Args:
        ndds: NDDs in the instance
        max_cycle: the cycle cap
        max_chain: the chain cap
        timelimit: the Gurobi timeout in seconds, or None for no timeout

    Returns:
        an OptSolution object
    """

    cycles = digraph.find_cycles(max_cycle)

    m = create_ip_model(timelimit)
    m.params.method = 2

    cycle_vars = [m.addVar(vtype=GRB.BINARY) for __ in cycles]
    m.update()
    
    vtx_to_vars = [[] for __ in digraph.vs]
    
    add_chain_vars_and_constraints(digraph, ndds, max_chain, m, vtx_to_vars, store_edge_positions=edge_success_prob!=1)

    for i, c in enumerate(cycles):
        for v in c:
            vtx_to_vars[v.id].append(cycle_vars[i])

    for l in vtx_to_vars:
        if len(l) > 0:
            m.addConstr(quicksum(l) <= 1)

    if max_chain==0:
        obj_expr = quicksum(failure_aware_cycle_score(c, digraph, edge_success_prob) * var for c, var in zip(cycles, cycle_vars))
    else:
        if edge_success_prob == 1:
            obj_expr = ( quicksum(cycle_score(c, digraph) * var for c, var in zip(cycles, cycle_vars)) +
                         quicksum(e.score * e.edge_var for ndd in ndds for e in ndd.edges) +
                         quicksum(e.score * var for e in digraph.es for var in e.grb_vars) )
        else:
            obj_expr = ( quicksum(failure_aware_cycle_score(c, digraph, edge_success_prob) * var for c, var in zip(cycles, cycle_vars)) +
                         quicksum(e.score*edge_success_prob * e.edge_var for ndd in ndds for e in ndd.edges) +
                         quicksum(e.score*edge_success_prob**(pos+1) * var
                                for e in digraph.es for var, pos in zip(e.grb_vars, e.grb_var_positions)) )

    m.setObjective(obj_expr, GRB.MAXIMIZE)
    m.optimize()

    return OptSolution(ip_model=m,
                       cycles=[c for c, v in zip(cycles, cycle_vars) if v.x > 0.5],
                       chains=[] if max_chain==0 else kidney_utils.get_optimal_chains(digraph, ndds, edge_success_prob),
                       digraph=digraph,
                       edge_success_prob=edge_success_prob)

###################################################################################################
#                                                                                                 #
#                                        Cycle formulation                                        #
#                                                                                                 #
###################################################################################################

def optimise_ccf(digraph, ndds, max_cycle, max_chain, timelimit, edge_success_prob=1):
    """Optimise using the cycle formulation (with one var per cycle and one var per chain).

    Args:
        ndds: NDDs in the instance
        max_cycle: the cycle cap
        max_chain: the chain cap
        timelimit: the Gurobi timeout in seconds, or None for no timeout

    Returns:
        an OptSolution object
    """

    cycles = digraph.find_cycles(max_cycle)
    chains = find_chains(digraph, ndds, max_chain, edge_success_prob)
        
    m = create_ip_model(timelimit)
    m.params.method = 2

    cycle_vars = [m.addVar(vtype=GRB.BINARY) for __ in cycles]
    chain_vars = [m.addVar(vtype=GRB.BINARY) for __ in chains]
    m.update()
    
    ndd_to_vars = [[] for __ in ndds]
    vtx_to_vars = [[] for __ in digraph.vs]
    
    for i, c in enumerate(cycles):
        for v in c:
            vtx_to_vars[v.id].append(cycle_vars[i])

    for i, c in enumerate(chains):
        ndd_to_vars[c.ndd_index].append(chain_vars[i])
        for v in c.vtx_indices:
            vtx_to_vars[v].append(chain_vars[i])

    # Each donor-patient pair and each each NDD is in at most one chosen cycle or chain
    for l in vtx_to_vars + ndd_to_vars:
        if len(l) > 0:
            m.addConstr(quicksum(l) <= 1)

    obj_expr = (quicksum(failure_aware_cycle_score(c, digraph, edge_success_prob) * var for (c, var) in zip(cycles, cycle_vars)) +
                quicksum(c.score * var for (c, var) in zip(chains, chain_vars)))
        
    m.setObjective(obj_expr, GRB.MAXIMIZE)
    m.optimize()

    return OptSolution(ip_model=m,
                       cycles=[c for c, v in zip(cycles, cycle_vars) if v.x > 0.5],
                       chains=[c for c, v in zip(chains, chain_vars) if v.x > 0.5],
                       digraph=digraph,
                       edge_success_prob=edge_success_prob)
