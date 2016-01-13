from kidney_digraph import *
from kidney_ndds import *
import kidney_utils

from gurobipy import *

class OptSolution(object):
    def __init__(self, ip_model, cycles, chains, digraph):
        self.ip_model = ip_model
        self.cycles = cycles
        self.chains = chains
        self.digraph = digraph
        self.total_score = (sum(c.score for c in chains) +
                sum(cycle_score(c, digraph) for c in cycles))

    def display(self):
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
        relabelled_cycles = [[old_to_new_vertices[v.id] for v in c] for c in self.cycles]
        relabelled_chains = [Chain(c.ndd_index,
                                   [old_to_new_vertices[i].id for i in c.vtx_indices],
                                   c.score)
                             for c in self.chains]
        return OptSolution(self.ip_model, relabelled_cycles, relabelled_chains, new_digraph)

def create_ip_model(time_limit):
    m = Model("kidney-mip")
    m.params.outputflag = 0
    m.params.mipGap = 0
    if time_limit is not None:
        m.params.timelimit = time_limit
    return m

def add_chain_vars_and_constraints(digraph, ndds, max_chain, m, vtx_to_vars):
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
            for i in range(max_chain-1):
                if dists_from_ndd[e.src.id] <= i+1:
                    edge_var = m.addVar(vtype=GRB.BINARY)
                    e.grb_vars.append(edge_var)
                    vtx_to_vars[e.dest.id].append(edge_var)
                    e.src.grb_vars_out[i].append(edge_var)
                    if i < max_chain-2:
                        e.dest.grb_vars_in[i+1].append(edge_var)

        m.update()

        # At each chain position, sum of edges into a vertex must be >= sum of edges out
        for i in range(max_chain-1):
            for v in digraph.vs:
                m.addConstr(quicksum(v.grb_vars_in[i]) >= quicksum(v.grb_vars_out[i]))

def add_unlimited_vars_and_constraints(digraph, ndds, m):
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
        e.dest.grb_vars_in.append(edge_var)

    m.update()

    for v in digraph.vs:
        if len(v.grb_vars_in) > 1:
            m.addConstr(quicksum(v.grb_vars_in) <= 1)

    # Sum of edges into a vertex must be >= sum of edges out
    for v in digraph.vs:
        m.addConstr(quicksum(v.grb_vars_in) >= quicksum(v.grb_vars_out))

def optimise_uuef(digraph, ndds, max_cycle, max_chain, timelimit):
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
                cycle_next_vv[e.src.id] = e.dest.id

    return OptSolution(ip_model=m,
                       cycles=kidney_utils.selected_edges_to_cycles(digraph, cycle_start_vv, cycle_next_vv),
                       chains=kidney_utils.get_optimal_chains(digraph, ndds),
                       digraph=digraph)
        

def add_hpief_prime_vars_and_constraints(max_cycle, digraph, vtx_to_in_edges, m):
    n = len(digraph.vs)

    E = [[{} for __ in digraph.es] for __ in range(max_cycle)]
    vars_and_edges = []
    
    ses = SortedEdgeStructure(digraph)
    for low_vtx in range(len(digraph.vs)-1):
        # Length of shortest path from low vertex to each vertex with a higher index
        # Default value is 999999999 (which represents infinity)
        shortest_path_from_lv = digraph.get_shortest_path_from_low_vtx(low_vtx, max_cycle-1) 
        shortest_path_to_lv = digraph.get_shortest_path_to_low_vtx(low_vtx, max_cycle-1) 

        for j in xrange(ses.list_starts[low_vtx], len(ses.edges)):
            e = ses.edges[j]
            if e.src.id > low_vtx:
                for pos in xrange(1, max_cycle):
                    if (
                                (e.dest.id==low_vtx or pos < max_cycle-1) and        
                                shortest_path_from_lv[e.src.id] <= pos and
                                shortest_path_to_lv[e.dest.id] < max_cycle - pos):
                        new_var = m.addVar(vtype=GRB.BINARY)
                        E[pos][e.id][low_vtx] = new_var
                        vars_and_edges.append((new_var, pos, e, low_vtx))
    m.update()

    # edge_vars_in[pos][v][low_v] is a list of variables corresponding to edges
    # at position pos that lead to v, in low_v's copy of the graph
    edge_vars_in = [[[[] for __ in range(n)] for __ in range(n)] for __ in range(max_cycle)]

    edge_vars_out = [[[[] for __ in range(n)] for __ in range(n)] for __ in range(max_cycle)]
    
    for pos in range(max_cycle):
        for edge_id in range(len(digraph.es)):
            for low_v_id in E[pos][edge_id]:
                edge_var = E[pos][edge_id][low_v_id]
                edge = digraph.es[edge_id]
                row = edge.src.id
                col = edge.dest.id
                if pos==1:
                    vtx_to_in_edges[row].append(edge_var)
                vtx_to_in_edges[col].append(edge_var)
                if col != low_v_id:
                    edge_vars_in[pos][col][low_v_id].append(edge_var)
                if row != low_v_id:
                    edge_vars_out[pos][row][low_v_id].append(edge_var)
    
    for l in vtx_to_in_edges:
        if len(l) > 0:
            m.addConstr(quicksum(l) <= 1)
    
    for pos in range(1, max_cycle-1):
        for v in range(n):
            for low_v_id in range(v+1):
                in_vars = edge_vars_in[pos][v][low_v_id]
                out_vars = edge_vars_out[pos+1][v][low_v_id]
                if len(in_vars) > 0 or len(out_vars) > 0:
                    m.addConstr(quicksum(in_vars) == quicksum(out_vars))

    return vars_and_edges

def optimise_hpief_prime(digraph, ndds, max_cycle, max_chain, timelimit):
    # This IP model is based on HPIEF, but does not include cycle-edges at position zero.
    
    m = create_ip_model(timelimit)
    m.params.method = 2
    m.params.presolve = 0

    # For each vertex v, a list of variables corresponding to in-edges to v
    vtx_to_in_edges = [[] for __ in digraph.vs]

    add_chain_vars_and_constraints(digraph, ndds, max_chain, m, vtx_to_in_edges)

    vars_and_edges = add_hpief_prime_vars_and_constraints(max_cycle, digraph, vtx_to_in_edges, m)

    obj_terms = []
    for (var, pos, edge, low_v_id) in vars_and_edges:
        if pos==1:
            edge0_score = digraph.adj_mat[low_v_id][edge.src.id].score
            obj_terms.append((edge0_score + edge.score) * var)
        else:
            obj_terms.append(edge.score * var)
    obj_expr = quicksum(obj_terms)
   
    if max_chain > 0:
        obj_expr += quicksum(e.score * e.edge_var for ndd in ndds for e in ndd.edges) 
        obj_expr += quicksum(e.score * var for e in digraph.es for var in e.grb_vars)
    
    m.setObjective(obj_expr, GRB.MAXIMIZE)
    m.optimize()

    cycle_start_vv = []
    cycle_next_vv = {}
    
    for (var, pos, edge, low_v_id) in vars_and_edges:
        if var.x > 0.1:
            cycle_next_vv[edge.src.id] = edge.dest.id
            if pos == 1:
                cycle_start_vv.append(low_v_id)
                cycle_next_vv[low_v_id] = edge.src.id
        
    return OptSolution(ip_model=m,
                       cycles=kidney_utils.selected_edges_to_cycles(digraph, cycle_start_vv, cycle_next_vv),
                       chains=[] if max_chain==0 else kidney_utils.get_optimal_chains(digraph, ndds),
                       digraph=digraph)

def optimise_picef(digraph, ndds, max_cycle, max_chain, timelimit):
    cycles = digraph.find_cycles(max_cycle)

    m = create_ip_model(timelimit)
    m.params.method = 2

    cycle_vars = [m.addVar(vtype=GRB.BINARY) for __ in cycles]
    m.update()
    
    vtx_to_vars = [[] for __ in digraph.vs]
    
    add_chain_vars_and_constraints(digraph, ndds, max_chain, m, vtx_to_vars)

    for i, c in enumerate(cycles):
        for v in c:
            vtx_to_vars[v.id].append(cycle_vars[i])

    for l in vtx_to_vars:
        if len(l) > 0:
            m.addConstr(quicksum(l) <= 1)

    if max_chain==0:
        obj_expr = quicksum(cycle_score(c, digraph) * var for (c, var) in zip(cycles, cycle_vars))
    else:
        obj_expr = ( quicksum(cycle_score(c, digraph) * var for (c, var) in zip(cycles, cycle_vars)) +
                     quicksum(e.score * e.edge_var for ndd in ndds for e in ndd.edges) +
                     quicksum(e.score * var for e in digraph.es for var in e.grb_vars) )

    m.setObjective(obj_expr, GRB.MAXIMIZE)
    m.optimize()

    return OptSolution(ip_model=m,
                       cycles=[c for c, v in zip(cycles, cycle_vars) if v.x > 0.5],
                       chains=[] if max_chain==0 else kidney_utils.get_optimal_chains(digraph, ndds),
                       digraph=digraph)


def optimise_ccf(digraph, ndds, max_cycle, max_chain, timelimit):
    cycles = digraph.find_cycles(max_cycle)
    chains = find_chains(digraph, ndds, max_chain)
        
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

    obj_expr = (quicksum(cycle_score(c, digraph) * var for (c, var) in zip(cycles, cycle_vars)) +
                quicksum(c.score * var for (c, var) in zip(chains, chain_vars)))
        
    m.setObjective(obj_expr, GRB.MAXIMIZE)
    m.optimize()

    return OptSolution(ip_model=m,
                       cycles=[c for c, v in zip(cycles, cycle_vars) if v.x > 0.5],
                       chains=[c for c, v in zip(chains, chain_vars) if v.x > 0.5],
                       digraph=digraph)

def optimise_relabelled(formulation_fun, digraph, ndds, max_cycle, max_chain, timelimit):
    # Sort vertices in descending order of (indegree + outdegree)
    in_degs = [0] * len(digraph.vs)
    for e in digraph.es:
        in_degs[e.dest.id] += 1

    sorted_vertices = sorted(digraph.vs,
                             key=lambda v: len(v.edges) + in_degs[v.id],
                             reverse=True)
    
    relabelled_digraph = digraph.induced_subgraph(sorted_vertices)
    relabelled_digraph.create_adj_mat()

    # old_to_new_vtx[i] is the vertex in the new graph corresponding to vertex
    # i in the original digraph
    old_to_new_vtx = [None] * len(digraph.vs)
    for i, v in enumerate(sorted_vertices):
        old_to_new_vtx[v.id] = relabelled_digraph.vs[i]

    relabelled_ndds = create_relabelled_ndds(ndds, old_to_new_vtx)

    opt_result = formulation_fun(relabelled_digraph, relabelled_ndds,
                                 max_cycle, max_chain, timelimit)

    return opt_result.relabelled_copy(sorted_vertices, digraph)
