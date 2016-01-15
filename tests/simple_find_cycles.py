def simple_find_cycles(digraph, max_length):
    # A simple, hopefully-bug-free cycle finding algorithm, to check
    # the results of faster algorithms against.

    # Returns a list of lists of vertices
    
    cycles = []
    vtx_used = [False] * len(digraph.vs)

    def cycle(digraph, low_vtx_id, current_list):
        last_vtx = current_list[-1]
        if digraph.edge_exists(last_vtx, current_list[0]):
            cycles.append(current_list[:])
        if len(current_list) < max_length:
            for edge in last_vtx.edges:
                v = edge.tgt
                if v.id > low_vtx_id and not vtx_used[v.id]:
                    vtx_used[v.id] = True
                    current_list.append(v)
                    cycle(digraph, low_vtx_id, current_list)
                    vtx_used[v.id] = False
                    del current_list[-1]

    for i in range(0, len(digraph.vs)):
        vtx_used[i] = True
        cycle(digraph, i, [digraph.vs[i]])
        vtx_used[i] = False
    return cycles


