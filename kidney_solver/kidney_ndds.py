class Ndd:
    def __init__(self):
        self.edges = []
    def add_edge(self, ndd_edge):
        self.edges.append(ndd_edge)

class NddEdge:
    def __init__(self, target_v, score):
        self.target_v = target_v
        self.score = score

def create_relabelled_ndds(ndds, old_to_new_vtx):
    """Creates a copy of a n array of NDDs, with target vertices changed.

    If a target vertex in an original NDD had ID i, then the new target vertex
    will be old_to_new_vtx[i].
    """
    new_ndds = [Ndd() for ndd in ndds]
    for i, ndd in enumerate(ndds):
        for edge in ndd.edges:
            new_ndds[i].add_edge(NddEdge(old_to_new_vtx[edge.target_v.id], edge.score))

    return new_ndds

def read_ndds(lines, digraph):
    ndds = []
    ndd_count, edge_count = [int(x) for x in lines[0].split()]
    ndds = [Ndd() for _ in range(ndd_count)]
    for line in lines[1:edge_count+1]:
        tokens = [t for t in line.split()]
        src_id = int(tokens[0])
        dest_id = int(tokens[1])
        score = float(tokens[2])
        ndds[src_id].add_edge(NddEdge(digraph.vs[dest_id], score))
    return ndds

class Chain(object):
    def __init__(self, ndd_index, vtx_indices, score):
        self.ndd_index = ndd_index
        self.vtx_indices = vtx_indices
        self.score = score

    def __repr__(self):
        return ("Chain NDD{} ".format(self.ndd_index) +
                        " ".join(str(v) for v in self.vtx_indices) +
                        " with score " + str(self.score))

    def __cmp__(self, other):
        # Compare on NDD ID, then chain length, then on score, then
        # lexicographically on vtx indices
        if self.ndd_index < other.ndd_index:
            return -1
        elif self.ndd_index > other.ndd_index:
            return 1
        elif len(self.vtx_indices) < len(other.vtx_indices):
            return -1
        elif len(self.vtx_indices) > len(other.vtx_indices):
            return 1
        elif self.score < other.score:
            return -1
        elif self.score > other.score:
            return 1
        else:
            for i, j in zip(self.vtx_indices, other.vtx_indices):
                if i < j:
                    return -1
                elif i > j:
                    return 1
        return 0
            
def find_chains(digraph, ndds, max_chain):
    def find_chains_recurse(vertices, score):
        chains.append(Chain(ndd_idx, vertices[:], score))
        if len(vertices) < max_chain:
            for e in digraph.vs[vertices[-1]].edges:
                if e.dest.id not in vertices:
                    vertices.append(e.dest.id)
                    find_chains_recurse(vertices, score+e.score)
                    del vertices[-1]
    chains = []
    if max_chain == 0:
        return chains
    for ndd_idx, ndd in enumerate(ndds):
        for e in ndd.edges:
            vertices = [e.target_v.id]
            find_chains_recurse(vertices, e.score)
    return chains

