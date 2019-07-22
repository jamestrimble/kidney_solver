"""This module has objects for non-directed donors (NDDs) and chains initiated by NDDs.

In contrast to many kidney-exchange formulations, we do not include in the main directed
graph a vertex for each NDD. Rather, the program maintains a separate array
of Ndd objects. Each of these Ndd objects maintains a list of outgoing edges, and
each of these edges points towards a vertex (which represents a donor-patient pair)
in the directed graph.
"""

from .kidney_digraph import KidneyReadException

class Ndd:
    """A non-directed donor"""
    def __init__(self):
        self.edges = []
    def add_edge(self, ndd_edge):
        """Add an edge representing compatibility with a patient who appears as a
        vertex in the directed graph."""
        self.edges.append(ndd_edge)

class NddEdge:
    """An edge pointing from an NDD to a vertex in the directed graph"""
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
    """Reads NDDs from an array of strings in the .ndd format."""

    ndds = []
    ndd_count, edge_count = [int(x) for x in lines[0].split()]
    ndds = [Ndd() for _ in range(ndd_count)]

    # Keep track of which edges have been created already so that we can
    # detect duplicates
    edge_exists = [[False for v in digraph.vs] for ndd in ndds]

    for line in lines[1:edge_count+1]:
        tokens = [t for t in line.split()]
        src_id = int(tokens[0])
        tgt_id = int(tokens[1])
        score = float(tokens[2])
        if src_id < 0 or src_id >= ndd_count:
            raise KidneyReadException("NDD index {} out of range.".format(src_id))
        if tgt_id < 0 or tgt_id >= digraph.n:
            raise KidneyReadException("Vertex index {} out of range.".format(tgt_id))
        if edge_exists[src_id][tgt_id]:
            raise KidneyReadException(
                    "Duplicate edge from NDD {0} to vertex {1}.".format(src_id, tgt_id))
        ndds[src_id].add_edge(NddEdge(digraph.vs[tgt_id], score))
        edge_exists[src_id][tgt_id] = True

    if lines[edge_count+1].split()[0] != "-1" or len(lines) < edge_count+2:
        raise KidneyReadException("Incorrect edge count")

    return ndds

class Chain(object):
    """A chain initiated by an NDD.
    
    Data members:
        ndd_index: The index of the NDD
        vtx_indices: The indices of the vertices in the chain, in order
        score: the chain's score
    """

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
            
def find_chains(digraph, ndds, max_chain, edge_success_prob=1):
    """Generate all chains with up to max_chain edges."""

    def find_chains_recurse(vertices, score):
        chains.append(Chain(ndd_idx, vertices[:], score))
        if len(vertices) < max_chain:
            for e in digraph.vs[vertices[-1]].edges:
                if e.tgt.id not in vertices:
                    vertices.append(e.tgt.id)
                    find_chains_recurse(vertices, score+e.score*edge_success_prob**len(vertices))
                    del vertices[-1]
    chains = []
    if max_chain == 0:
        return chains
    for ndd_idx, ndd in enumerate(ndds):
        for e in ndd.edges:
            vertices = [e.target_v.id]
            find_chains_recurse(vertices, e.score*edge_success_prob)
    return chains

