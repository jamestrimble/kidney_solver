class Edge(object):
    def __init__(self, src, tgt, weight):
        self.src = src
        self.tgt = tgt
        self.weight = weight

def read_instance(lines):
    lines = lines[:]  # Create a working copy of input

    n_pairs, pair_edge_count = [int(x) for x in lines[0].split()]
    lines = lines[1:]

    pair_edges = []
    for line in lines[:pair_edge_count]:
        src, tgt, weight = line.split()
        pair_edges.append(Edge(int(src), int(tgt), float(weight)))
    lines = lines[pair_edge_count:]

    assert lines[0].split()[0] == "-1"
    lines = lines[1:]

    if len(lines) == 0:
        n_ndds = 0
        ndd_edge_count = 0
    else:
        n_ndds, ndd_edge_count = [int(x) for x in lines[0].split()]
        lines = lines[1:]

        ndd_edges = []
        for line in lines[:ndd_edge_count]:
            src, tgt, weight = line.split()
            ndd_edges.append(Edge(int(src) + n_pairs, int(tgt), float(weight)))
        lines = lines[ndd_edge_count:]

    assert len(lines) == 1
    assert lines[0].split()[0] == "-1"

    return n_pairs, n_ndds, pair_edges, ndd_edges
