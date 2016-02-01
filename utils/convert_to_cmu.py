import sys

class Edge(object):
    def __init__(self, src, tgt, weight):
        self.src = src
        self.tgt = tgt
        self.weight = weight

if len(sys.argv) > 1 and sys.argv[1] in ["-h", "--help"]:
    print "Convert concatenated .input and .ndds format from standard input to CMU .input format."
else:
    lines = sys.stdin.readlines()

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

    dummy_edges = []
    for i in range(n_pairs):
        for j in range(n_pairs, n_pairs+n_ndds):
            dummy_edges.append(Edge(i, j, 0))

    all_edges = pair_edges + ndd_edges + dummy_edges

    print "{}\t{}".format(n_pairs + n_ndds, len(all_edges))
    
    for edge in all_edges:
        # Print source, target, weight, (1 if target is an NDD else 0)
        print "{}\t{}\t{}\t{}\t0".format(edge.src, edge.tgt, edge.weight,
                                        1 if edge.tgt >= n_pairs else 0)

    print "-1\t-1\t-1"
