import sys
import instance_reader

if len(sys.argv) > 1 and sys.argv[1] in ["-h", "--help"]:
    print("Convert concatenated .input and .ndds format from standard input to CMU .input format.")
else:
    lines = sys.stdin.readlines()

    n_pairs, n_ndds, pair_edges, ndd_edges = instance_reader.read_instance(lines)

    dummy_edges = []
    for i in range(n_pairs):
        for j in range(n_pairs, n_pairs+n_ndds):
            dummy_edges.append(instance_reader.Edge(i, j, 0))

    all_edges = pair_edges + ndd_edges + dummy_edges

    print("{}\t{}".format(n_pairs + n_ndds, len(all_edges)))
    
    for edge in all_edges:
        # Print source, target, weight, (1 if target is an NDD else 0)
        print("{}\t{}\t{}\t{}\t0".format(edge.src, edge.tgt, edge.weight,
                                        1 if edge.tgt >= n_pairs else 0))

    print("-1\t-1\t-1")
