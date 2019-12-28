import sys
import instance_reader

if len(sys.argv) > 1 and sys.argv[1] in ["-h", "--help"]:
    print("Convert concatenated .input and .ndds format from standard input")
    print("to the matrix TSV format used by Kristiaan Glorie.") 
    print()
    print("Required parameters: cycle_cap chain_cap")
    print("Note that cycle_cap is number of *vertices* in longest permitted cycle.")

else:
    cycle_cap = sys.argv[1]
    chain_cap = sys.argv[2]

    lines = sys.stdin.readlines()

    n_pairs, n_ndds, pair_edges, ndd_edges = instance_reader.read_instance(lines)

    dummy_edges = []
    for i in range(n_pairs+n_ndds):
        for j in range(n_pairs, n_pairs+n_ndds):
            dummy_edges.append(instance_reader.Edge(i, j, 0))

    all_edges = pair_edges + ndd_edges + dummy_edges

    adj_mat = [[0]*(n_pairs+n_ndds) for i in range(n_pairs+n_ndds)]
    
    for edge in all_edges:
        adj_mat[edge.src][edge.tgt] = 1

    print(cycle_cap + "\t" + chain_cap)
    for i, row in enumerate(adj_mat):
        print("\t".join(str(x) for x in row) + "\t" + \
                ("P" if i < n_pairs else "A"))
