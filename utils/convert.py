# This script is for converting from the PrefLib .wmd format.
# Caution: it may not be robust to small changes in the file format.

# PrefLib seems to have modified the files to a slightly different format in 2022,
# so I've tried to modify the script to accept either the old or the new format.

import sys

class WmdException(Exception):
    pass

if len(sys.argv) > 1 and sys.argv[1] in ["-h", "--help"]:
    print("Convert .wmd from standard input to combined .input and .ndds format.")
else:
    first_line = sys.stdin.readline()
    new_file_format = first_line.startswith('#')

    n_pairs = 0
    n_ndds = 0

    if new_file_format:
        while True:
            line = sys.stdin.readline()
            if "# NUMBER ALTERNATIVES" in line:
                n = int(line.split()[3])
            if "# NUMBER EDGES" in line:
                m = int(line.split()[3])
                break
    else:
        n, m = [int(x) for x in first_line.split(",")]

    for i in range(n):
        s = sys.stdin.readline()
        if "Pair" in s:
            if n_ndds > 0:
                raise WmdException("Didn't expect a pair to appear after NDDs")
            n_pairs += 1
        else:
            n_ndds += 1


    pair_edges = [[] for _ in range(n_pairs)]
    ndd_edges = [[] for _ in range(n_ndds)]

    for i in range(m):
        tokens = sys.stdin.readline().split(',')
        src = int(tokens[0])
        tgt = int(tokens[1])
        if new_file_format:
            wt = float(tokens[2])
        else:
            wt = int(tokens[2])
        if new_file_format:
            src -= 1
            tgt -= 1
        if tgt < n_pairs:   # Discard edges to NDDs
            if src < n_pairs:
                pair_edges[src].append((tgt, wt))
            else:
                ndd_edges[src-n_pairs].append((tgt, wt))

    def write(n_agents, edges):
        n_edges = len([e for l in edges for e in l])
        print("{}\t{}".format(n_agents, n_edges))
        for i in range(n_agents):
            for edge in edges[i]:
                print("{}\t{}\t{}".format(i, edge[0], edge[1]))
        print("-1\t-1\t-1")

    write(n_pairs, pair_edges)
    write(n_ndds, ndd_edges)
