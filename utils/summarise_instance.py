"""Read an instance in .input (optionally with .ndds) format from standard
input, and write a summary of vertex and arc counts to standard output.
"""

import sys

def summarise(lines):
    count_lines = [line.split() for line in lines if len(line.split())==2]
    if len(count_lines) > 2:
        print("Error: two many lines with two tokens", file=sys.stderr)
        sys.exit(1)
    elif len(count_lines) == 0:
        print("Error: not enough lines with two tokens", file=sys.stderr)
        sys.exit(1)
    elif len(count_lines) == 2:
        n_ndds, n_arcs_from_ndd = (int(x) for x in count_lines[1])
    else:
        n_ndds, n_arcs_from_ndd = 0, 0

    n_pairs, n_arcs_between_pairs = (int(x) for x in count_lines[0])

    print("n_pairs,n_arcs_between_pairs,n_ndds,n_arcs_from_ndd,n_vertices,n_arcs")
    print("{},{},{},{},{},{}".format(n_pairs, n_arcs_between_pairs,
            n_ndds, n_arcs_from_ndd,
            n_pairs + n_ndds,
            n_arcs_between_pairs + n_arcs_from_ndd))

if __name__=="__main__":
    summarise(sys.stdin.readlines())
