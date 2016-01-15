from kidney_solver.kidney_digraph import *
import kidney_solver.kidney_ip as k_ip
import kidney_solver.kidney_ndds as k_ndds
import kidney_solver.kidney_utils as k_utils

def read_with_ndds(basename):
    with open(basename + ".input") as f:
        lines = f.readlines()
    d = read_digraph(lines)
    with open(basename + ".ndds") as f:
        lines = f.readlines()
    ndds = k_ndds.read_ndds(lines, d)
    return d, ndds

def test_check_validity():
    pass #TODO
