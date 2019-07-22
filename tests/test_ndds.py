import nose.tools

from kidney_solver.kidney_digraph import *
import kidney_solver.kidney_ndds as k_ndds
from .test_ip import read_with_ndds

@nose.tools.raises(KidneyReadException)
def test_raises_exception_on_duplicate_edge():
    d, ndds = read_with_ndds("test-fixtures/duplicate-ndd-edge")

def test_chain_finding():
    d, ndds = read_with_ndds("test-fixtures/chain-finding")
    chains = k_ndds.find_chains(d, ndds, 3)
    lengths = sorted(len(c.vtx_indices) for c in chains)
    scores = sorted(c.score for c in chains)
    assert lengths == [1, 2, 3]
    assert scores == [1.2, 2.7, 3.7]

    chains = k_ndds.find_chains(d, ndds, 1)
    lengths = sorted(len(c.vtx_indices) for c in chains)
    scores = sorted(c.score for c in chains)
    assert lengths == [1]
    assert scores == [1.2]

    chains = k_ndds.find_chains(d, ndds, 0)
    assert len(chains) == 0
