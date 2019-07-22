import time

import nose.tools

from kidney_solver.kidney_digraph import *
from .simple_find_cycles import simple_find_cycles

def read(filename):
    with open(filename) as f:
        lines = f.readlines()
    return read_digraph(lines)

def test_cycle_score():
    d = Digraph(4)
    d.add_edge(1.5, d.vs[0], d.vs[1])
    d.add_edge(1, d.vs[1], d.vs[2])
    d.add_edge(1, d.vs[2], d.vs[3])
    d.add_edge(1, d.vs[3], d.vs[0])
    c = [d.vs[i] for i in range(4)]
    assert cycle_score(c, d) == 4.5

def test_shortest_path():
    d = read("test-fixtures/small1.input")
    spl = d.calculate_shortest_path_lengths(d.vs[1], 4)
    assert spl == [1, 0, 2, 2, 1, 2, 3, 4, 999999999]
    print(spl)

    min_vtx = 1
    spl = d.calculate_shortest_path_lengths(d.vs[min_vtx], 4,
            lambda v: (e.tgt for e in v.edges if e.tgt.id > min_vtx))
    assert spl == [999999999, 0, 999999999, 2, 1, 2, 3, 4, 999999999]
    print(spl)

    min_vtx = 7
    spl = d.calculate_shortest_path_lengths(d.vs[min_vtx], 4,
            lambda v: (e.tgt for e in v.edges if e.tgt.id > min_vtx))
    print(spl)
    assert spl == [999999999, 999999999, 999999999, 999999999, 999999999, 999999999, 999999999, 0, 1]

    min_vtx = 1
    spl = d.calculate_shortest_path_lengths(d.vs[min_vtx], 4,
            lambda v: [d.vs[v.id+1]] if v.id<4 else [])
    assert spl == [999999999, 0, 1, 2, 3, 999999999, 999999999, 999999999, 999999999]
    print(spl)


def test_find_cycles():
    d = read("test-fixtures/100.input")
    max_cycle = 5
    start = time.time()
    cycles = d.find_cycles(max_cycle)
    print(time.time() - start)
    print(len(cycles))
    start = time.time()
    slow_cycles = simple_find_cycles(d, max_cycle)
    print(time.time() - start)
    assert len(cycles) > 100
    assert len(slow_cycles) == len(cycles)

@nose.tools.raises(KidneyReadException)
def test_raises_exception_on_self_loop():
    d = read("test-fixtures/self-loop.input")

@nose.tools.raises(KidneyReadException)
def test_raises_exception_on_duplicate_edge():
    d = read("test-fixtures/duplicate-edge.input")

@nose.tools.raises(KidneyReadException)
def test_raises_exception_on_incorrect_edge_count_1():
    d = read("test-fixtures/incorrect_edge_count1.input")

@nose.tools.raises(KidneyReadException)
def test_raises_exception_on_incorrect_edge_count_2():
    d = read("test-fixtures/incorrect_edge_count2.input")

@nose.tools.raises(KidneyReadException)
def test_raises_exception_on_index_out_of_range_1():
    d = read("test-fixtures/out-of-range1.input")

@nose.tools.raises(KidneyReadException)
def test_raises_exception_on_index_out_of_range_2():
    d = read("test-fixtures/out-of-range2.input")

@nose.tools.raises(KidneyReadException)
def test_raises_exception_on_index_out_of_range_3():
    d = read("test-fixtures/out-of-range3.input")

@nose.tools.raises(KidneyReadException)
def test_raises_exception_on_index_out_of_range_4():
    d = read("test-fixtures/out-of-range4.input")
