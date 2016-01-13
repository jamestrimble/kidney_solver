from collections import deque

def cycle_score(cycle, digraph):
    val = 0
    cycle_len = len(cycle)
    for i in range(cycle_len):
        v1 = cycle[i]
        v2 = cycle[(i+1) % cycle_len]
        edge = digraph.adj_mat[v1.id][v2.id]
        val += edge.score
    return val

class Vertex:
    def __init__(self, id, fail_prob):
        self.id = id
        self.fail_prob = fail_prob
        self.edges = []

    def _add_edge(self, edge):
        self.edges.append(edge)
        
    def __str__(self):
        return ("V" + str(self.id) + "  with edges to " +
                ", ".join([str(e.dest.id) for e in self.edges]))

class Edge:
    def __init__(self, id, fail_prob, score, src, dest):
        self.id = id
        self.fail_prob = fail_prob
        self.score = score
        self.src = src   # source vertex
        self.dest = dest # destination vertex

    def __str__(self):
        return ("V" + str(self.src.id) + "-V" + str(self.dest.id))

class Digraph:
    def __init__(self):
        self.vs = []
        self.es = []

    def add_vertex(self, fail_prob):
        id = len(self.vs)
        v = Vertex(id, fail_prob)
        self.vs.append(v)
        return id

    def add_edge(self, fail_prob, score, source, dest):
        id = len(self.es)
        v_s = self.vs[source.id]
        v_d = self.vs[dest.id]
        e = Edge(id, fail_prob, score, v_s, v_d)
        self.es.append(e)
        v_s.edges.append(e)
        return e
    
    def create_relabelled_copy(self, new_vtx_labels):
        old_vtx_labels = [None] * len(self.vs)
        for i in range(len(self.vs)):
            old_vtx_labels[new_vtx_labels[i]] = i
        new_digraph = Digraph()
        for i in range(len(self.vs)):
            old_vtx = self.vs[old_vtx_labels[i]]
            new_digraph.add_vertex(old_vtx.fail_prob)
        for e in self.es:
            new_src = new_digraph.vs[new_vtx_labels[e.src.id]]
            new_dest = new_digraph.vs[new_vtx_labels[e.dest.id]]
            new_digraph.add_edge(e.fail_prob, e.score, new_src, new_dest)
        return new_digraph

    def find_cycles(self, max_length):
        """Find cycles of length up to max_length in the digraph.

        Return value: a list of cycles, in which each cycle is represented
        as a list of vertices, with the first vertex not repeated.
        """

        cycles = []
        vtx_used = [False] * len(self.vs)

        def cycle(current_list):
            last_vtx = current_list[-1]
            if self.edge_exists(last_vtx, current_list[0]):
                cycles.append(current_list[:])
            if len(current_list) < max_length:
                for e in last_vtx.edges: 
                    v = e.dest
                    if (len(current_list) + shortest_paths_to_low_vtx[v.id] <= max_length
                                and not vtx_used[v.id]):
                        current_list.append(v)
                        vtx_used[v.id] = True
                        cycle(current_list)
                        vtx_used[v.id] = False
                        del current_list[-1]

        # Adjacency lists for transpose graph
        transp_adj_lists = [[] for v in self.vs]
        for edge in self.es:
            transp_adj_lists[edge.dest.id].append(edge.src)

        for v in self.vs:
            shortest_paths_to_low_vtx = self.calculate_shortest_path_lengths(
                    v, max_length - 1,
                    lambda u: (w for w in transp_adj_lists[u.id] if w.id > v.id))
            vtx_used[v.id] = True
            cycle([v])
            vtx_used[v.id] = False
        return cycles
    
    def get_shortest_path_from_low_vtx(self, low_vtx, max_path):
        """ Returns an array of path lengths. For each v > low_vtx, if the shortest
            path from low_vtx to v is shorter than max_path, then element v of the array
            will be the length of this shortest path. Otherwise, element v will be
            999999999."""
        return self.calculate_shortest_path_lengths(self.vs[low_vtx], max_path,
                    adj_list_accessor=lambda v: (e.dest for e in v.edges if e.dest.id >= low_vtx))

    def get_shortest_path_to_low_vtx(self, low_vtx, max_path):
        """ Returns an array of path lengths. For each v > low_vtx, if the shortest
            path to low_vtx from v is shorter than max_path, then element v of the array
            will be the length of this shortest path. Otherwise, element v will be
            999999999."""
        def adj_list_accessor(v):
            for i in range(low_vtx, len(self.vs)):
                if self.adj_mat[i][v.id]:
                    yield self.vs[i]
            
        return self.calculate_shortest_path_lengths(self.vs[low_vtx], max_path,
                    adj_list_accessor=adj_list_accessor)

    def calculate_shortest_path_lengths(self, from_v, max_dist,
                adj_list_accessor=lambda v: (e.dest for e in v.edges)):
        """Calculate the length of the shortest path from vertex from_v to each
        vertex with a greater or equal index, using paths containing
        only vertices indexed greater than or equal to from_v.

        Return value: a list of distances of length equal to the number of vertices.
        If the shortest path to a vertex is greater than max_dist, the list element
        will be 999999999.

        Arguments:
        from_v -- The starting vertex
        max_dist -- The maximum distance we're interested in
        adj_list_accessor -- An function taking a vertex and returning an
                iterable of out-edge targets
        """
        # Breadth-first search
        q = deque([from_v])
        distances = [999999999] * len(self.vs)
        distances[from_v.id] = 0

        while q:
            v = q.popleft()
            if distances[v.id] == max_dist:
                break
            for w in adj_list_accessor(v):
                if distances[w.id] == 999999999:
                    distances[w.id] = distances[v.id] + 1
                    q.append(w)

        return distances

    def edge_exists(self, v1, v2):
        return self.adj_mat[v1.id][v2.id] is not None
                    
    def create_adj_mat(self):
        n = len(self.vs)
        self.adj_mat = [[None for x in range(n)] for x in range(n)]
        for e in self.es:
            self.adj_mat[e.src.id][e.dest.id] = e

    def induced_subgraph(self, vertices):
        subgraph = Digraph()
        for v in vertices:
            subgraph.add_vertex(v.fail_prob)
#            subgraph.vs[-1].id_in_original_digraph = v.id
        for (i,v) in enumerate(vertices):
            for (j,w) in enumerate(vertices):
                e = self.adj_mat[v.id][w.id]
                if e is not None:
                    new_src = subgraph.vs[i]
                    new_dest = subgraph.vs[j]
                    subgraph.add_edge(e.fail_prob, e.score, new_src, new_dest)
        return subgraph

    def __str__(self):
        return "\n".join([str(v) for v in self.vs])
        
def read_digraph_without_prob(lines):
    DUMMY_PROB = 0
    digraph = Digraph()
    vtx_count, edge_count = [int(x) for x in lines[0].split()]
    for i in range(vtx_count):
        digraph.add_vertex(DUMMY_PROB)
    for line in lines[1:edge_count+1]:
        tokens = [x for x in line.split()]
        src_id = int(tokens[0])
        dest_id = int(tokens[1])
        score = float(tokens[2])
            
        digraph.add_edge(DUMMY_PROB, score, digraph.vs[src_id], digraph.vs[dest_id])
    return digraph   

