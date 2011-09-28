from numpy import inf
import networkx as nx
def graph_from_array(arr):
    import networkx as nx
    g  = nx.Graph()
    #node_number = 0
    for i in range(arr.shape[0]-1):
        for j in range(arr.shape[1]-1):
            if arr[i,j] == 0:
                continue
            if arr[i,j+1]!=0:
                g.add_edge((i,j),(i,j+1))
            if arr[i+1,j]!=0:
                g.add_edge((i,j),(i+1,j))
    return g

def graph_from_coords(instance,coords):
# TODO possible speed up is to use start, stop parameters of tuple.index() to reduce search
    graph  = nx.Graph()
    tuple_of_coords = tuple(zip(coords[0],coords[1]))
    for idx in range(len(tuple_of_coords)-1): #-1 so i dont check the item after the last
        if tuple_of_coords[idx][1]+1 == tuple_of_coords[idx+1][1]:
            graph.add_edge(idx,idx+1,neightbour_in_same='row')
        try:
            graph.add_edge(idx,
                       tuple_of_coords.index((tuple_of_coords[idx][0]+1, tuple_of_coords[idx][1])),
                       neighbour_in_same='column')
        except ValueError:
            print 'No node below node: ',(idx)
    return graph, tuple_of_coords


def digraph_from_coords(instance,coords):
# TODO possible speed up is to use start, stop parameters of tuple.index() to reduce search
    graph  = nx.DiGraph()
    tuple_of_coords = tuple(zip(coords[0],coords[1]))
    for idx in range(len(tuple_of_coords)-1): #-1 so i dont check the item after the last
        if tuple_of_coords[idx][1]+1 == tuple_of_coords[idx+1][1]:
            graph.add_edge(idx,idx+1,weight=-instance.t0,neightbour_in_same='row')
            graph.add_edge(idx+1,idx,weight=-instance.t0,neightbour_in_same='row')
        try:
            graph.add_edge(idx,
                           tuple_of_coords.index((tuple_of_coords[idx][0]+1, tuple_of_coords[idx][1])),
                           weight=-instance.t0,
                           neighbour_in_same='column')
            graph.add_edge(tuple_of_coords.index((tuple_of_coords[idx][0]+1, tuple_of_coords[idx][1])),
                           idx,
                           weight=-instance.t0,
                           neighbour_in_same='column')
        except ValueError:
            print 'No node below node: ',(idx)
    return graph, tuple_of_coords
"""
Breadth First Search.
D. Eppstein, May 2007.
modified by
J. Siegl, Sep 2011.
"""

def BreadthFirstLevels(graph,root,locked_nodes=(),end=None,level=None,max_nodes=inf):
    # TODO speed up
    """
    Generate a sequence of bipartite directed graphs, each consisting
    of the edges from level i to level i+1 of G. Edges that connect
    vertices within the same level are not included in the output.
    The vertices in each level can be listed by iterating over each
    output graph.
    """
    visited = set(locked_nodes)
    currentLevel = set(root)
    nodes_distributed=0
    count = 0
    if end is not None:
        end = set(end)
    else:
        end = set()
    while currentLevel:
        #from pudb import set_trace; set_trace()
        if (end & currentLevel):
            rest_of_nodes = set(graph.nodes())-visited
            nodes_distributed += len(rest_of_nodes)
            yield rest_of_nodes
            break
        for v in currentLevel:                                  #combine in visited = set(currentLevel)
            visited.add(v)
        nextLevel = set()
        #levelGraph = dict([(v,set()) for v in currentLevel])    #not needed
        for v in currentLevel:
            for w in graph[v]:
                if w not in visited:
                    #levelG]ra2ph[v].add(w)                        #not needed
                    nextLevel.add(w)
        nodes_distributed += len(currentLevel)
        yield currentLevel
        if nodes_distributed >= max_nodes:break
        if count == level: break
        currentLevel = nextLevel
        count +=1

def colorarray_from_levelstructure(instance,levelstructure):
    from numpy import zeros
    colorarray = zeros(instance.wafer.shape)
    color = 1
    for level in levelstructure:
        for node in level:
            colorarray[instance.tuple_of_coords[node]] = color
        color+=1
    return colorarray

def bisect(graph,Ni,nodes_left,nodes_to_bisect,nodes_right,locked_nodes=set()):
    from matplotlib.cbook import flatten
    from numpy import floor
#return when at the end of the recursive bisection and no more levels to bisect
    if Ni == 1: return [set(sorted(nodes_to_bisect))]
#calculate how many levels are in any of the two parts. Caution: Integer division!
    Ni1 = Ni/2
    Ni2 = Ni-Ni/2
#lock all nodes except those to bisect
    locked_set = set(graph.nodes())-nodes_to_bisect
    nodes_i1_bfs = set(flatten(BreadthFirstLevels(graph,root=nodes_left,locked_nodes=locked_set,level=Ni1)))-nodes_left
    locked_set |= nodes_i1_bfs
    nodes_i2_bfs = set(flatten(BreadthFirstLevels(graph,root=nodes_right,locked_nodes=locked_set,level=Ni2)))-nodes_right
    locked_set |= nodes_i2_bfs
    max_nodes_i1 = floor(Ni1*len(nodes_to_bisect)/float(Ni))
    max_nodes_i2 = len(nodes_to_bisect)-max_nodes_i1
    if len(nodes_i1_bfs) <= len(nodes_i2_bfs):
        nodes_i1 = set(flatten(BreadthFirstLevels(graph,root=nodes_i1_bfs,locked_nodes=locked_set,max_nodes=max_nodes_i1)))
        nodes_i2 = nodes_i2_bfs | (nodes_to_bisect-nodes_i1)
    else:
        nodes_i2 = set(flatten(BreadthFirstLevels(graph,root=nodes_i2_bfs,locked_nodes=locked_set,max_nodes=max_nodes_i2)))
        nodes_i1 = nodes_i1_bfs | (nodes_to_bisect-nodes_i2)

    return bisect(graph,Ni1,nodes_left,nodes_i1,nodes_i2) + bisect(graph,Ni2,nodes_i1,nodes_i2,nodes_right)

def spingraph_from_graph(instance,graph):
    even_graph = nx.relabel_nodes(graph, lambda x:x*2)
    odd_graph = nx.relabel_nodes(graph, lambda x:2*x+1)
    union_graph  = nx.union(even_graph, odd_graph)
    for contact in instance.contacts:
        list_of_nodenames = list(contact.names)
        contact.names=set()
        for contact_node in list_of_nodenames:
            contact.names.add(contact_node*2)
            contact.names.add(2*contact_node+1)
    # from pudb import set_trace; set_trace()
    for spin_down_node in xrange(1,union_graph.order(),2):
        spin_up_node = spin_down_node -1
        for spin_down_node_neighbour in union_graph[spin_down_node].keys():
            if spin_down_node_neighbour % 2 ==0:
                continue
            if spin_down_node_neighbour < spin_down_node:             # is either top or left neighbour
                if spin_down_node_neighbour == spin_down_node-2:      # is left neighbour
                    union_graph.add_edge(spin_up_node,spin_down_node_neighbour,weight=-instance.tso)
                    union_graph.add_edge(spin_down_node_neighbour,spin_up_node,weight=-instance.tso)
                else:
                    union_graph.add_edge(spin_up_node,spin_down_node_neighbour,weight=+1j*instance.tso)
                    union_graph.add_edge(spin_down_node_neighbour,spin_up_node,weight=-1j*instance.tso)
            if spin_down_node_neighbour > spin_down_node:             # is either right or bottom neighbour
                if spin_down_node_neighbour == spin_down_node+2:      # is right neighbour
                    union_graph.add_edge(spin_up_node,spin_down_node_neighbour,weight=instance.tso)
                    union_graph.add_edge(spin_down_node_neighbour,spin_up_node,weight=instance.tso)
                else:
                    union_graph.add_edge(spin_up_node,spin_down_node_neighbour,weight=-1j*instance.tso)
                    union_graph.add_edge(spin_down_node_neighbour,spin_up_node,weight=+1j*instance.tso)
    instance.graph = union_graph
