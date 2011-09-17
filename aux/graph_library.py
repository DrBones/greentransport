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
    import networkx as nx
    graph  = nx.Graph()
    tuple_of_coords = tuple(zip(coords[0],coords[1]))
    for idx in range(len(tuple_of_coords)-1): #-1 so i dont check the item after the last
        graph.add_edge(idx,idx,weight=4*instance.t0)
        if tuple_of_coords[idx][1]+1 == tuple_of_coords[idx+1][1]:
            graph.add_edge(idx,idx+1,neightbour_in_same='row')
        try:
            graph.add_edge(idx,
                       tuple_of_coords.index((tuple_of_coords[idx][0]+1, tuple_of_coords[idx][1])),
                       neighbour_in_same='column')
        except ValueError:
            print 'No node below node: ',(idx)
    return graph, tuple_of_coords

"""
Breadth First Search.
D. Eppstein, May 2007.
"""

def BreadthFirstLevels(G,root,end=None,level=None):
    # TODO speed up
    """
    Generate a sequence of bipartite directed graphs, each consisting
    of the edges from level i to level i+1 of G. Edges that connect
    vertices within the same level are not included in the output.
    The vertices in each level can be listed by iterating over each
    output graph.
    """
    visited = set()
    currentLevel = set(root)
    count = 0
    if end is not None:
        end = set(end)
    else:
        end = set()
    while currentLevel:
        #from pudb import set_trace; set_trace()
        for v in currentLevel:                                  #combine in visited = set(currentLevel)
            visited.add(v)
        nextLevel = set()
        #levelGraph = dict([(v,set()) for v in currentLevel])    #not needed
        for v in currentLevel:
            for w in G[v]:
                if w not in visited:
                    #levelGraph[v].add(w)                        #not needed
                    nextLevel.add(w)
        yield currentLevel
        if (end & nextLevel):
            yield  set(G.nodes())-visited
            break
        currentLevel = nextLevel
        if count == level: break
        count +=1

def colorarray_from_levelset(instance,levelset):
    from numpy import zeros
    colorarray = zeros(instance.wafer.shape)
    color = 0
    for i in levelset:
        for j in i:
            colorarray[instance.tuple_of_coords[j]] = color
        color+=1
    return colorarray

def bisect(graph,Ni,vleft,vright):
    if Ni == 1: return graph
    Ni1 = Ni/2
    Ni2 = Ni-Ni/2
    Gi1 = BreadthFirstLevels(graph,vleft,level=Ni1)
    Gi2 = BreadthFirstLevels(graph,vright,level=Ni2)


def blocktridiagonalize(instance):
    from matplotlib.cbook import flatten
    instance.generate_graph()
    #instance.add_contacts_to_graph()
    levelSet = BreadthFirstLevels(instance.graph,instance.contacts[0].names,instance.contacts[1].names)
    N = len(levelSet)
    Vleft  = instance.contacts[0].names
    Vright = instance.contacts[1].names
    G1 = instance.graph.subgraph(set(instance.graph.nodes())- Vleft - Vright)
    bisect(G1,N-2,Vleft,Vright)
    return permuted_nodelist


