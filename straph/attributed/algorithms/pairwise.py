from collections import defaultdict

############################################################
#           Pairwise Paths Algorithms for source nodes     #
############################################################

#######################################
#           F-Algorithm               #
#######################################

#####################################################################
#           Pairwise Paths Algorithms for temporal source nodes     #
#####################################################################

def F_Algorithm(S, L_functions):
    """
    An implementation of F-Algorithm (add ref).
    Pairwise Algorithm to compute temporal paths in Stream Graph.

    :param S: A Stream Graph
    :param L_functions: Functions according to the path's type (supported :
        - foremost path
        - shortest foremost path
        - shortest path
        - fastest shortest path
        - fastest path
        - shortest fastest path)
    :return:
    """
    # Initialisation
    F = {}
    for n, np in zip(S.nodes, S.node_presence):
        for t0, t1 in zip(np[::2], np[1::2]):
            s = (t0, t1, n)
            F[s] = L_functions['initialisation'](s)
    temporal_adjacency_list = defaultdict(set)  # currently present links in the form of adjacency list
    nodes_to_udpate = set()  # Reached nodes so far
    E = S.augmented_ordered_links()
    for i in E:
        c = i[0]
        if c == 1:  # Arrival
            _, t0, t1, u, v = i
            l = (t0, t1, u, v)
            temporal_adjacency_list[u].add((t1, v))
            temporal_adjacency_list[v].add((t1, u))
            nodes_to_udpate.add(u)
            nodes_to_udpate.add(v)
            F_update(F, [l], temporal_adjacency_list, nodes_to_udpate, L_functions)
        else:
            _, t1, u, v = i
            temporal_adjacency_list[u].remove((t1, v))
            temporal_adjacency_list[v].remove((t1, u))
    return F


def F_update(F, batch_arrival, temporal_adjacency_list, nodes_to_update, L_functions):
    """
    Proceed to BFS_Update W_n for every encounter node so far (i.e. in nodes_to_update).
    :param F:
    :param batch_arrival:
    :param temporal_adjacency_list:
    :param nodes_to_update:
    :param L_functions:
    :return:
    """
    for n in nodes_to_update:
        bfs_update(F[n], n, batch_arrival, temporal_adjacency_list, L_functions)

