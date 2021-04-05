# Copyright (C) 2017-2021 Léo Rannou - Sorbonne Université/LIP6 - Thales
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from joblib import Parallel, delayed

import straph as sg


########################################################
#           Weakly Connected Components : DFS          #
########################################################


def compute_wcc_dfs(S, free_memory=False):
    """
    Compute the Weakly Connected Components of a ``StreamGraph`` using a Depth First Search procedure.

    :param S:  A ``StreamGraph`` object.
    :param free_memory: Optional parameter to free some memore. WARNING: \
    It does impact the original ``StreamGraph`` object
    :return:
    """
    components = []
    Neighborhood = S.neighborhood_with_node_presence()
    if free_memory:
        S.nodes = []
        S.node_presence = []
        S.links = []
        S.link_presence = []

    unvisited = set(Neighborhood.keys())
    while len(unvisited) > 0:
        v = unvisited.pop()
        current_component, visited = sg.DFS_iterative(v, Neighborhood)
        unvisited -= visited
        components.append(current_component)

    return components


########################################################
#      Weakly Connected Components as Substreams       #
########################################################


def compute_wcc_as_substreams(S, n_jobs=-1):
    """
    Return the weakly connected components a ``StreamGraph``
    as substreams (a stream graph induced by the component/cluster)

    :param S: A ``StreamGraph`` object
    :param n_jobs: Number of cores to allocate for a parallel computation.
    :return: A list of ``StreamGraph`` objects, one for each WCC.
    """
    list_WCC = S.weakly_connected_components()
    # 1. attribuer à chaque wcc ses events.
    list_sub_events = [[] for _ in list_WCC]
    seg_node_to_wcc = {}
    for i, c in enumerate(list_WCC):
        for t0, t1, n in c:
            list_sub_events[i].append((2, t0, t1, n))
            list_sub_events[i].append((-2, t1, n))
            seg_node_to_wcc[(t0, t1, n)] = i

    node_to_current_wcc = {}
    for e in S.ordered_arrivals():
        c = e[0]
        if c == 1:
            _, t0, t1, u, v = e
            id_wcc = node_to_current_wcc[u]
            list_sub_events[id_wcc].append((1, t0, t1, u, v))
            list_sub_events[id_wcc].append((-1, t1, u, v))
        elif c == 2:
            _, t0, t1, u = e
            node_to_current_wcc[u] = seg_node_to_wcc[(t0, t1, u)]

    def para_sg_from_events(ev, j):
        return sg.stream_graph_from_events_list(ev, j)

    list_substreams = Parallel(n_jobs=n_jobs, mmap_mode='r+')(
        delayed(para_sg_from_events)(ev, i) for i, ev in enumerate(list_sub_events))
    return list_substreams


#############################################################
#           Weakly Connected Components : Streaming         #
#############################################################

def find_node(u, dict_components):
    # First step : Find the node's root (component)
    p = u
    while dict_components[p] != p:
        p = dict_components[p]
    # Second step : Update the childs according to the root
    v = u
    while dict_components[v] != v:
        tmp = dict_components[v]
        dict_components[v] = p
        v = tmp
    return p


def link_components(u, v, dict_components, rank):
    # Choose the "biggest component" to append the node
    if rank[u] > rank[v]:
        dict_components[v] = u
    else:
        dict_components[u] = v
        if rank[u] == rank[v] and u != v:
            rank[v] += 1


def compute_wcc_streaming(S, reformat=True, free_memory=False):
    """
    Compute the Weakly Connected Components of a ``StreamGraph`` in a streaming fashion with the Union-Find algorithm.

    :param S:  A ``StreamGraph`` object.
    :param reformat: If False, output a dictionary associating each root node to its child node (the other member \
    of its wcc). If True, output WCC as a list of clusters.
    :param free_memory: Optional parameter to free some memore. \
    WARNING: It does impact the original ``StreamGraph`` object
    :return: Depends on the 'reformat' parameter. By default, a list of clusters.
    """
    node_to_wcc = {}
    rank = {}
    node_to_segmented_node = {}
    E = S.ordered_arrivals(free_memory=free_memory)
    for e in E:
        c = e[0]
        if c == 2:
            _, t0, t1, n = e
            rank[(t0, t1, n)] = 0
            node_to_wcc[(t0, t1, n)] = (t0, t1, n)
            node_to_segmented_node[n] = (t0, t1, n)
        elif c == 1:
            _, t0, t1, u, v = e
            u, v = node_to_segmented_node[u], node_to_segmented_node[v]

            root_u = find_node(u, node_to_wcc)
            root_v = find_node(v, node_to_wcc)
            link_components(root_u, root_v,
                            node_to_wcc,
                            rank)
    if reformat:
        def reformat_components(dict_components):
            k = 0
            dict_roots = {}
            for w in dict_components:
                if dict_components[w] == w:
                    dict_roots[w] = k
                    k += 1
            components = [[] for _ in range(k)]
            for w in dict_components:
                p = w
                while dict_components[p] != p:
                    p = dict_components[p]
                components[dict_roots[p]].append(w)
            return components

        wcc = reformat_components(node_to_wcc)

        return wcc
    return node_to_wcc
