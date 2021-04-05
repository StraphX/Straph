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

import copy
from collections import defaultdict

from straph import dags as sdag
from straph.components import ConnectedComponent


def compute_stable_connected_components(S, format="object_with_links", stable_dag=False, isolated_nodes=True,
                                        streaming_output=None, free_memory=False, node_list=None):
    """
    Compute the Stable Connected Components (StCC) of a Stream Graph.

    :param node_list:
    :param free_memory:
    :param streaming_output:
    :param isolated_nodes:
    :param S: A ``StreamGraph`` object
    :param format: Format of the output can be "cluster" or "object" or "object_with_links"
    :param stable_dag: Boolean, true if we want to output the Condensation DAG, false otherwise
    :return:
    """
    node_2_status = {}  # Dictionary associating a node to his current status : (current degree, number current comp)
    tmp_components = []  # List of current strongly connected components (object)
    final_components = []  # Clusters : [(t0,t1,u)]...
    cnt_scc_id = 0
    if stable_dag:
        stable_dag = sdag.StableDag()
    else:
        stable_dag = None
    #
    id_wcc = S.id

    E = S.ordered_batch_links(free_memory=free_memory)

    if streaming_output:
        opt = open(streaming_output, 'w')
    else:
        opt = None

    for batch in E:
        c = batch[0][0]
        if c == 1:  # ARRIVAL
            cnt_scc_id = process_batch_link_arrival(batch, node_2_status, tmp_components,
                                                    final_components,
                                                    cnt_scc_id,
                                                    stable_dag=stable_dag,
                                                    format=format,
                                                    streaming_output=opt,
                                                    node_list=node_list)

        else:  # DEPARTURE

            cnt_scc_id = process_batch_link_departure(batch, node_2_status, tmp_components,
                                                      final_components,
                                                      cnt_scc_id,
                                                      stable_dag=stable_dag,
                                                      format=format,
                                                      streaming_output=opt,
                                                      node_list=node_list)
    # Add isolated Nodes
    if isolated_nodes:
        for c in S.isolated_nodes(node_list=node_list):
            if format == "cluster":
                final_components.append([c])
                if stable_dag:
                    stable_dag.add_node([c])
            elif format == "object" or format == "object_with_links":
                c = StableConnectedComponent(id=cnt_scc_id,
                                             times=[c[0], c[1]],
                                             nodes={c[2]})
                final_components.append(c)
                if stable_dag:
                    stable_dag.add_node(c)
            elif format == "streaming":
                c = (c[0], c[1], 1)
                if streaming_output:
                    opt.write(str(c[0]) + ";" + str(c[1]) + ";" + str(1))
                    opt.write("\n")
                else:
                    final_components.append(c)
            cnt_scc_id += 1

    if stable_dag:
        stable_dag.set_id(id_wcc)
        return final_components, stable_dag
    else:
        return final_components


def merge_scc(l, node_2_status, tmp_components,
              final_components,
              stable_dag=None,
              cnt_scc_id=None,
              format="cluster", streaming_output=None,
              node_list=None):
    t0, t1, u, v = l
    id_comp_u = node_2_status[u][1]
    id_comp_v = node_2_status[v][1]

    if len(tmp_components[id_comp_v].nodes) > len(tmp_components[id_comp_u].nodes):
        #  If a component is bigger than another we merge into the bigger one.
        id_comp_u, id_comp_v = id_comp_v, id_comp_u

    comp_1 = tmp_components[id_comp_u]
    comp_2 = tmp_components[id_comp_v]

    if comp_1.times[0] != t0:
        cnt_scc_id = close_component(comp_1, t0, final_components, cnt_scc_id, stable_dag=stable_dag,
                                     format=format, streaming_output=streaming_output,
                                     node_list=node_list)
        # predecessor_in_dag[n_comp_1] += [cnt_scc_id - 1]  # previous closed component
        comp_1.set_begin_time(t0)

    if comp_2.times[0] != t0:
        cnt_scc_id = close_component(comp_2, t0, final_components, cnt_scc_id, stable_dag=stable_dag,
                                     format=format, streaming_output=streaming_output,
                                     node_list=node_list)
        # predecessor_in_dag[n_comp_1] += [cnt_scc_id - 1]  # previous closed component (n_comp_1 because we merge)

    for n in comp_2.nodes:
        node_2_status[n][1] = id_comp_u  # Actualize referencement before deletion of 2nd comp

    comp_1.merge(comp_2)
    comp_1.add_link(l)  # Add the current link
    node_2_status[u][0] += 1
    node_2_status[v][0] += 1

    tmp_components[id_comp_v] = None
    return cnt_scc_id


def update_scc(node_to_update, node_to_add, l, node_2_status, tmp_components,
               final_components,
               stable_dag=None,
               cnt_scc_id=None,
               format="cluster", streaming_output=None,
               node_list=None
               ):
    """

    :param node_list:
    :param streaming_output:
    :param format:
    :param final_components:
    :param node_to_update:
    :param node_to_add:
    :param l:
    :param node_2_status:
    :param tmp_components:
    :param stable_dag:
    :param cnt_scc_id:
    :return:
    """
    t0, t1 = l[0], l[1]
    id_current_comp = node_2_status[node_to_update][1]
    current_comp = tmp_components[id_current_comp]
    if current_comp.times[0] != t0:
        cnt_scc_id = close_component(current_comp, t0, final_components, cnt_scc_id,
                                     stable_dag=stable_dag,
                                     format=format, streaming_output=streaming_output,
                                     node_list=node_list)
        # predecessor_in_dag[n_current_comp] += [cnt_scc_id - 1]  # previous closed component
        current_comp.set_begin_time(t0)  # Input a new begining time
    current_comp.add_node(node_to_add)  # Add the node to the comp
    current_comp.add_link(l)  # Actualize the component with the new link
    node_2_status[node_to_add] = [1, id_current_comp]
    node_2_status[node_to_update][0] += 1
    return cnt_scc_id


def create_scc(l, node_2_status, tmp_components, format="cluster"):
    """
    Create a Strongly Connected Component from the link *l*

    :param format:
    :param l:
    :param node_2_status:
    :param tmp_components:
    :return:
    """
    t0, t1, u, v = l
    new_id_comp = len(tmp_components)
    node_2_status[u] = [1, new_id_comp]
    node_2_status[v] = [1, new_id_comp]
    if format == "object_with_links":
        lks = [[t0, t1, u, v]]
    else:
        lks = None
    tmp_components.append(StableConnectedComponent(times=[t0, t0],
                                                   nodes={u, v},
                                                   active_links={(u, v)},
                                                   links=lks))


def process_batch_link_arrival(batch, node_2_status, tmp_components, final_components,
                               cnt_scc_id, stable_dag=None, format="cluster", streaming_output=None,
                               node_list=None):
    for b in batch:
        t0, t1, u, v = b[1:]
        l = (t0, t1, u, v)
        if u not in node_2_status and v not in node_2_status:
            create_scc(l, node_2_status, tmp_components, format=format)

        elif u in node_2_status and v not in node_2_status:
            cnt_scc_id = update_scc(u, v, l, node_2_status, tmp_components,
                                    final_components,
                                    stable_dag=stable_dag,
                                    cnt_scc_id=cnt_scc_id,
                                    format=format, streaming_output=streaming_output,
                                    node_list=node_list
                                    )

        elif u not in node_2_status and v in node_2_status:
            cnt_scc_id = update_scc(v, u, l, node_2_status, tmp_components,
                                    final_components,
                                    stable_dag=stable_dag,
                                    cnt_scc_id=cnt_scc_id,
                                    format=format, streaming_output=streaming_output,
                                    node_list=node_list
                                    )

        elif node_2_status[u][1] != node_2_status[v][1]:
            cnt_scc_id = merge_scc(l, node_2_status, tmp_components,
                                   final_components,
                                   stable_dag=stable_dag,
                                   cnt_scc_id=cnt_scc_id,
                                   format=format, streaming_output=streaming_output,
                                   node_list=node_list
                                   )
        else:
            # We close the component and add the link to the new one :)
            t0, t1, u, v = l
            id_comp_uv = node_2_status[u][1]

            comp_uv = tmp_components[id_comp_uv]
            if comp_uv.times[0] != t0:
                cnt_scc_id = close_component(comp_uv, t0, final_components, cnt_scc_id, stable_dag=stable_dag,
                                             format=format, streaming_output=streaming_output,
                                             node_list=node_list)
                # predecessor_in_dag[n_comp_1] += [cnt_scc_id - 1]  # previous closed component
                comp_uv.set_begin_time(t0)

            comp_uv.add_link(l)  # Add the current link
            node_2_status[u][0] += 1
            node_2_status[v][0] += 1

    return cnt_scc_id


def process_batch_link_departure(batch, node_2_status, tmp_components,
                                 final_components, cnt_scc_id, stable_dag=None,
                                 format="cluster",
                                 streaming_output=None,
                                 node_list=None):
    """

    :param node_list:
    :param streaming_output:
    :param batch:
    :param node_2_status:
    :param tmp_components:
    :param final_components:
    :param stable_dag:
    :param format:
    :param cnt_scc_id:
    :return:
    """
    id_comp_to_split = set()
    id_comp_to_close = set()
    nodes_to_remove = set()
    t1 = batch[0][1]
    for l in batch:
        u, v = l[2], l[3]
        node_2_status[u][0] -= 1
        node_2_status[v][0] -= 1
        id_comp = node_2_status[u][1]
        comp = tmp_components[id_comp]
        comp.remove_link((u, v))
        # By default we split the component
        if node_2_status[u][0] == 0 or node_2_status[v][0] == 0:
            # If it's a node's departure, there is several cases:
            # 1. No more links in the components (it's empty)
            if not comp.active_links:
                cnt_scc_id = close_component(comp, t1, final_components, cnt_scc_id, stable_dag=stable_dag,
                                             format=format, streaming_output=streaming_output,
                                             node_list=node_list)
                tmp_components[id_comp] = None
                id_comp_to_split.discard(id_comp)
                id_comp_to_close.discard(id_comp)
                del node_2_status[u]
                del node_2_status[v]
            # 2. A node left but there is still some nodes inside (and other departure to come)
            else:
                if node_2_status[u][0] == 0:
                    id_comp_to_close.add(id_comp)
                    nodes_to_remove.add(u)
                    del node_2_status[u]
                if node_2_status[v][0] == 0:
                    id_comp_to_close.add(id_comp)
                    nodes_to_remove.add(v)
                    del node_2_status[v]
        else:
            id_comp_to_split.add(id_comp)

    for id_comp in id_comp_to_split:
        comp = tmp_components[id_comp]
        if comp.active_links:
            R = comp.split()
            if R:
                # We close the current component :)
                cnt_scc_id = close_component(comp, t1, final_components, cnt_scc_id, stable_dag=stable_dag,
                                             format=format, streaming_output=streaming_output,
                                             node_list=node_list)
                tmp_components[id_comp] = None
                id_comp_to_close.discard(id_comp)
                for C in R:
                    # New components
                    # assert is_connected(C)
                    C.set_begin_time(t1)  # set new begin time
                    new_id_comp = len(tmp_components)
                    for n in C.nodes:
                        node_2_status[n][1] = new_id_comp
                        # predecessor_in_dag_tmp[len(tmp_components)] += [cnt_scc_id - 1]  # previous closed component
                    tmp_components.append(C)  # to the antecedent of news comp
            else:
                # No Split, but we close it anyway :)
                id_comp_to_close.add(id_comp)

    # Id comp to close is necessary.
    for id_comp in id_comp_to_close:
        comp = tmp_components[id_comp]
        cnt_scc_id = close_component(comp, t1, final_components, cnt_scc_id, stable_dag=stable_dag,
                                     format=format, streaming_output=streaming_output,
                                     node_list=node_list)
        comp.nodes -= nodes_to_remove
        if comp.nodes:
            comp.set_begin_time(t1)  # A node left but other are still presents.
        else:
            raise ValueError("Starfoullah")

    return cnt_scc_id


def close_component(comp,
                    t,
                    final_components,
                    cnt_scc_id,
                    stable_dag=None,
                    format="cluster",
                    streaming_output=None,
                    node_list=None):
    """
    Close current component

    :param streaming_output:
    :param node_list:
    :param comp:
    :param t:
    :param final_components:
    :param cnt_scc_id:
    :param stable_dag:
    :param format:
    :return:
    """
    c = None
    if format == "object" or format == "object_with_links":
        copy_comp = copy.copy(comp)
        copy_comp.set_end_time(t)  # Put an end time to the previous component
        copy_comp.id = cnt_scc_id
        if node_list is None:
            final_components.append(copy_comp)
        elif set(node_list) & set(copy_comp.nodes):
            final_components.append(copy_comp)
        c = copy_comp
    elif format == "cluster":

        if node_list is None:
            c = [(comp.times[0], t, n) for n in comp.nodes]
            final_components.append(c)
        elif set(node_list) & set(comp.nodes):
            c = [(comp.times[0], t, n) for n in comp.nodes]
            final_components.append(c)
    elif format == "streaming":
        n_nodes = len(comp.nodes)
        if node_list is None:
            if streaming_output:
                streaming_output.write(str(comp.times[0]) + ";" + str(t) + ";" + str(n_nodes))
                streaming_output.write("\n")
            else:
                c = (comp.times[0], t, n_nodes)
                final_components.append(c)
        elif set(node_list) & set(comp.nodes):
            if streaming_output:
                streaming_output.write(str(comp.times[0]) + ";" + str(t) + ";" + str(n_nodes))
                streaming_output.write("\n")
            else:
                c = (comp.times[0], t, n_nodes)
                final_components.append(c)
    if stable_dag:
        if node_list is None:
            stable_dag.add_node(c)
        elif set(node_list) & set(c.nodes):
            stable_dag.add_node(c)
    cnt_scc_id += 1
    return cnt_scc_id


#############################
#   General Functions       #
#############################

def algo_kcores_batagelj(a_l, degrees, core_ordering=False):
    """
    Compute k_cores of a static graph from its adjacency list and nodes degrees
    References
    ----------
    [1] An O(m) Algorithm for Cores Decomposition of Networks
    Vladimir Batagelj and Matjaz Zaversnik, 2003.
    http://arxiv.org/abs/cs.DS/0310049

    :param a_l:
    :param degrees:
    :param core_ordering:
    :return: The core number of each node in the graph
    """
    ordering = None
    if core_ordering:
        ordering = {}
    bin = [0]
    sorted_nodes = sorted(degrees, key=degrees.get)
    curr_degree = 0
    for i, v in enumerate(sorted_nodes):
        if degrees[v] > curr_degree:
            bin.extend([i] * (degrees[v] - curr_degree))
            curr_degree = degrees[v]
    node_pos = {v: pos for pos, v in enumerate(sorted_nodes)}
    cores = degrees
    for v in sorted_nodes:
        for w in a_l[v]:
            if cores[w] > cores[v]:
                a_l[w].discard(v)
                pos = node_pos[w]
                bin_start = bin[cores[w]]
                node_pos[w] = bin_start
                node_pos[sorted_nodes[bin_start]] = pos
                sorted_nodes[bin_start], sorted_nodes[pos] = sorted_nodes[pos], sorted_nodes[bin_start]
                bin[cores[w]] += 1
                cores[w] -= 1
        if core_ordering:
            ordering[v] = len(ordering)
    if core_ordering:
        return cores, ordering
    return cores


def get_graph_from_ordering(t0, t1, links, ordering):
    a_l = defaultdict(set)
    nodes = set()
    for l in links:
        u = (t0, t1, l[0])
        v = (t0, t1, l[1])
        nodes.add(u)
        nodes.add(v)
        if ordering[v] < ordering[u]:
            a_l[v].add(u)
        elif ordering[u] < ordering[v]:
            a_l[u].add(v)
    return a_l


def algo_kcliques_KCList(k, a_l, node_label, C=None, R=None):
    """
    Compute k_cliques of a static graph from its core odering dag.
    References
    ---------
    [2] Listing k-cliques in Sparse Real-World Graphs
    Maximilien Danisch, et al., 2018
    https://dl.acm.org/citation.cfm?id=3186125

    :param k: The parameter k, number of nodes in the considered cliques
    :param a_l: Adjacency list of the core ordering dag
    :param node_label: label of each node
    :param C: List of current cliques
    :param R: List of completed cliques
    :return: The k_cliques of the graph.
    """
    if R is None:
        R = []
    if C is None:
        C = []
    if k == 2:
        for u in node_label:
            if a_l[u]:
                for v in a_l[u]:
                    if v in node_label:
                        if node_label[v] == k:
                            final_C = copy.copy(C) + [u, v]
                            R.append(final_C)

    else:
        for u in node_label:
            # INDUCED NODES
            new_node_label = defaultdict(int)
            if a_l[u]:
                for v in a_l[u]:
                    if v in node_label:
                        if node_label[v] == k:
                            new_node_label[v] = k - 1
                new_C = copy.copy(C) + [u]
                if new_node_label:
                    algo_kcliques_KCList(k - 1, a_l, new_node_label, C=new_C, R=R)
    return R


def neighborhood_and_degrees_from_links(t0, t1, links):
    a_l = defaultdict(set)
    for l in links:
        u = (t0, t1, l[0])
        v = (t0, t1, l[1])
        a_l[u].add(v)
        a_l[v].add(u)
    degrees = {n: len(a_l[n]) for n in a_l}
    return a_l, degrees


class StableConnectedComponent(ConnectedComponent):
    def __init__(self,
                 id=None,
                 times=None,
                 nodes=None,
                 active_links=None,
                 links=None,
                 clusters=None
                 ):
        """
        A basic constructor for a ``StableConnectedComponent`` object

        :param id: an identifier
        :param times: [beginning time, ending time]
        :param nodes: A set of nodes present in the component
        :param active_links: A set of links present in the component (Only useful during construction)
        :param links: a list of 'segmented' links
        :param clusters: a nodes partition of the set of *nodes* of the stable connected component (e.g. a community)
        """
        super().__init__(id, times, nodes, active_links, links)
        self.clusters = clusters

    def __repr__(self):
        rep = "Id StCC :" + str(self.id) + " time window :" + str(self.times) + "\n"
        rep += "Nodes :" + str(self.nodes) + "\n"
        rep += "Links :" + str(self.links) + "\n"
        return rep

    def __copy__(self):
        comp_copy = StableConnectedComponent()
        comp_copy.times = copy.copy(self.times)
        if self.links:
            comp_copy.links = [copy.copy(l) for l in self.links]
        else:
            comp_copy.links = None
        comp_copy.nodes = copy.copy(self.nodes)
        if self.clusters is not None:
            comp_copy.clusters = [copy.copy(c) for c in self.clusters]
        return comp_copy

    def split(self):
        R = super().split()  # [(current_nodes,current_links,set_links)]
        return [StableConnectedComponent(nodes=c[0],
                                         links=c[1],
                                         active_links=c[2]) for c in R]

    def core_number(self):
        # Compute cores with Batagelj
        L = defaultdict(list)
        t0, t1 = self.times

        if self.size() == 2:
            u, v = self.nodes
            L[1] = [(t0, t1, u), (t0, t1, v)]
            return L

        a_l = self.to_adjacency_list()
        degrees = {n: len(a_l[n]) for n in a_l}
        cores = algo_kcores_batagelj(a_l, degrees)  # Algo Batelgej
        for k, v in cores.items():
            L[v].append((t0, t1, k))
        return L

    def k_core(self, k):
        L = self.core_number()
        return L[k]

    def k_clique(self, k):
        # Compute cliques with KCList
        L = []
        t0, t1 = self.times
        if self.size() == 2:
            return L
        a_l = self.to_adjacency_list()
        degrees = {n: len(a_l[n]) for n in a_l}

        cores, core_ordering = algo_kcores_batagelj(a_l, degrees, core_ordering=True)

        a_l = self.to_adjacency_list(ordering=core_ordering)

        node_label = defaultdict(int, {n: k for n in degrees})
        cliques = algo_kcliques_KCList(k, a_l, node_label, R=[])
        for c in cliques:
            L.append([(t0, t1, u) for u in c])
        return L

    def all_cliques(self):
        # Compute cliques with KCList
        L = defaultdict(list)
        t0, t1 = self.times
        if self.size() == 2:
            return L
        a_l = self.to_adjacency_list()
        degrees = {n: len(a_l[n]) for n in a_l}
        cores, core_ordering = algo_kcores_batagelj(a_l, degrees, core_ordering=True)
        a_l = self.to_adjacency_list(ordering=core_ordering)
        max_core_number = max(cores.values())
        K = 3

        while K <= max_core_number + 1:
            node_label = defaultdict(int, {n: K for n in degrees})
            cliques = algo_kcliques_KCList(K, a_l, node_label, R=[])
            for c in cliques:
                L[K].append([(t0, t1, u) for u in c])
            K += 1
        return L
