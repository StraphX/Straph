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

from straph import dags as cdag
from straph.components import ConnectedComponent


def compute_strongly_connected_components(S, format="object_with_links", condensation_dag=False, isolated_nodes=True,
                                          streaming_output=None, free_memory=False):
    """
    Compute Strongly Connected Components (SCC) of a ``StreamGraph``.

    :param free_memory:
    :param streaming_output:
    :param isolated_nodes:
    :param S: A Stream Graph
    :param format: Format of the output can be "cluster" or "scc_object"
    :param condensation_dag: Boolean, true if we want to output the Condensation DAG, false otherwise
    :return:
    """
    node_2_status = {}  # Dictionary associating a node to his current status : (current degree, number current comp)
    tmp_components = []  # List of current strongly connected components (object)
    final_components = []  # Clusters : [(t0,t1,u)]...
    cnt_scc_id = 0
    # Condensation DAG
    if condensation_dag:
        condensation_dag = cdag.CondensationDag()
    else:
        condensation_dag = None
    #
    id_wcc = S.id

    E = S.ordered_batch_links(free_memory=free_memory)

    if streaming_output:
        opt = open(streaming_output, 'w')
    else:
        opt = None

    for batch in E:
        # print("\n Batch :",batch)
        c = batch[0][0]
        if c == 1:  # ARRIVAL
            cnt_scc_id = process_batch_link_arrival(batch, node_2_status, tmp_components,
                                                    final_components,
                                                    cnt_scc_id,
                                                    condensation_dag=condensation_dag,
                                                    format=format,
                                                    streaming_output=opt)

        else:  # DEPARTURE
            cnt_scc_id = process_batch_link_departure(batch, node_2_status, tmp_components,
                                                      final_components,
                                                      cnt_scc_id,
                                                      condensation_dag=condensation_dag,
                                                      format=format,
                                                      streaming_output=opt)

    # Add isolated Nodes
    if isolated_nodes:
        for c in S.isolated_nodes():
            if format == "cluster":
                final_components.append([c])
                if condensation_dag:
                    condensation_dag.add_node([c])
            elif format == "object" or format == "object_with_links":
                c = StronglyConnectedComponent(id=cnt_scc_id,
                                               times=[c[0], c[1]],
                                               nodes={c[2]})
                final_components.append(c)
                if condensation_dag:
                    condensation_dag.add_node(c)
            elif format == "streaming":
                c = (c[0], c[1], 1)
                if streaming_output:
                    opt.write(str(c[0]) + ";" + str(c[1]) + ";" + str(1))
                    opt.write("\n")
                else:
                    final_components.append(c)
            cnt_scc_id += 1

    if condensation_dag:
        condensation_dag.set_id(id_wcc)
        return final_components, condensation_dag
    else:
        return final_components


def merge_scc(l, node_2_status, tmp_components,
              final_components,
              condensation_dag=None,
              cnt_scc_id=None,
              format="cluster", streaming_output=None):
    t0, t1, u, v = l
    id_comp_u = node_2_status[u][1]
    id_comp_v = node_2_status[v][1]

    if len(tmp_components[id_comp_v].nodes) > len(tmp_components[id_comp_u].nodes):
        #  If a component is bigger than another we merge into the bigger one.
        id_comp_u, id_comp_v = id_comp_v, id_comp_u

    comp_1 = tmp_components[id_comp_u]
    comp_2 = tmp_components[id_comp_v]

    if comp_1.times[0] != t0:
        cnt_scc_id = close_component(comp_1, t0, final_components, cnt_scc_id, condensation_dag=condensation_dag,
                                     format=format, streaming_output=streaming_output)
        # predecessor_in_dag[n_comp_1] += [cnt_scc_id - 1]  # previous closed component
        comp_1.set_begin_time(t0)

    if comp_2.times[0] != t0:
        cnt_scc_id = close_component(comp_2, t0, final_components, cnt_scc_id, condensation_dag=condensation_dag,
                                     format=format, streaming_output=streaming_output)
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
               condensation_dag=None,
               cnt_scc_id=None,
               format="cluster", streaming_output=None
               ):
    """

    :param streaming_output:
    :param format:
    :param final_components:
    :param node_to_update:
    :param node_to_add:
    :param l:
    :param node_2_status:
    :param tmp_components:
    :param condensation_dag:
    :param cnt_scc_id:
    :return:
    """
    t0, t1 = l[0], l[1]
    id_current_comp = node_2_status[node_to_update][1]
    current_comp = tmp_components[id_current_comp]
    if current_comp.times[0] != t0:
        cnt_scc_id = close_component(current_comp, t0, final_components, cnt_scc_id,
                                     condensation_dag=condensation_dag,
                                     format=format, streaming_output=streaming_output)
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
    tmp_components.append(StronglyConnectedComponent(times=[t0, t0],
                                                     nodes={u, v},
                                                     active_links={(u, v)},
                                                     links=lks))


def process_batch_link_arrival(batch, node_2_status, tmp_components, final_components,
                               cnt_scc_id, condensation_dag=None, format="cluster", streaming_output=None):
    for b in batch:
        t0, t1, u, v = b[1:]
        l = (t0, t1, u, v)
        if u not in node_2_status and v not in node_2_status:
            create_scc(l, node_2_status, tmp_components, format=format)

        elif u in node_2_status and v not in node_2_status:
            cnt_scc_id = update_scc(u, v, l, node_2_status, tmp_components,
                                    final_components,
                                    condensation_dag=condensation_dag,
                                    cnt_scc_id=cnt_scc_id,
                                    format=format, streaming_output=streaming_output
                                    )

        elif u not in node_2_status and v in node_2_status:
            cnt_scc_id = update_scc(v, u, l, node_2_status, tmp_components,
                                    final_components,
                                    condensation_dag=condensation_dag,
                                    cnt_scc_id=cnt_scc_id,
                                    format=format, streaming_output=streaming_output
                                    )

        elif node_2_status[u][1] != node_2_status[v][1]:
            cnt_scc_id = merge_scc(l, node_2_status, tmp_components,
                                   final_components,
                                   condensation_dag=condensation_dag,
                                   cnt_scc_id=cnt_scc_id,
                                   format=format, streaming_output=streaming_output
                                   )
        else:
            node_2_status[u][0] += 1
            node_2_status[v][0] += 1
            current_comp = tmp_components[node_2_status[u][1]]
            current_comp.add_link(l)

    return cnt_scc_id


def process_batch_link_departure(batch, node_2_status, tmp_components,
                                 final_components, cnt_scc_id, condensation_dag=None,
                                 format="cluster",
                                 streaming_output=None):
    """

    :param streaming_output:
    :param batch:
    :param node_2_status:
    :param tmp_components:
    :param final_components:
    :param condensation_dag:
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
                cnt_scc_id = close_component(comp, t1, final_components, cnt_scc_id, condensation_dag=condensation_dag,
                                             format=format, streaming_output=streaming_output)
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
                cnt_scc_id = close_component(comp, t1, final_components, cnt_scc_id, condensation_dag=condensation_dag,
                                             format=format, streaming_output=streaming_output)
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

    # Id comp to close is necessary.
    for id_comp in id_comp_to_close:
        comp = tmp_components[id_comp]
        cnt_scc_id = close_component(comp, t1, final_components, cnt_scc_id, condensation_dag=condensation_dag,
                                     format=format, streaming_output=streaming_output)
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
                    condensation_dag=None,
                    format="cluster",
                    streaming_output=None):
    """
    Close current component

    :param streaming_output:
    :param comp:
    :param t:
    :param final_components:
    :param cnt_scc_id:
    :param condensation_dag:
    :param format:
    :return:
    """
    c = None
    if format == "object" or format == "object_with_links":
        copy_comp = copy.copy(comp)
        copy_comp.set_end_time(t)  # Put an end time to the previous component
        copy_comp.id = cnt_scc_id
        final_components.append(copy_comp)
        c = copy_comp
    elif format == "cluster":
        c = [(comp.times[0], t, n) for n in comp.nodes]
        final_components.append(c)
    elif format == "streaming":
        n_nodes = len(comp.nodes)
        if streaming_output:
            streaming_output.write(str(comp.times[0]) + ";" + str(t) + ";" + str(n_nodes))
            streaming_output.write("\n")
        else:
            c = (comp.times[0], t, n_nodes)
            final_components.append(c)
    if condensation_dag:
        condensation_dag.add_node(c)
    cnt_scc_id += 1
    return cnt_scc_id


class StronglyConnectedComponent(ConnectedComponent):
    def __init__(self,
                 id=None,
                 times=None,
                 nodes=None,
                 active_links=None,
                 links=None
                 ):
        """
        A basic constructor for a connected component object

        :param id: an identifier (nojoke)
        :param times: [beginning time, ending time]
        :param nodes: A set of nodes present in the component
        :param active_links: A set of links present in the component (Only useful during components computation)
        :param links: a list of 'segmented' links
        """

        super().__init__(id, times, nodes, active_links, links)

    def __repr__(self):
        rep = "Id SCC :" + str(self.id) + " time window :" + str(self.times) + "\n"
        rep += "Nodes :" + str(self.nodes) + "\n"
        rep += "Links :" + str(self.links) + "\n"
        return rep

    def __copy__(self):
        comp_copy = StronglyConnectedComponent()
        comp_copy.times = copy.copy(self.times)
        if self.links:
            comp_copy.links = [copy.copy(l) for l in self.links]
        else:
            comp_copy.links = None
        comp_copy.nodes = copy.copy(self.nodes)
        return comp_copy

    def split(self):
        R = super().split()  # (current_nodes,current_links,active_links)
        return [StronglyConnectedComponent(nodes=c[0],
                                           links=c[1],
                                           active_links=c[2]) for c in R]

    def get_stable_components(self, format="object"):
        """
        Compute the stable connected components included in the current ``StronglyConnectedComponent``

        :param format:
        :return: A list of ``StableConnectedComponents`` objects
        """
        stable_components = []
        if self.links and len(self.get_interactions_times()) > 1:
            interact_times = self.get_interactions_times()
            time_to_pos = {t: i for t, i in zip(interact_times, range(len(interact_times)))}
            inter_nodes = [set() for _ in range(len(interact_times) - 1)]
            inter_links = [[] for _ in range(len(interact_times) - 1)]
            for l in self.links:
                t0, t1, u, v = l
                t0 = max(t0, self.times[0])
                t1 = min(t1, self.times[1])
                for i in range(time_to_pos[t0], time_to_pos[t1]):
                    inter_nodes[i].add(u)
                    inter_nodes[i].add(v)
                    inter_links[i].append((t0, t1, u, v))
            if format == "object" or format == "object_with_links":
                for j in range(len(interact_times) - 1):
                    c = StronglyConnectedComponent(id=(self.id, j),
                                                   times=(interact_times[j], interact_times[j + 1]),
                                                   nodes=set(
                                                       [u for u in inter_nodes[time_to_pos[interact_times[j]]]]),
                                                   links=[l for l in inter_links[time_to_pos[interact_times[j]]]]
                                                   )
                    stable_components.append(c)
            if format == "cluster":
                for j in range(len(interact_times) - 1):
                    c = [(interact_times[j], interact_times[j + 1], u) for u in
                         inter_nodes[time_to_pos[interact_times[j]]]]
                    stable_components.append(c)
        else:
            if format == "object" or format == "object_with_links":
                if self.links:
                    stable_components = [StronglyConnectedComponent(id=self.id,
                                                                    times=self.times,
                                                                    nodes=self.nodes,
                                                                    links=self.links)]
                else:
                    stable_components = [StronglyConnectedComponent(id=self.id,
                                                                    times=self.times,
                                                                    nodes=self.nodes)]
            if format == "cluster":
                stable_components = [[(self.times[0], self.times[1], u) for u in self.nodes]]
        return stable_components
