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
import matplotlib.pyplot as plt
import random
from collections import defaultdict

from straph.utils import get_cmap


class ConnectedComponent:
    def __init__(self,
                 id=None,
                 times=None,
                 nodes=None,
                 active_links=None,
                 links=None
                 ):
        """
        A basic constructor for a ``ConnectedComponent`` object

        :param id : A unique identifier for the component
        :param times : Time window of the component [begin time, end time]
        :param nodes : A set of nodes present in the component
        :param active_links : A set of links present in the component (Only useful during components computation)
        :param links : a list of 'segmented' links
        """
        self.id = id
        self.times = times
        self.nodes = nodes
        self.active_links = active_links
        self.links = links

    def __copy__(self):
        t = copy.copy(self.times)
        if self.links:
            l = [copy.copy(l) for l in self.links]
        else:
            l = None
        n = copy.copy(self.nodes)
        return ConnectedComponent(times=t,
                                  nodes=n,
                                  links=l)

    def set_id(self, id):
        self.id = id

    def clear(self):
        """
        Clear the component

        :return:
        """
        self.nodes.clear()
        self.links.clear()
        self.times.clear()
        self.active_links.clear()

    def size(self):
        """
        Return the size of the component, its number of nodes

        :return: The component size (an integer)
        """
        if self.nodes:
            return len(self.nodes)
        else:
            return len(self.get_nodes())

    def to_adjacency_list(self, ordering=None):
        """
        Return the adjacency list corresponding to the links present in the component.
        An assumption is made, the component must be equivalent to a static graph as we return the
        adjacency list of the aggregated component.

        :param ordering: Dictionary providing an ordering for the component's nodes, the associated adjacency list \
        will be directed. A directed edge between *u* and *v* will only be created if ordering[u] < ordering[v].
        :return: The component adjacency list (a dictionnary).
        """

        al = defaultdict(set)
        if self.links:
            if ordering:
                for l in self.links:
                    _, _, u, v = l
                    if ordering[v] < ordering[u]:
                        al[v].add(u)
                    elif ordering[u] < ordering[v]:
                        al[u].add(v)
            else:
                for l in self.links:
                    _, _, u, v = l
                    al[u].add(v)
                    al[v].add(u)
        return al

    def surface(self):
        """
        The component's surface is defined as the sum of its links presence.

        :return: The component's surface (a float).
        """
        return sum([t1 - t0 for t0, t1, _, _ in self.links])

    def set_begin_time(self, t):
        self.times = [t, t]
        if self.links:
            self.links = [l for l in self.links if l[1] >= t and (l[2], l[3]) in self.active_links]

    def set_end_time(self, t):
        self.times[1] = t

    def add_node(self, n):
        self.nodes.add(n)

    def add_link(self, link):
        u, v = link[2], link[3]
        if u not in self.nodes:
            self.add_node(u)
        if v not in self.nodes:
            self.add_node(v)
        self.active_links.add((u, v))
        if self.links:
            self.links.append(list(link))

    def merge(self, comp):
        """
        Merge the component *comp* into the current one.

        :param comp: Another ``ConnectedComponent`` object
        :return: Nothing, it modify inplace the current component.
        """
        self.nodes |= comp.nodes
        self.active_links |= comp.active_links
        if self.links:
            self.links += comp.links

    def remove_link(self, link):
        self.active_links.discard(link)

    def remove_node(self, n):
        self.nodes.discard(n)

    def get_nodes(self):
        """
        Retrieve the nodes in the component based on its set of present links (``active_links``)

        :return: A set of nodes
        """
        return {n for l in self.active_links for n in l}

    def split(self):
        """
        Split the component based on its ``active_links``.
        This method is only used in connected component algorithms.

        :return: A list containing the (nodes, links, set_links) of the connected components stemming \
        from the current component.
        """
        # CUSTOM BFS
        R = []
        component_2_set_links = []
        node_2_component = {}
        for l in self.active_links:
            n1, n2 = l
            if n1 not in node_2_component and n2 not in node_2_component:
                n_comp = len(component_2_set_links)
                node_2_component[n1] = n_comp
                node_2_component[n2] = n_comp
                component_2_set_links.append({l})

            elif n1 in node_2_component and n2 not in node_2_component:
                n_comp = node_2_component[n1]
                node_2_component[n2] = n_comp

                component_2_set_links[n_comp].add(l)

            elif n1 not in node_2_component and n2 in node_2_component:
                n_comp = node_2_component[n2]
                node_2_component[n1] = n_comp

                component_2_set_links[n_comp].add(l)

            elif node_2_component[n1] != node_2_component[n2]:
                n_comp_1, n_comp_2 = node_2_component[n1], node_2_component[n2]
                if len(component_2_set_links[n_comp_2]) > len(component_2_set_links[n_comp_1]):
                    # If a component is bigger than another we merge into the bigger one.
                    n_comp_1, n_comp_2 = n_comp_2, n_comp_1

                for e in component_2_set_links[n_comp_2]:
                    node_2_component[e[0]] = n_comp_1
                    node_2_component[e[1]] = n_comp_1

                component_2_set_links[n_comp_1] |= component_2_set_links[n_comp_2]
                component_2_set_links[n_comp_1].add(l)

                component_2_set_links[n_comp_2] = None

            else:
                component_2_set_links[node_2_component[n1]].add(l)

        if sum([1 for el in component_2_set_links if el is not None]) > 1:
            for set_links in component_2_set_links:
                if set_links:
                    if self.links:
                        current_links = [l for l in self.links if (l[2], l[3]) in set_links]
                    else:
                        current_links = None
                    current_nodes = {n for l in set_links for n in l}
                    R.append((current_nodes, current_links, set_links))
        return R

    def get_interactions_times(self):
        """
        Return the set of interactions times in the component.

        :return: A set.
        """
        interact_times = set()
        for t0, t1, _, _ in self.links:
            t0 = max(t0, self.times[0])
            t1 = min(t1, self.times[1])
            interact_times.add(t0)
            interact_times.add(t1)
        return sorted(interact_times)

    def adjacency_list_at_t(self, t):
        """
        Return the adjacency list of the induced instant graph at time *t* of the component.

        :param t: A time instantin the component's time window
        :return: A dictionary
        """

        a_l = defaultdict(set)
        for t0, t1, u, v in self.links:
            if t0 <= t <= t1:
                a_l[u].add(u)
                a_l[v].add(v)
        return a_l

    def degrees(self):
        """
        Return the nodes degree in the aggregated graph induced by the component.

        :return: A dictionary associating each node to its degree
        """

        a_l = self.to_adjacency_list()
        return {n: len(a_l[n]) for n in a_l}

    def random_path(self, source, destination, rand_time=False):
        """
        Return a random path between the source and the destination inside the component

        :param source: A source node
        :param destination: A destination node
        :param rand_time: Force the path to occur at a random instant in the component (optional).
        :return: A random path in the component
        """
        if rand_time:
            t = random.random() * (self.times[1] - self.times[0]) + self.times[0]  # Random time
        else:
            t = self.times[0]
        a_l = self.adjacency_list_at_t(t)
        path_queue = [(source, [source])]
        # BFS
        while path_queue:
            (v, path) = path_queue.pop(0)
            if destination in a_l[v]:
                return t, path + [destination]
            else:
                path_queue += [(w, path + [w]) for w in a_l[v]]
        return

    def random_path_ss(self, source):
        """
        Return random paths between the source and all the other nodes in the component.

        :param source: A souce node
        :return:
        """
        assert source in self.nodes
        t = random.random() * (self.times[1] - self.times[0]) + self.times[0]  # Random time
        a_l = self.adjacency_list_at_t(t)
        path_queue = [(source, [source])]
        targets = set(self.nodes)
        targets.discard(source)
        # BFS
        while targets:
            (v, path) = path_queue.pop(0)
            to_discard = set()
            for w in a_l[v]:
                if w in targets:
                    to_discard.add(w)
                    yield t, path + [w]
            targets -= to_discard
            path_queue += [(w, path + [w]) for w in a_l[v]]
        return

    def random_path_pw(self):
        """
        Return random paths between every node (pairwise)

        :return:
        """
        paths = []
        for n in self.nodes:
            paths.append(self.random_path_ss(n))
        return paths

    def plot(self, plot_links=True, title=None):
        """
        Display the ``ConnectedComponent`` using ``matplotlib``

        :param plot_links: To display links in the component
        :param title: Figure title
        :return:  A ``matplotlib`` figure
        """
        lnodes = len(self.nodes)
        c_map = get_cmap(lnodes)
        dict_colors = {n: c_map(i) for n, i in zip(self.nodes, range(lnodes))}
        fig = plt.figure()
        # Plot Clusters
        nodes_list = list(self.nodes)

        for p in self.nodes:
            coln = dict_colors[p]
            plt.hlines([nodes_list.index(p)], xmin=self.times[0], linewidth=2,
                       xmax=self.times[1], colors=coln, alpha=1)
        # Plot Links
        if plot_links is True:
            for t0, t1, u, v in self.links:
                id1 = nodes_list.index(u)
                id2 = nodes_list.index(v)
                idmax = max(id1, id2)
                idmin = min(id1, id2)
                eps = random.choice([1, -1]) * (random.random() / 5)
                plt.hlines([(idmax + idmin) / 2 + eps], xmin=[t0], xmax=[t1],
                           colors='k',
                           linewidth=1.7,
                           alpha=0.5)
                plt.vlines([t0],
                           ymin=idmin, ymax=idmax,
                           linewidth=1.4, alpha=0.15)
                # ARRIVALS
                plt.plot([t0], [idmin], color='#004d00', marker='^', alpha=1, markersize=7)
                plt.plot([t0], [idmax], color='#004d00', marker='v', alpha=1, markersize=7)
        plt.yticks(range(lnodes), self.nodes,
                   fontname='Ubuntu', fontsize=12, color='#666699')
        plt.xticks(range(int(self.times[0]), int(self.times[1]) + 1),
                   fontname='Ubuntu', fontsize=12, color='#476b6b')
        plt.ylabel("Nodes",
                   fontname='Ubuntu', fontsize=12, color='#666699')
        plt.xlabel("t",
                   fontname='Ubuntu', fontsize=12, color='#476b6b')
        if title:
            plt.title(title, fontname='Ubuntu', fontsize=14)
        else:
            plt.title("Stream Graph", fontname='Ubuntu', fontsize=14)
        # Get rid of the frame
        for spine in plt.gca().spines.values():
            spine.set_visible(False)
        plt.tick_params(top='off', bottom='on', right='off', left='on', labelleft='on', labelbottom='on')
        return fig

    def is_connected(self):
        """
        Return True if the component is connected, False otherwise.

        :return: A boolean
        """
        set_comp = {}
        node_to_id_comp = {}
        for l in self.links:
            _, _, u, v = l
            if u in node_to_id_comp and v in node_to_id_comp:
                id_u_comp = node_to_id_comp[u]
                id_v_comp = node_to_id_comp[v]
                if id_u_comp != id_v_comp:
                    # MERGE u comp into v comp
                    u_comp = set_comp[id_u_comp]
                    v_comp = set_comp[id_v_comp]
                    u_comp |= v_comp
                    for n in v_comp:
                        node_to_id_comp[n] = id_u_comp
                    set_comp[id_v_comp] = None
            elif u in node_to_id_comp and v not in node_to_id_comp:
                id_u_comp = node_to_id_comp[u]
                u_comp = set_comp[id_u_comp]
                node_to_id_comp[v] = id_u_comp
                u_comp.add(v)

            elif u not in node_to_id_comp and v in node_to_id_comp:
                id_v_comp = node_to_id_comp[v]
                v_comp = set_comp[id_v_comp]
                node_to_id_comp[u] = id_v_comp
                v_comp.add(u)
            else:
                new_comp = {u, v}
                node_to_id_comp[u], node_to_id_comp[v] = len(set_comp), len(set_comp)
                set_comp[len(set_comp)] = new_comp
        if sum([1 for c in set_comp if set_comp[c] is not None]) == 1:
            return True
        return False
