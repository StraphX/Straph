# Copyright (C) 2017-2020 Léo Rannou - Sorbonne Université/LIP6 - Thales
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

import math
import matplotlib.collections as mcol
import matplotlib.patches as mpatch
import matplotlib.pyplot as plt
import networkx as nx
import numpy
import random
from collections import defaultdict

from straph.utils import get_cmap


def merge_dags(dag_list):
    """
    Merge every other streams into the first one, we assume there's no connection between them
    and that links aren't constructed yet

    :param dag_list: A list of ``Dag`` objects
    :return:
    """
    scc1 = dag_list[0]
    min_t, max_t = None, None
    if scc1.times:
        min_t, max_t = scc1.times[0], scc1.times[1]
    for scc in dag_list[1:]:
        scc1.c_nodes += scc.c_nodes
        if scc.times:
            min_t = min(min_t, scc.times[0])
            max_t = max(max_t, scc.times[1])
    if scc1.times:
        scc1.times = [min_t, max_t]
    return scc1


def add_arrow(e, ax, direction='right', size=25, color=None, alpha=0.5):
    """
    Thanks to : https://stackoverflow.com/questions/34017866/arrow-on-a-line-plot-with-matplotlib
    add an arrow to a line.

    :param e: ((a_x,_y),(b_x,b_y))
    :param ax: matplotlib axes
    :param direction: 'left' or 'right'
    :param size: size of the arrow in fontsize points
    :param color: color of arrow (should be coherent with line color, or not )
    :param alpha:
    :return:
    """
    if direction == 'left':
        ax.annotate('',
                    xy=(e[0][0], e[0][1]), xycoords='data',
                    xytext=(e[1][0], e[1][1]), textcoords='data',
                    arrowprops=dict(arrowstyle="<-", connectionstyle="arc3", color=color, alpha=alpha),
                    size=size,
                    )
    else:
        ax.annotate('',
                    xy=(e[1][0], e[1][1]), xycoords='data',
                    xytext=(e[0][0], e[0][1]), textcoords='data',
                    arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color=color, alpha=alpha),
                    size=size,
                    )


class Dag:
    def __init__(self,
                 id=None,
                 times=None,
                 c_nodes=None,
                 c_links=None,
                 id_comp_to_comp=None,
                 node_to_id_comp=None,
                 segmented_node_to_id_comp=None,
                 adj_list=None,
                 ):
        """
        A basic constructor for a ``Dag`` object.

        :param id:
        :param times:
        :param c_nodes : A list of nodes identifiers (each node represents a component:
        a set of nodes, a begin time, an end time)
        :param c_links : A list of directed link (each link represent two adjacent components)
        :param id_comp_to_comp:
        :param node_to_id_comp:
        :param segmented_node_to_id_comp:
        :param adj_list:
        """
        self.id = id
        self.times = times
        # Do not touch !
        if c_nodes is None:
            self.c_nodes = []
        else:
            self.c_nodes = c_nodes
        if c_links is None:
            self.c_links = []
        else:
            self.c_links = c_links
        self.node_to_id_comp = node_to_id_comp
        self.segmented_node_to_id_comp = segmented_node_to_id_comp
        self.id_comp_to_comp = id_comp_to_comp
        self.adj_list = adj_list

    def __repr__(self):
        rep = "\nId wcc :" + str(self.id)
        rep += "\nNodes :" + repr(self.c_nodes)
        rep += "\nLinks :" + str(self.c_links)
        return rep

    def describe(self):
        """
        Print a Short description of the current ``Dag`` object.

        :return: Nothing
        """
        print("DAG Id:", self.id)
        print("Nb of Nodes (comp) : ", len(self.c_nodes))
        print("Nb of stream nodes : ", len(set([n for c in self.c_nodes for n in c.nodes])))
        print("Nb of links (betw. comp): ", len(self.c_links))
        stream_links = set()
        degree_in = {n: 0 for n in self.id_comp_to_comp}
        degree_out = {n: 0 for n in self.id_comp_to_comp}
        for c in self.c_nodes:
            if c.links:
                for l in c.links:
                    stream_links.add(tuple(l))
        print("Nb of stream links : ", len(stream_links))

        for l in self.c_links:
            degree_out[l[0]] += 1
            degree_in[l[1]] += 1

        print("Mean In degree :", sum([i for i in degree_in.values()]) / len(degree_in),
              " mean out degree :", sum([i for i in degree_out.values()]) / len(degree_out))

    def set_index_node_to_id_comp(self):
        node_to_id_comp = defaultdict(list)
        # order comp by their starting time
        for c in sorted(self.c_nodes, key=lambda x: x.times[0]):
            for n in c.nodes:
                node_to_id_comp[n].append(c.id)
        self.node_to_id_comp = node_to_id_comp

    def set_index_segmented_node_to_id_comp(self, segmented_nodes):
        if not self.node_to_id_comp:
            self.set_index_node_to_id_comp()
        segmented_node_to_id_comp = defaultdict(list)
        for n in segmented_nodes:
            for id_scc in self.node_to_id_comp[n[2]]:
                c = self.id_comp_to_comp[id_scc]
                if n[0] <= c.times[0] <= n[1]:
                    segmented_node_to_id_comp[n].append(id_scc)
        self.segmented_node_to_id_comp = segmented_node_to_id_comp

    def set_id_comp_to_comp(self):
        id_comp_to_comp = {}
        if type(self.c_nodes[0]) == list:
            cnt_id_scc = 0
            for c in self.c_nodes:
                id_comp_to_comp[cnt_id_scc] = c
                cnt_id_scc += 1
        else:
            for c in self.c_nodes:
                id_comp_to_comp[c.id] = c
        self.id_comp_to_comp = id_comp_to_comp

    def set_id(self, id):
        if id is None:
            self.id = 0
        else:
            self.id = id

    def set_index_id_comp_to_comp(self, index):
        self.id_comp_to_comp = index

    def add_node(self, n):
        self.c_nodes.append(n)

    def add_nodes(self, l):
        self.c_nodes += l

    def add_link(self, l):
        self.c_links.append(l)

    def size(self):
        return len(self.c_nodes)

    def refactor(self):
        return [[(cc.times[0], cc.times[1], n) for n in cc.nodes] for cc in self.c_nodes]

    def compute_links_inplace(self):

        # chrono = time.time()
        # IF there's is some links, discard them.
        if self.c_links:
            self.c_links = []
        dict_begin_time_to_cnodes = defaultdict(list)
        dict_end_time_to_cnodes = defaultdict(list)

        if type(self.c_nodes[0]) == list:
            for i, c in self.id_comp_to_comp.items():
                dict_begin_time_to_cnodes[c[0][0]].append(i)
                dict_end_time_to_cnodes[c[0][1]].append(i)
            # For each end_time find a begin_time (if it exists) then compute the intersection of nodes
            # if not empty add to dag links
            set_links = set()
            a_l = defaultdict(set)
            for te in dict_end_time_to_cnodes:
                if te in dict_begin_time_to_cnodes:
                    # Double for loop in order to match pairwise
                    for cn_parent_id in dict_end_time_to_cnodes[te]:
                        for cn_child_id in dict_begin_time_to_cnodes[te]:
                            if cn_child_id != cn_parent_id:
                                cn_child_nodes = set([c[2] for c in self.id_comp_to_comp[cn_child_id]])
                                cn_parent_nodes = set([c[2] for c in self.id_comp_to_comp[cn_parent_id]])
                                if not cn_child_nodes.isdisjoint(cn_parent_nodes):
                                    set_links.add((cn_parent_id, cn_child_id))
                                    a_l[cn_parent_id].add(cn_child_id)
        else:
            for cn in self.c_nodes:
                dict_begin_time_to_cnodes[cn.times[0]].append(cn.id)
                dict_end_time_to_cnodes[cn.times[1]].append(cn.id)
            # For each end_time find a begin_time (if it exists) then compute the intersection of nodes
            # if not empty add to dag links
            set_links = set()
            a_l = defaultdict(set)
            for te in dict_end_time_to_cnodes:
                if te in dict_begin_time_to_cnodes:
                    # Double for loop in order to match pairwise
                    for cn_parent_id in dict_end_time_to_cnodes[te]:
                        for cn_child_id in dict_begin_time_to_cnodes[te]:
                            if cn_child_id != cn_parent_id:
                                cn_child = self.id_comp_to_comp[cn_child_id]
                                cn_parent = self.id_comp_to_comp[cn_parent_id]
                                if not cn_child.nodes.isdisjoint(cn_parent.nodes):
                                    set_links.add((cn_parent_id, cn_child_id))
                                    a_l[cn_parent_id].add(cn_child_id)
        #
        # print("[DAG] compute link inplace :", time.time() - chrono)
        # chrono = time.time()
        # Remove useless links:
        to_remove = set()
        for l in set_links:
            i, j = l
            for k in a_l[i]:
                if k in a_l:
                    if j in a_l[k]:
                        to_remove.add((i, j))
                        break
        set_links -= to_remove
        # print("[DAG] remove useless links :", time.time() - chrono)
        self.c_links = list(set_links)

    def adjacency_list(self):
        if not self.adj_list:
            a_l = defaultdict(list)
            for l in self.c_links:
                # ONLY IF the destination is accesible from begin time
                # c = self.index_id_scc_to_nodes[l[1]]
                # if c.times[0] >= threshold:
                a_l[l[0]].append(l[1])
            self.adj_list = a_l
        return self.adj_list

    #############################
    #       Plot Functions      #
    #############################

    def plot(self, colors=None, title=None,
             node_to_label=None, fontsize=26, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
        rectangles = []
        labels = []
        unique_nodes = set()
        min_t = math.inf
        max_t = -math.inf
        icol = 0
        if not colors:
            colors = get_cmap(len(self.c_nodes))
        for cn in self.c_nodes:
            if len(cn.nodes) != 1:
                if cn.times[1] == cn.times[0]:
                    rectangles.append(mpatch.Rectangle((cn.times[0], min(cn.nodes)),
                                                       width=0.2 / (self.times[1] - self.times[0]),
                                                       height=max(cn.nodes) - min(cn.nodes) + 0.3,
                                                       edgecolor='k',
                                                       facecolor=colors(icol),
                                                       alpha=0.4
                                                       ))

                else:
                    rectangles.append(mpatch.Rectangle((cn.times[0], min(cn.nodes)),
                                                       width=cn.times[1] - cn.times[0],
                                                       height=max(cn.nodes) - min(cn.nodes) + 0.3,
                                                       edgecolor='k',
                                                       facecolor=colors(icol),
                                                       alpha=0.3
                                                       ))
            else:
                if cn.times[1] == cn.times[0]:
                    rectangles.append(mpatch.Rectangle((cn.times[0], min(cn.nodes)),
                                                       width=0.2 / (self.times[1] - self.times[0]),
                                                       height=0.3,
                                                       edgecolor='k',
                                                       facecolor=colors(icol),
                                                       alpha=0.5
                                                       ))
                else:
                    rectangles.append(mpatch.Rectangle((cn.times[0], min(cn.nodes)),
                                                       width=cn.times[1] - cn.times[0],
                                                       height=0.3,
                                                       edgecolor='k',
                                                       facecolor=colors(icol),
                                                       alpha=0.5
                                                       ))
            icol += 1
            labels.append("\n\n\n".join(reversed([str(n) for n in cn.nodes])))
            unique_nodes |= cn.nodes
            min_t = min(min_t, cn.times[0])
            max_t = max(max_t, cn.times[1])
        icol = 0
        for r, l in zip(rectangles, labels):
            eps = random.choice([1, -1]) * 0.1
            ax.add_artist(r)
            rx, ry = r.get_xy()
            cx = rx + r.get_width() / 2.0 + eps
            cy = ry + r.get_height() / 2.0
            ax.annotate(str(icol), (cx, cy), color='k',
                        fontsize=13, ha='center', va='center')
            icol += 1
        icol = 0
        rectangles = []
        for cn in self.c_nodes:
            if len(cn.nodes) != 1:
                for n in cn.nodes:
                    if cn.times[1] == cn.times[0]:
                        rectangles.append(mpatch.Rectangle((cn.times[0], n),
                                                           width=0.2 / (self.times[1] - self.times[0]),
                                                           height=0.3,
                                                           edgecolor='k',
                                                           facecolor=colors(icol),
                                                           alpha=0.5
                                                           ))

                    else:
                        rectangles.append(mpatch.Rectangle((cn.times[0], n),
                                                           width=cn.times[1] - cn.times[0],
                                                           height=0.3,
                                                           edgecolor='k',
                                                           facecolor=colors(icol),
                                                           alpha=0.5
                                                           ))
            icol += 1
        for r in rectangles:
            ax.add_artist(r)
        ax.set_xlim((min_t - 0.3, max_t + 0.3))
        ax.set_ylim((-0.3, max(unique_nodes) + 0.3))

        # Set xticks
        xticks = numpy.linspace(int(self.times[0]), int(self.times[1]), 11)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks,
                           fontname='Garamond', fontsize=fontsize, color='#476b6b')

        # Set yticks
        yticks = numpy.linspace(min(unique_nodes), int(max(unique_nodes)), min(len(unique_nodes), 25), dtype=int)
        ax.set_yticks(yticks)

        if node_to_label is not None:
            ax.set_yticklabels([node_to_label[i] for i in yticks], fontname='Garamond',
                               fontsize=fontsize,
                               color='#666699')
        else:
            ax.set_yticklabels(yticks, fontname='Garamond', fontsize=fontsize,
                               color='#666699')

        # Set axes label
        ax.set_ylabel("Nodes", fontname='Garamond', fontsize=fontsize, color='#666699')
        ax.set_xlabel("t", fontname='Garamond', fontsize=fontsize, color='#476b6b')
        if title:
            ax.set_title(title, fontname='Garamond', fontsize=fontsize)
        for place, spine in plt.gca().spines.items():
            if place != 'bottom':
                spine.set_visible(False)
            else:
                spine.set_bounds(self.times[0], self.times[1])
                spine.set_color('#476b6b')

        ax.tick_params(top=False, bottom=True, right=False, left=True, labelleft=True, labelbottom=True)

        return ax

    def plot_as_nx(self, label=True):
        fig = plt.figure()
        g_adjacency_list = {i: [] for i, _ in enumerate(self.c_nodes)}
        print("adj_list : ", g_adjacency_list)
        if self.c_links:
            for l in self.c_links:
                n1 = l[0]
                n2 = l[1]
                g_adjacency_list[n1].append(n2)
            G_glob = nx.from_dict_of_lists(g_adjacency_list)
            pos = nx.kamada_kawai_layout(G_glob)
            if label:
                nx.draw_networkx_labels(G_glob, pos, font_size=15)
            nx.draw_networkx_nodes(G_glob, pos, node_size=10,
                                   node_color="#339966", alpha=0.5)
            nx.draw_networkx_edges(G_glob, pos, edge_color='#2d5986',
                                   alpha=0.5, width=1, arrows=True)
        return fig

    def plot_custom(self, label=True, arrow=True, fontsize=30, path=None):
        fig = plt.figure()
        ax = plt.gca()
        pos = defaultdict(list)
        min_x = math.inf
        max_x = -math.inf
        max_y = 0
        min_y = 0
        for i, n in enumerate(self.c_nodes):
            t0, t1 = n.times
            y = numpy.mean(list(n.nodes))
            x = (t0 + t1) / 2
            min_x = min(min_x, t0)
            max_x = max(max_x, t1)
            max_y = max(max_y, y)
            pos[i] = [x, y]

        g_adjacency_list = self.adjacency_list()

        # Plot nodes and links : plot nodes in increasing order
        # adjust 'y' depending on the number of neighbors
        E, E_path = [], []
        if path is not None:
            # If we choose to plot a path on the DAG, we need a second edge_collections
            for n, _ in sorted(pos.items(), key=lambda j: j[1][0]):
                for u in g_adjacency_list[n]:
                    if (n, u) in path or (u, n) in path:
                        E_path.append((pos[n], pos[u]))
                    else:
                        E.append((pos[n], pos[u]))

        else:
            for n, _ in sorted(pos.items(), key=lambda j: j[1][0]):
                for u in g_adjacency_list[n]:
                    E.append((pos[n], pos[u]))
        xy = numpy.asarray([pos[v] for v in pos])
        ax.scatter(xy[:, 0], xy[:, 1],
                   c="#339966",
                   s=180,
                   marker='o',
                   alpha=1
                   )
        if label:
            for n in pos:
                ax.annotate(n, pos[n], fontsize=28)

        edge_collections = mcol.LineCollection(E, colors=['#2d5986'], linewidths=3.5, alpha=0.5)
        if arrow is True:
            for e in E:
                add_arrow(e, ax, direction='right', color='#2d5986')
        ax.add_collection(edge_collections)

        if path:
            edge_collections_path = mcol.LineCollection(E_path, colors=['#8f246b'], linewidths=5, alpha=1)
            if arrow is True:
                for e in E_path:
                    add_arrow(e, ax, direction='right', color='#8f246b')
            ax.add_collection(edge_collections_path)

        ax.set_ylim((min_y - 3, max_y + 3))
        ax.set_xlim((min_x, max_x))

        ax.set_xticks(numpy.linspace(int(min_x), int(max_x), 11))
        ax.set_xlabel("t", fontname='Garamond', fontsize=fontsize, color='#476b6b')
        for place, spine in plt.gca().spines.items():
            if place != 'bottom':
                spine.set_visible(False)
            else:
                spine.set_bounds(self.times[0], self.times[1])
                spine.set_color('#476b6b')
        ax.tick_params(right=False, left=False, labelleft=False,
                       labelbottom=True,
                       labelsize=fontsize,
                       colors='#476b6b')
        return fig
