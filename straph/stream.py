# Copyright (C) 2017-2020 Léo Rannou - Sorbonne Université/LIP6 - Thalès
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

# OVERALL Todos:
# TODO : Harmoniser les noms des fonctions, augmented ou segmented dès que temporal node au lieu de node
#  TODO : Mettre toutes les fonctions concernant les WCC dans la classe correspondante

import copy, csv
import itertools
# import logging
import os
import pathlib
import pickle
import random
import time
import math
import json

import dateutil.parser as du

from collections import Counter, defaultdict

from sortedcollections import SortedList, SortedSet

from warnings import warn

from string import ascii_uppercase

import matplotlib.collections as mcol
import matplotlib.patches as mpatch
import matplotlib.dates as md
import matplotlib.pyplot as plt
import datetime as dt
import msgpack
import networkx as nx
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from matplotlib import animation
from matplotlib import cm
# for p in sys.path:
#     print(p)

from straph.dag import condensation_dag as cdag
from straph.utils import get_cmap, profile_shit
from straph import components as comp
from straph.paths import paths as ap
from straph import generators as gsg
from straph import parser
from straph import EulerTourForest as etf

plt.rcParams["figure.figsize"] = (15, 10)


def postprocess_clusters(dict_clusters):
    '''
    Merge adjacent clusters withe the same value
    :param dict_clusters:  Dictionnary : prop -> clusters
    :return:
    '''
    chrono = time.time()
    for k, cluster in dict_clusters.items():
        # 1. Associate elements in v to a node : node -> Times clusters
        #  2. For each node: delete duplicate times
        # 3. Update dictionary

        d = defaultdict(Counter)
        # d2 = defaultdict(set)
        for t0, t1, u in cluster:
            # if t0 != t1:
            d[u][t0] += 1
            d[u][t1] += 1
            # else:   #  Process instantaneous comps separately
            #     d2[u].add(t0)

        d3 = defaultdict(list)
        #  Remove duplicate elements
        for n in d:
            for t in d[n]:
                if d[n][t] == 1 or d[n][t] == 3:
                    d3[n].append(t)

        dict_clusters[k] = []
        for n, l in d3.items():
            v = sorted(l)
            # dict_clusters[k] += [(t, t, n) for t in d2[n]]
            dict_clusters[k] += [(t0, t1, n) for t0, t1 in zip(v[::2], v[1::2])]
    print("Postprocess clusters :", time.time() - chrono)
    return dict_clusters


def clusters_to_signals(S, dict_clusters, datetime=True):
    '''

    :param S:
    :param dict_clusters:
    :return:
    '''
    # times_index = pd.Index(sorted([v[0] for l in dict_clusters.values() for v in l]
    #                               + [v[1] for l in dict_clusters.values() for v in l]))

    nodes_to_time_to_value = defaultdict(dict)
    #  as intervals can overlap we need to define an eps for the signal !!
    eps = 10 ** -6

    for k, v in dict_clusters.items():
        for clust in v:
            t0, t1, n = clust
            if t0 == t1:
                nodes_to_time_to_value[n][t0] = k
            else:
                nodes_to_time_to_value[n][t0 + eps] = k
                nodes_to_time_to_value[n][t1 - eps] = k

    #  If node aren't present we input -1.
    for n, np in zip(S.nodes, S.node_presence):
        np = [S.times[0]] + np + [S.times[1]]
        for t0, t1 in zip(np[::2], np[1::2]):
            if t0 != t1:
                nodes_to_time_to_value[n][t0 + eps] = -1
                nodes_to_time_to_value[n][t1 - eps] = -1
    node_to_series = {}

    for n in nodes_to_time_to_value:
        dates = [dt.datetime.fromtimestamp(ts) for ts in nodes_to_time_to_value[n]]
        Index = pd.to_datetime(dates)
        node_to_series[n] = pd.Series(list(nodes_to_time_to_value[n].values()), index=Index)
        node_to_series[n].sort_index(inplace=True)
    return node_to_series


##############################################################

def discretize_time_interval(t, bins):
    """
    Discretize a continous interval according to the bins parameter

    :param t: The interval to discretize
    :param bins: The bins used to discretize (gap between each value)
    :return: The discretized interval
    """
    if t[1] - t[0] <= 0:
        raise ValueError("This time interval is not defined")
    # to deal with float precision (bins/10)
    return np.arange(t[0], t[1] + bins / 10, bins)


def get_bounds_for_evenly_spaced_intervals(t, bins):
    """
    Get bounds for one or multiple evenly spaced discrete intervals

    :param t: The interval to get bounds for
    :param bins: The bins used to discretize the intervals
    :return: Bounds of each intervals
    """
    intervals = []
    tmp = [t[0]]
    for i in range(1, len(t)):
        # To deal with float precision (bins/10)
        if (t[i] - t[i - 1]) > (bins + bins / 10):
            tmp.append(t[i - 1])
            intervals.append(tmp)
            tmp = [t[i]]
    tmp.append(t[-1])
    intervals.append(tmp)
    return intervals


def read_stream_graph(path_links, path_nodes=False, node_label=True,
                      path_weights=None, path_trips=None
                      ):
    """
    tb : time of arrival (b: begin)
    te : time of departure (e: end)
    Input format :
        node file:
            id_node1 tb_0 te_0 tb_1 te_1 ... tb_n1 te_n1
            id_node2 tb_0 te_0 ...
            ...
        link file:
            id_node1 id_node2 tb_0 te_0 tb_1 te_1 ... tb_l1 te_l1
            id_node3 id_node4 tb_0 te_0 tb_1 te_1 ... tb_l2 te_l2
            ...

    :param path_links: path to store nodes and their time of presence
    :param path_nodes: path to store links and their time of presence
    :return:
    """
    nodes = []
    node_presence = []
    links = []
    link_presence = []

    if node_label:
        id_to_label = {}
        label_to_id = defaultdict(lambda: len(label_to_id))

    if path_nodes:
        with open(path_nodes, 'r') as file_input:
            for line in file_input:
                line = line.strip().split(" ")
                if len(line) > 1:
                    if node_label:
                        n_label = str(line[0])
                        n = label_to_id[n_label]
                        id_to_label[n] = n_label
                    else:
                        n = int(line[0])
                    nodes.append(n)
                    np = [float(t) for t in line[1:]]
                    node_presence.append(np)
    with open(path_links, 'r') as file_input:
        for line in file_input:
            line = line.strip().split(" ")
            if len(line) > 2:
                if node_label:
                    u_label = str(line[0])
                    v_label = str(line[1])
                    # assert u_label in label_to_id
                    # assert v_label in label_to_id
                    u = label_to_id[u_label]
                    v = label_to_id[v_label]
                else:
                    u = int(line[0])
                    v = int(line[1])
                links.append((u, v))
                if not path_nodes:
                    if u not in nodes:
                        nodes.append(u)
                    if v not in nodes:
                        nodes.append(v)
                lp = [float(t) for t in line[2:]]
                link_presence.append(lp)

    weights = []
    if path_weights:
        with open(path_weights, 'r') as file_input:
            for line in file_input:
                line = line.strip().split(" ")
                if len(line) > 2:
                    w = [float(t) for t in line[2:]]
                    weights.append(w)
    trips = []
    if path_trips:
        with open(path_trips, 'r') as file_input:
            for line in file_input:
                line = line.strip().split(" ")
                if len(line) > 2:
                    tp = [float(t) for t in line[2:]]
                    trips.append(tp)

    S = stream_graph(times=[min([t for n in node_presence for t in n]),
                            max([t for n in node_presence for t in n])],
                     nodes=nodes,
                     node_presence=node_presence,
                     links=links,
                     link_presence=link_presence,
                     node_to_label=(id_to_label if id_to_label else None),
                     weights=(weights if weights else None),
                     trips=(trips if trips else None))
    return S


def sum_presence(np):
    return sum([t1 - t0 for t0, t1 in zip(np[::2], np[1::2])])


def algo_kcores_batagelj(a_l, degrees):
    """
    Compute k_cores of a static graph from its adjacency list and nodes degrees

    :return:
    References
    ----------
    [1] An O(m) Algorithm for Cores Decomposition of Networks
    Vladimir Batagelj and Matjaz Zaversnik, 2003.
    http://arxiv.org/abs/cs.DS/0310049
    """
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
    return cores


def read_stream_graph_from_sgf(input, node_label=False):
    unpacker = msgpack.Unpacker(open(input, 'rb'), use_list=False)
    node_presence = {}
    link_presence = {}
    min_t = math.inf
    max_t = -math.inf
    id_to_label = {}
    if node_label:
        label_to_id = defaultdict(lambda: len(label_to_id))

    for i in unpacker:
        if i[0] == 1:
            t0, t1 = i[1], i[2]
            min_t = min(t0, min_t)
            max_t = max(t1, max_t)
            if node_label:
                u_label, v_label = str(i[3]), str(i[4])
                u, v = label_to_id[u_label], label_to_id[v_label]
                id_to_label[u], id_to_label[v] = u_label, v_label
            else:
                u, v = int(i[3]), int(i[4])

            if u in node_presence:
                if node_presence[u][-1] >= t0 >= node_presence[u][-2]:
                    node_presence[u][-1] = max(t1, node_presence[u][-1])
                elif t0 > node_presence[u][-1]:
                    node_presence[u] += [t0, t1]
            else:
                node_presence[u] = [t0, t1]
            if v in node_presence:
                if node_presence[v][-1] >= t0 >= node_presence[v][-2]:
                    node_presence[v][-1] = max(t1, node_presence[v][-1])
                elif t0 > node_presence[v][-1]:
                    node_presence[v] += [t0, t1]
            else:
                node_presence[v] = [t0, t1]

            if (u, v) in link_presence:
                if link_presence[(u, v)][-1] >= t0 >= link_presence[(u, v)][-2]:
                    link_presence[(u, v)][-1] = max(t1, link_presence[(u, v)][-1])
                elif t0 > link_presence[(u, v)][-1]:
                    link_presence[(u, v)] += [t0, t1]
            else:
                link_presence[(u, v)] = [t0, t1]

    S = stream_graph(times=[min_t, max_t],
                     nodes=[n for n in node_presence],
                     node_presence=[k for k in node_presence.values()],
                     links=[l for l in link_presence],
                     link_presence=[k for k in link_presence.values()],
                     node_to_label=(id_to_label if id_to_label else None)
                     )

    print("\n AFTER READING FROM SGF:")
    print("max nodes :", max(S.nodes))
    print("max node in links :", max([max(i) for i in S.links]))
    return S


def read_stream_graph_from_json(path_nodes, path_links):
    nodes = []
    node_presence = []
    links = []
    link_presence = []
    nodes_to_id = {}
    id_to_node = {}
    nodes_to_label = {}
    with open(path_nodes, 'r') as file_input:
        nodes_json = json.load(file_input)
        for js in nodes_json["nodes"]:
            n_id = int(js["id"])
            if 'label' in js:
                n_label = str(js["label"])
            else:
                n_label = n_id
            n = len(nodes)
            nodes_to_id[n] = n_id
            id_to_node[n_id] = n
            nodes_to_label[n] = n_label
            nodes.append(n)
            np = []
            for i in js["intervals"]:
                np += [i["t0"], i["t1"]]
            np = sorted(np)
            # Remove Duplicatas
            id_to_remove = set()
            for i in range(len(np) - 1):
                if np[i] == np[i + 1]:
                    id_to_remove.add(i)
                    id_to_remove.add(i + 1)
            np = [np[i] for i in range(len(np)) if i not in id_to_remove]
            node_presence.append(np)
        times = nodes_json["timeExtent"]

    with open(path_links, 'r') as file_input:
        links_json = json.load(file_input)
        for js in links_json["links"]:
            u_id, v_id = int(js['node1']), int(js['node2'])
            links.append((id_to_node[u_id], id_to_node[v_id]))
            lp = []
            for i in js["intervals"]:
                lp += [i["t0"], i["t1"]]
            link_presence.append(lp)
    # print("nodes :",nodes)
    # print("node presence :",node_presence)
    # print("links :",links)
    # print("link presence :",link_presence)

    S = stream_graph(times=times,
                     nodes=nodes,
                     node_presence=node_presence,
                     links=links,
                     link_presence=link_presence,
                     node_to_label=nodes_to_label,
                     node_to_id=nodes_to_id)
    return S


class stream_graph:
    def __init__(self,
                 id=None,
                 times=None,
                 nodes=None,
                 node_to_label=None,
                 node_to_id=None,
                 node_presence=None,
                 links=None,
                 link_presence=None,
                 bins=None,
                 card_E=None,
                 card_W=None,
                 weights=None,
                 trips=None):
        """
        A basic constructor for a stream graph object
        :param id: A parameter to identify a stream graph.
        :param times: Continous interval of time during the stream graph exists
        :param nodes: A list of nodes present in  the stream graph
        :param node_presence : List of lists in the same order as the nodes. Each list contains
        succescively the time of apparition and time of disparition of the node.
        :param links : A list of links present in the stream graph
        :param link_presence : same as node_presence
        :param bins : If we are in the discrete case or if later we want to discretize the stream
                     graph
        """
        # Continuous intervals or discrete times ?
        self.id = id
        self.times = times
        self.nodes = nodes
        self.node_to_label = node_to_label
        self.node_to_id = node_to_id
        self.node_presence = node_presence
        self.links = links
        self.link_presence = link_presence
        self.bins = bins
        self.card_E = card_E
        self.card_W = card_W
        self.weights = weights
        self.trips = trips

    # for continuous stream graphs
    # Easily optimisable, suffit de comparer les éléments de liste pairwise puis de double check if TRUE
    def check_integrity(self):
        """
        Check node presence and link presence for overlapping time windows
        Check that a link presence include both node presence
        :return: True if the structure is coherent
                 An error with the problematic link/node
        """

        for n, np in zip(self.links, self.link_presence):
            assert list(sorted(np)) == np
            for t0, t1 in zip(np[::2], np[1::2]):
                for t2, t3 in zip(np[::2], np[1::2]):
                    if (t0, t1) != (t2, t3):
                        if t0 <= t3 and t2 <= t1:
                            raise ValueError("Integrity compromised on node : "
                                             + str(n) +
                                             ". Check node presence (probably overlapping intervals) !\n")

        for l, lp in zip(self.links, self.link_presence):
            assert list(sorted(lp)) == lp
            id1 = self.nodes.index(l[0])
            id2 = self.nodes.index(l[1])
            for lt0, lt1 in zip(lp[::2], lp[1::2]):
                check_1 = False
                check_2 = False
                for nt0, nt1 in zip(self.node_presence[id1][::2],
                                    self.node_presence[id1][1::2]):
                    if nt0 <= lt0 and nt1 >= lt1:
                        check_1 = True
                        continue
                for nt0, nt1 in zip(self.node_presence[id2][::2],
                                    self.node_presence[id2][1::2]):
                    if nt0 <= lt0 and nt1 >= lt1:
                        check_2 = True
                        continue
                if not check_1 or not check_2:
                    raise ValueError("Integrity compromised on link : "
                                     + str(l) + " at time :" + str((lt0, lt1)) +
                                     ". Check node presence !\n" + "Node " + str(l[0]) +
                                     " presence:" + str(self.node_presence[l[0]]) + "\nNode " + str(l[1]) +
                                     " presence:" + str(self.node_presence[l[1]]))
        print("Integrity check ok !")
        return True

    def add_node(self, node, node_presence):
        """
        Add a node to the stream graph, if not already present :
        this new node will be represented by a an integer, correponding to the length of self.nodes
        :param node: Label of the node to add
        :param node_presence: Presence times of the added node [b,e,b',e',...]
        :return: The id corresponding to the added node
        """

        if self.node_to_label is not None:
            label_to_node = {v: k for k, v in self.node_to_label.items()}
            if node in label_to_node:
                new_node = label_to_node[node]
                assert self.node_presence[new_node][-1] < node_presence[0]
                self.node_presence[new_node] += node_presence
            else:
                new_node = len(self.nodes)
                self.node_to_label[new_node] = node
                self.nodes.append(new_node)
                self.node_presence.append(node_presence)
        else:
            if node < len(self.nodes):
                new_node = node
                assert self.node_presence[node][-1] < node_presence[0]
                self.node_presence[node] += node_presence
            else:
                new_node = len(self.nodes)
                self.nodes.append(new_node)
                self.node_presence.append(node_presence)
        return new_node

    def add_link(self, link, link_presence):
        """
        Add a new link to the stream graph. If the extremities are not already present we add them to the stream graph.
        Their presence times will be the one the link.
        :param link: (A,B) where 'A' and 'B' are the nodes' labels
        :param link_presence: Presence times of the added link [b,e,b',e',...]
        :return:
        """
        u, v = link
        if self.node_to_label is not None:
            label_to_node = {v: k for k, v in self.node_to_label.items()}
            if u in label_to_node:
                new_u = label_to_node[u]
            else:
                new_u = self.add_node(u, link_presence)
            if v in label_to_node:
                new_v = label_to_node[v]
            else:
                new_v = self.add_node(v, link_presence)
        else:
            if u in self.nodes:
                new_u = u
            else:
                new_u = self.add_node(u, link_presence)
            if v in self.nodes:
                new_v = v
            else:
                new_v = self.add_node(v, link_presence)

        if (new_u, new_v) in self.links:
            assert self.link_presence[self.links.index((new_u, new_v))][-1] < link_presence[0]
            self.link_presence[self.links.index((new_u, new_v))] += link_presence
        elif (new_v, new_u) in self.links:
            assert self.link_presence[self.links.index((new_v, new_u))][-1] < link_presence[0]
            self.link_presence[self.links.index((new_v, new_u))] += link_presence
        else:
            new_link = (new_u, new_v)
            self.links.append((new_link))
            self.link_presence.append(link_presence)

    #####################################
    #       Events Representation       #
    #####################################

    def event_times(self, all=False):
        """
        Return distinct event times in the Stream Graph
        :return:
        """
        if all:
            event_times = sorted([t for np in self.node_presence for t in np] + \
                                 [t for lp in self.link_presence for t in lp])
            return event_times
        else:
            event_t_nodes = set([t for np in self.node_presence for t in np])
            event_t_links = set([t for lp in self.link_presence for t in lp])
            return event_t_links | event_t_nodes

    def number_of_event_times(self):
        """
        Return the number of distinct event times in the Stream Graph
        :return:
        """
        return len(self.event_times())

    def ordered_events(self, weights_or_trips=False):
        """
        :return:
        """
        links = []
        if weights_or_trips:
            if self.weights and self.trips:
                for l, lp, we, tr in zip(self.links, self.link_presence, self.weights, self.trips):
                    for t0, t1, w, d in zip(lp[::2], lp[1::2], we, tr):
                        u, v = l
                        links.append((1, t0, t1, u, v, w, d))  # code each link, 1 for a beginning, -1 for an ending
                        links.append((-1, t1, u, v, w, d))

            elif self.weights:
                for l, lp, we in zip(self.links, self.link_presence, self.weights):
                    for t0, t1, w in zip(lp[::2], lp[1::2], we):
                        u, v = l
                        links.append((1, t0, t1, u, v, w, 0))  # code each link, 1 for a beginning, -1 for an ending
                        links.append((-1, t1, u, v, w, 0))

            elif self.trips:
                for l, lp, tr in zip(self.links, self.link_presence, self.trips):
                    for t0, t1, d in zip(lp[::2], lp[1::2], tr):
                        u, v = l
                        links.append((1, t0, t1, u, v, 1, d))  # code each link, 1 for a beginning, -1 for an ending
                        links.append((-1, t1, u, v, 1, d))
            else:
                for l, lp in zip(self.links, self.link_presence):
                    for t0, t1 in zip(lp[::2], lp[1::2]):
                        u, v = l
                        links.append((1, t0, t1, u, v, 1, 0))  # code each link, 1 for a beginning, -1 for an ending
                        links.append((-1, t1, u, v, 1, 0))
        else:
            for l, lp in zip(self.links, self.link_presence):
                for t0, t1 in zip(lp[::2], lp[1::2]):
                    u, v = l
                    links.append((1, t0, t1, u, v))  # code each link, 1 for a beginning, -1 for an ending
                    links.append((-1, t1, u, v))
        nodes = []
        for n, np in zip(self.nodes, self.node_presence):
            for t0, t1 in zip(np[::2], np[1::2]):
                nodes.append((2, t0, t1, n))  #  code a node arrival with a 2
                nodes.append((-2, t1, n))  #  code a node departure with a -2
        events = sorted(links + nodes, key=lambda x: (x[1], -x[0]))
        return events

    def ordered_batch_events(self):
        """
        Return an ordered of batch events (event happening at the same time).
        In each batch there can be (in the same order): node arrival, link arrival, link departure, node departure.
        :return:
        """

        chrn = time.time()
        batchs = []
        E = self.ordered_events()
        t_old = E[0][1]
        c_old = E[0][0]
        current_batch = []
        for e in E:
            c, t = e[0], e[1]
            if t == t_old and c == c_old:
                current_batch.append(e)
            else:
                batchs.append(current_batch)
                current_batch = [e]
                t_old = t
                c_old = c

        if current_batch:
            batchs.append(current_batch)
        # print("[TIME] Ordered Batch Events :", time.time() - chrn)
        return batchs

    #######################################
    #       Arrivals Events               #
    #######################################

    def ordered_arrivals(self, free_memory=False):
        """
        :return:
        """

        nodes = []
        for n, np in zip(self.nodes, self.node_presence):
            for t0, t1 in zip(np[::2], np[1::2]):
                nodes.append((2, t0, t1, n))  #  code a node arrival with a 2
        if free_memory:
            self.nodes = []
            self.node_presence = []

        links = []
        for l, lp in zip(self.links, self.link_presence):
            for t0, t1 in zip(lp[::2], lp[1::2]):
                u, v = l
                links.append((1, t0, t1, u, v))  # code each link, 1 for a beginning, -1 for an ending

        if free_memory:
            self.links = []
            self.link_presence = []

        events = sorted(links + nodes, key=lambda x: (x[1], -x[0]))
        return events

    #####################################
    #       Links Representation        #
    #####################################

    def augmented_ordered_links(self, weights_or_trips=False):
        """
        :return:
        """
        links = []
        if weights_or_trips:
            if self.weights and self.trips:
                for l, lp, we, tr in zip(self.links, self.link_presence, self.weights, self.trips):
                    for t0, t1, w, d in zip(lp[::2], lp[1::2], we, tr):
                        nu, nv = self.get_according_node_presence_from_link(l, t0, t1)
                        if nu and nv:
                            links.append(
                                (1, t0, t1, nu, nv, w, d))  # code each link, 1 for a beginning, -1 for an ending
                            links.append((-1, t1, nu, nv, w, d))

            elif self.weights:
                for l, lp, we in zip(self.links, self.link_presence, self.weights):
                    for t0, t1, w in zip(lp[::2], lp[1::2], we):
                        nu, nv = self.get_according_node_presence_from_link(l, t0, t1)
                        if nu and nv:
                            links.append(
                                (1, t0, t1, nu, nv, w, 0))  # code each link, 1 for a beginning, -1 for an ending
                            links.append((-1, t1, nu, nv, w, 0))

            elif self.trips:
                for l, lp, tr in zip(self.links, self.link_presence, self.trips):
                    for t0, t1, d in zip(lp[::2], lp[1::2], tr):
                        nu, nv = self.get_according_node_presence_from_link(l, t0, t1)
                        if nu and nv:
                            links.append(
                                (1, t0, t1, nu, nv, 1, d))  # code each link, 1 for a beginning, -1 for an ending
                            links.append((-1, t1, nu, nv, 1, d))
            else:
                for l, lp in zip(self.links, self.link_presence):
                    for t0, t1 in zip(lp[::2], lp[1::2]):
                        nu, nv = self.get_according_node_presence_from_link(l, t0, t1)
                        if nu and nv:
                            links.append(
                                (1, t0, t1, nu, nv, 1, 0))  # code each link, 1 for a beginning, -1 for an ending
                            links.append((-1, t1, nu, nv, 1, 0))
        else:
            for l, lp in zip(self.links, self.link_presence):
                for t0, t1 in zip(lp[::2], lp[1::2]):
                    nu, nv = self.get_according_node_presence_from_link(l, t0, t1)
                    if nu and nv:
                        links.append((1, t0, t1, nu, nv))  # code each link, 1 for a beginning, -1 for an ending
                        links.append((-1, t1, nu, nv))
        links = sorted(links, key=lambda x: (x[1], -x[0]))
        return links

    def ordered_batch_links(self, free_memory=False):
        chrono_obl = time.time()
        links = []
        batchs = []
        for l, lp in zip(self.links, self.link_presence):
            u, v = l
            for t0, t1 in zip(lp[::2], lp[1::2]):
                links.append((1, t0, t1, u, v))  # code each link, 1 for a beginning, -1 for an ending
                links.append((-1, t1, u, v))
        if free_memory:
            self.nodes = []
            self.node_presence = []
            self.links = []
            self.link_presence = []

        links = sorted(links, key=lambda x: (x[1], -x[0]))
        t_old = links[0][1]
        c_old = links[0][0]
        current_batch = []
        for e in links:
            c, t = e[0], e[1]
            if t == t_old and c == c_old:
                current_batch.append(e)
            else:
                batchs.append(current_batch)
                current_batch = [e]
                t_old = t
                c_old = c
        if current_batch:
            batchs.append(current_batch)
        # print("[TIME] Ordered Batch Links :", time.time() - chrono_obl)
        return batchs

    def ordered_links(self):
        """
        :return:
        """
        links = []
        if self.weights and self.trips:
            for l, lp, we, tr in zip(self.links, self.link_presence, self.weights, self.trips):
                u, v = l
                for t0, t1, w, d in zip(lp[::2], lp[1::2], we, tr):
                    links.append((1, t0, t1, u, v, w, d))  # code each link, 1 for a beginning, -1 for an ending
                    links.append((-1, t1, u, v, w, d))
            links = sorted(links, key=lambda x: (x[1], -x[0]))

        elif self.weights:
            for l, lp, we in zip(self.links, self.link_presence, self.weights):
                u, v = l
                for t0, t1, w in zip(lp[::2], lp[1::2], we):
                    links.append((1, t0, t1, u, v, w))  # code each link, 1 for a beginning, -1 for an ending
                    links.append((-1, t1, u, v, w))
            links = sorted(links, key=lambda x: (x[1], -x[0]))

        elif self.trips:
            for l, lp, tr in zip(self.links, self.link_presence, self.trips):
                u, v = l
                for t0, t1, d in zip(lp[::2], lp[1::2], tr):
                    links.append((1, t0, t1, u, v, d))  # code each link, 1 for a beginning, -1 for an ending
                    links.append((-1, t1, u, v, d))
            links = sorted(links, key=lambda x: (x[1], -x[0]))
        else:
            for l, lp in zip(self.links, self.link_presence):
                u, v = l
                for t0, t1 in zip(lp[::2], lp[1::2]):
                    links.append((1, t0, t1, u, v))  # code each link, 1 for a beginning, -1 for an ending
                    links.append((-1, t1, u, v))
            links = sorted(links, key=lambda x: (x[1], -x[0]))
        return links

    ####################################################################
    #               WRITERS                                            #
    ####################################################################

    def write_to_sg(self, output_file):
        """
        tb : time of arrival (b: begin)
        te : time of departure (e: end)
        Output format :
            node file:
                id_node1 tb_0 te_0 tb_1 te_1 ... tb_n1 te_n1
                id_node2 tb_0 te_0 ...
                ...
            link file:
                id_node1 id_node2 tb_0 te_0 tb_1 te_1 ... tb_l1 te_l1
                id_node3 id_node4 tb_0 te_0 tb_1 te_1 ... tb_l2 te_l2
                ...
        :param path_links: path to store nodes and their time of presence
        :param path_nodes: path to store links and their time of presence
        :return:
        """
        if self.node_to_label:
            with open(output_file + '_nodes.sg', 'w') as file_output:
                for n, np in zip(self.nodes, self.node_presence):
                    file_output.write(str(self.node_to_label[n]) + " ")
                    for t in np:
                        file_output.write(str(t) + " ")
                    file_output.write("\n")
            with open(output_file + '_links.sg', 'w') as file_output:
                for l, lp in zip(self.links, self.link_presence):
                    file_output.write(str(self.node_to_label[l[0]]) + " " + str(self.node_to_label[l[1]]) + " ")
                    for t in lp:
                        file_output.write(str(t) + " ")
                    file_output.write("\n")
        else:
            with open(output_file + '_nodes.sg', 'w') as file_output:
                for n, np in zip(self.nodes, self.node_presence):
                    file_output.write(str(n) + " ")
                    for t in np:
                        file_output.write(str(t) + " ")
                    file_output.write("\n")
            with open(output_file + '_links.sg', 'w') as file_output:
                for l, lp in zip(self.links, self.link_presence):
                    file_output.write(str(l[0]) + " " + str(l[1]) + " ")
                    for t in lp:
                        file_output.write(str(t) + " ")
                    file_output.write("\n")

    def write_to_sgf(self, output_file):
        """
        tb : time of arrival (b: begin)
        te : time of departure (e: end)
        Output format :
            link file:
                1 (id_node1,id_node2) tb_0 te_0
                -1 (id_node3,id_node4) te_0
                ...
        :param output_file: path to store links and their time of presence
        :return:
        """
        packer = msgpack.Packer()
        links = []

        if self.node_to_label:
            for l, lp in zip(self.links, self.link_presence):
                for t0, t1 in zip(lp[::2], lp[1::2]):
                    links.append((1, t0, t1, self.node_to_label[l[0]],
                                  self.node_to_label[l[1]]))  # code each link, 1 for a beginning, -1 for an ending
                    links.append((-1, t1, self.node_to_label[l[0]], self.node_to_label[l[1]]))
        else:
            for l, lp in zip(self.links, self.link_presence):
                for t0, t1 in zip(lp[::2], lp[1::2]):
                    links.append((1, t0, t1, l[0], l[1]))  # code each link, 1 for a beginning, -1 for an ending
                    links.append((-1, t1, l[0], l[1]))
        # print("Links :",links)
        links = sorted(links, key=lambda x: (x[1], -x[0]))
        # print("Sorted Links :",links)
        with open(output_file + "_ordered_links.sgf", 'wb') as file_output:
            for l in links:
                file_output.write(packer.pack(l))

    def write_to_csv(self, output_name):
        """
        Write the stream graph to CSV format(node1;node2;start_time;duration).
        :param output_name:
        :return:
        """
        links = []
        if self.node_to_label:
            for l, lp in zip(self.links, self.link_presence):
                for t0, t1 in zip(lp[::2], lp[1::2]):
                    links.append((t0, t1 - t0, self.node_to_label[l[0]],
                                  self.node_to_label[l[1]],))
        else:
            for l, lp in zip(self.links, self.link_presence):
                for t0, t1 in zip(lp[::2], lp[1::2]):
                    links.append((t0, t1 - t0, l[0], l[1],))
        links = sorted(links, key=lambda x: (x[0], x[1]))
        with open(output_name, 'w', newline='') as file_output:
            stream_writer = csv.writer(file_output, delimiter=';')
            for l in links:
                stream_writer.writerow(l)

    def write_to_lsf(self, output_name):
        links = []
        if self.node_to_label:
            for l, lp in zip(self.links, self.link_presence):
                for t0, t1 in zip(lp[::2], lp[1::2]):
                    links.append((t0, t1, self.node_to_label[l[0]],
                                  self.node_to_label[l[1]]))
        else:
            for l, lp in zip(self.links, self.link_presence):
                for t0, t1 in zip(lp[::2], lp[1::2]):
                    links.append((t0, t1, l[0], l[1]))
        links = sorted(links, key=lambda x: (x[2], x[3]))
        links = [str(l[0]) + " " + str(l[1]) + " " + str(l[2]) + " " + str(l[3])
                 for l in links]
        with open(output_name + ".lsf", 'w', newline='') as otp:
            otp.write("alpha " + str(self.times[0]) + "\n")
            otp.write("omega " + str(self.times[1]) + "\n")
            for l in links:
                otp.write(l + "\n")

    def write_to_json(self, output_name, index_pos=None):
        """
        Write to a JSON format. For the tool stream-graph-visualisation.
        :param output_name: Name prefix of json files. Will be stored under "output_name_node_activity.json"
                            and "output_name_link_presence.json"
        :return:
        """
        nodes_json = self.node_activity_to_json(index_pos)
        links_json = self.link_presence_to_json()
        with open(output_name + "_node_activity.json", 'w') as file_output:
            json.dump(nodes_json, file_output)
        with open(output_name + "_link_presence.json", 'w') as file_output:
            json.dump(links_json, file_output)

    def node_activity_to_json(self, node_to_position=None, sort_type=None):
        """
        For the tool stream-graph-visualisation.
        :return: node activity in JSON
        """

        if self.node_to_id is None:
            self.node_to_id = {n: n for n in self.nodes}

        degrees_partition = self.degrees_partition()
        isolated_nodes = self.get_isolated_nodes()
        degrees_partition[0] = isolated_nodes
        # print("Degrees partition :",degrees_partition)
        if node_to_position is None:
            if self.node_to_id:
                node_to_position = {self.node_to_id[n]: n for n in self.nodes}
            else:
                node_to_position = {n: n for n in self.nodes}
        if sort_type == "arrival":
            node_to_position = self.node_position_by_arrival()
        elif sort_type == "increasing_degree":
            node_to_position = self.node_position_by_increasing_degree(degrees_partition)
        elif sort_type == "peak_degree_arrival":
            node_to_position = self.node_position_by_peak_degree_arrival(degrees_partition)
        # print("SORT TYPE :",sort_type)
        # print(" NODE TO POSITION :",node_to_position)
        max_d = 0
        min_d = 0
        nodes_json = {}
        if self.node_to_label:
            nodes_json["nodes"] = [{"id": (self.node_to_id[n] if self.node_to_id else n),
                                    "label": str(self.node_to_label[n]),
                                    "pos": node_to_position[(self.node_to_id[n] if self.node_to_id else n)],
                                    "intervals": []}
                                   for n in self.nodes]
            for degree, nodes_set in degrees_partition.items():
                for n in nodes_set:
                    t0, t1, u = n
                    if u in self.nodes:  # important pour la visualisation si le noeud est dans les liens mais
                        # pas la liste des noeuds
                        max_d = max(degree, max_d)
                        min_d = min(degree, min_d)
                        nodes_json["nodes"][u]["intervals"].append({"t0": t0, "t1": t1, "degree": degree})
        else:
            nodes_json["nodes"] = [{"id": (self.node_to_id[n] if self.node_to_id else n),
                                    "label": n,
                                    "pos": node_to_position[(self.node_to_id[n] if self.node_to_id else n)],
                                    "intervals": []}
                                   for n in self.nodes]
            for degree, nodes_set in degrees_partition.items():
                for n in nodes_set:
                    t0, t1, u = n
                    nodes_json["nodes"][u]["intervals"].append({"t0": t0, "t1": t1, "degree": degree})

        nodes_json["timeExtent"] = [self.times[0], self.times[1]]
        nodes_json["degreeExtent"] = [min_d, max_d]
        return nodes_json

    def link_presence_to_json(self, link_to_position=None, sort_type=None):
        """
        For the tool stream-graph-visualisation.
        :return: link presence in JSON
        """
        links_json = {"links": [], "timeExtent": [self.times[0], self.times[1]]}
        max_duration_link = 0
        min_duration_link = self.times[1] - self.times[0]

        if link_to_position is None:
            if self.node_to_id:
                link_to_position = {(self.node_to_id[l[0]], self.node_to_id[l[1]]): i for i, l in
                                    enumerate(self.links)}
            else:
                link_to_position = {l: i for i, l in
                                    enumerate(self.links)}
        if sort_type == "duration":
            link_to_position = self.link_position_by_duration()
        elif sort_type == "arrival":
            link_to_position = self.link_position_by_arrival()
        # print("SORT TYPE :",sort_type)
        # print("LINK TO POSITION :",link_to_position)

        for l, lp in zip(self.links, self.link_presence):
            link_presence = []
            duration = 0
            for t0, t1 in zip(lp[::2], lp[1::2]):
                link_presence.append({'t0': t0, 't1': t1})
                duration += t1 - t0
            max_duration_link = max(max_duration_link, duration)
            min_duration_link = min(min_duration_link, duration)
            links_json["links"].append({"node1": (str(self.node_to_id[l[0]]) if self.node_to_id else l[0]),
                                        "node2": (str(self.node_to_id[l[1]]) if self.node_to_id else l[1]),
                                        "intervals": link_presence,
                                        "duration": duration,
                                        "pos": (link_to_position[(self.node_to_id[l[0]], self.node_to_id[l[1]])]
                                                if self.node_to_id else link_to_position[l])})

        links_json["durationExtent"] = [min_duration_link, max_duration_link]
        return links_json

    def describe(self):
        """
        TODO : Add multiple information
        :return:
        """
        print("Nb of Nodes : ", len(self.nodes))
        print("Nb of segmented nodes : ", sum([len(item) / 2 for item in self.node_presence]))
        print("Nb of links : ", len(self.links))
        print("Nb of segmented links : ", sum([len(item) / 2 for item in self.link_presence]))
        print("Nb of event times : ", self.number_of_event_times())

    ####################################################################
    #               Static Graph Properties                            #
    ####################################################################

    def graph_property(self, node_property_function_nx,
                       format="cluster",
                       stable_components=None,
                       postprocess=True,
                       n_jobs=-1,
                       datetime=True):
        if stable_components is None:
            stable_components = self.stable_connected_components(format="object_with_links")
        if format != "cluster" and format != "signal":
            raise KeyError("Read the API documentation #USERISLIMITED")

        def para_prop(sc):
            prop_to_clusters = defaultdict(list)
            a_l = sc.to_adjacency_list()
            if a_l:
                G = nx.from_dict_of_lists(a_l)
                prop = node_property_function_nx(G)
                for n in prop:
                    prop_to_clusters[prop[n]].append((sc.times[0], sc.times[1], n))
            else:  #  Isolated nodes (default value 0)
                for n in sc.nodes:
                    prop_to_clusters[0].append((sc.times[0], sc.times[1], n))
            return prop_to_clusters

        r = Parallel(n_jobs=n_jobs)(delayed(para_prop)(sc) for sc in stable_components)
        prop_to_clusters = defaultdict(list)
        for d in r:
            for k, v in d.items():
                prop_to_clusters[k] += v

        if postprocess == True:
            prop_to_clusters = postprocess_clusters(prop_to_clusters)

        if format == "cluster":
            return prop_to_clusters

        if format == "signal":
            return clusters_to_signals(self, prop_to_clusters, datetime=datetime)

    def graph_communities(self, community_function_nx,
                          stable_components=None, n_jobs=-1):
        # TODO: Not suported yet, must finish and parallelized this shit
        if stable_components is None:
            stable_components = self.stable_connected_components(format="object")
        communities_clusters = []
        for sc in stable_components:
            G = nx.from_dict_of_lists(sc.to_adjacency_list())
            communities = community_function_nx(G)
            for c in communities:
                cluster = [(sc.times[0], sc.times[1], n) for n in c]
            communities_clusters.append(cluster)
        return communities_clusters

    ####################################################################
    #               PLOT                                               #
    ####################################################################

    # For continuous stream graphs
    def plot(self, animated=False, clusters=None, links=True, title=None,
             link_pimping=False, repeat=False, saving_path=None, format='pdf',
             label_nodes_by_letters=False, fontsize=26, xtimestamp=False):
        """
        Display in an elegant way a small stream graph
        We can also specify a path to save the animated of ploted stream graph
        :param animated: If we want an animation with cursor which move according to the time
        :return: A matplotlib plot
        """
        # lnodes = len(self.nodes)
        lnodes = max(self.nodes) + 1
        if clusters:
            random.seed(3)
            c_map_cluster = get_cmap(len(clusters) + 2)
            #  list_colors = [c_map_cluster(c) for c in range(len(clusters))]
            # random.shuffle(list_colors)
        c_map = get_cmap(lnodes)
        l_colors = [c_map(i) for i in range(lnodes)]
        # random.shuffle(l_colors)
        fig = plt.figure()
        if animated:
            ax = plt.axes(xlim=(self.times[0] - 1,
                                self.times[1] + 1), ylim=(-0.3, lnodes))
        # Plot Clusters
        if clusters:
            for c in range(len(clusters)):
                for i in clusters[c]:
                    plt.hlines(i[2], xmin=i[0], linewidth=5,
                               xmax=i[1], color=c_map_cluster(c), alpha=0.7)
        else:
            #  PLOT NODES
            for p, t in zip(self.nodes, self.node_presence):
                coln = l_colors[p]
                # plt.axhline(p, linestyle='--', linewidth=0.7,
                #             color=coln, alpha=0.3)
                plt.hlines([p] * (len(t) // 2), xmin=t[::2], linewidth=3,
                           xmax=t[1::2], colors=coln, alpha=1)
        # Plot Links
        if links:
            for k, lp in zip(self.links, self.link_presence):
                id1 = k[0]
                id2 = k[1]
                coln1 = l_colors[id1]
                coln2 = l_colors[id2]
                idmax = max(id1, id2)
                idmin = min(id1, id2)
                eps = random.choice([1, -1]) * (random.random() / 4)
                if link_pimping:
                    plt.hlines([(idmax + idmin) / 2 + eps] * (len(lp) // 2), xmin=lp[::2], xmax=lp[1::2],
                               linestyles=(3, (3, 3)),
                               color=coln1,
                               linewidth=2.5,
                               alpha=0.6)
                    plt.hlines([(idmax + idmin) / 2 + eps] * (len(lp) // 2), xmin=lp[::2], xmax=lp[1::2],
                               linestyles=(0, (3, 3)),
                               color=coln2,
                               linewidth=2.5,
                               alpha=0.6)
                else:
                    plt.hlines([(idmax + idmin) / 2 + eps] * (len(lp) // 2), xmin=lp[::2], xmax=lp[1::2],
                               colors='k', linewidth=3, alpha=0.5)
                plt.vlines(lp[::2], ymin=idmin, ymax=idmax,
                           linewidth=2, alpha=0.15)
                # ARRIVALS
                plt.plot([lp[::2]], [idmin], color='#004d00', marker='^', alpha=1, markersize=7)
                plt.plot([lp[::2]], [idmax], color='#004d00', marker='v', alpha=1, markersize=7)
                # DEPARTURE
                # plt.plot([lp[1::2]],[idmin,],color='#800000',marker='v',alpha =1,markersize=7)
                # plt.plot([lp[1::2]],[idmax],color='#800000',marker='^',alpha =1,markersize=7)
                # Colored Markers :
                # markerstyle1 = dict(color=l_colors[idmax], marker='>',
                #     markersize=8)
                # markerstyle2 = dict(color=l_colors[idmin],  marker='>',
                #     markersize=8)
                #
                # plt.plot([lp[::2]], [idmin],alpha = 0.8, fillstyle ='full',**markerstyle1)
                # plt.plot([lp[::2]], [idmax],alpha = 0.8, fillstyle ='full',**markerstyle2)
                #
                # markerstyle1 = dict(color=l_colors[idmax], marker='<',
                #     markersize=8)
                # markerstyle2 = dict(color=l_colors[idmin],  marker='<',
                #     markersize=8)
                #
                # plt.plot([lp[1::2]], [idmin],alpha = 0.8, fillstyle ='full',**markerstyle1)
                # plt.plot([lp[1::2]], [idmax],alpha = 0.8, fillstyle ='full',**markerstyle2)
        ticks = np.linspace(int(self.nodes[0]), int(self.nodes[-1]), 25)

        if label_nodes_by_letters == True and lnodes <= 26:
            plt.yticks(ticks, [ascii_uppercase[i] for i in ticks], fontname='Ubuntu',
                       fontsize=fontsize, color='#666699')
        elif self.node_to_label:

            plt.yticks(ticks, [self.node_to_label[i] for i in ticks], fontname='Ubuntu',
                       fontsize=fontsize,
                       color='#666699')
        else:
            plt.yticks(ticks, [self.nodes[i] for i in ticks], fontname='Ubuntu', fontsize=fontsize,
                       color='#666699')

        if xtimestamp == True:
            timestamps = np.linspace(int(self.times[0]), int(self.times[1]), 11)
            dates = [dt.datetime.fromtimestamp(ts) for ts in timestamps]
            xfmt = md.DateFormatter('%Y-%m-%d %H:%M:%S')
            ax = plt.gca()
            ax.xaxis.set_major_formatter(xfmt)
            plt.xticks(ticks=np.linspace(int(self.times[0]), int(self.times[1]), 11),
                       labels=dates,
                       rotation=30, fontname='Ubuntu', fontsize=fontsize, color='#476b6b')
        else:
            plt.xticks(ticks=np.linspace(int(self.times[0]), int(self.times[1]), 11),
                       fontname='Ubuntu', fontsize=fontsize, color='#476b6b')

        plt.ylabel("Nodes", fontname='Ubuntu', fontsize=fontsize, color='#666699')
        plt.xlabel("t", fontname='Ubuntu', fontsize=fontsize, color='#476b6b')
        plt.ylim((-1, lnodes))
        if title:
            plt.title(title, fontname='Ubuntu', fontsize=fontsize)
        # else:
        #     plt.title("Stream Graph", fontname='Ubuntu', fontsize=14)
        # Get rid of the frame
        for place, spine in plt.gca().spines.items():
            if place != 'bottom':
                spine.set_visible(False)
            else:
                spine.set_bounds(self.times[0], self.times[1])
                spine.set_color('#476b6b')
        plt.tick_params(top=False, bottom=True, right=False, left=True, labelleft=True, labelbottom=True)
        if animated:

            interactions_times = set()
            for k, v in self.get_interaction_times().items():
                for i in v:
                    interactions_times.add(i)
            line, = ax.plot([], [], color='#476b6b', lw=8, alpha=0.3)

            def init():
                line.set_data([], [])
                return line,

            def animate(t):
                y = np.linspace(0, len(self.nodes), 1000)
                x = [t]
                ax.set_xlabel("Temps : {0}".format(t))
                line.set_data(x, y)
                return line,

            # anim = animation.FuncAnimation(fig, animate, init_func=init,
            #                                frames=np.linspace(self.times[0], self.times[1], 100), interval=500,
            #                                blit=True, repeat=repeat)
            anim = animation.FuncAnimation(fig, animate, init_func=init,
                                           frames=list(sorted(interactions_times)), interval=1000,
                                           blit=True, repeat=repeat)
            # anim.save('stream_graph_animation.mp4',
            #           writer=Writer(bitrate = -1,fps =1),
            #           dpi = 120)
            # plt.show()
            return anim
        else:
            plt.tight_layout()
            # plt.show()
        if saving_path:
            fig.savefig(saving_path + "." + format, format=format)
        if clusters:
            return c_map_cluster
        return fig

    def plot_dict_clusters(self, dict_clusters, title=None, links=True, saving_path=None,
                           label_nodes_by_letters=False, format='pdf', cmap='plasma', xtimestamp=False, fontsize=26):
        """
        Display in an elegant way a small stream graph
        We can also specify a path to save the animated of ploted stream graph
        :param animated: If we want an animation with cursor which move according to the time
        :return: A matplotlib plot
        """
        lnodes = len(self.nodes)
        # c_map_cluster = get_cmap(max(dict_clusters.keys()))
        norm_colors = plt.Normalize(vmin=min(dict_clusters.keys()), vmax=max(dict_clusters.keys()))
        c_map_cluster = cm.get_cmap(cmap, len(dict_clusters.keys()))
        fig = plt.figure()
        list_legend = []

        # Plot Nodes
        for p, t in zip(range(lnodes), self.node_presence):
            color_n = '#666699'  # c_map(p)
            plt.hlines([p] * (len(t) // 2), xmin=t[::2], linewidth=2,
                       xmax=t[1::2], colors=color_n, alpha=0.9)

        # # Deal with labels
        # if self.node_to_label and :
        #     label_to_node = {v:k for k,v in self.node_to_label.items()}
        #     dict_clusters = {k:[(t0,t1,label_to_node[n]) for t0,t1,n in v] for k,v in dict_clusters.items()}

        # Plot Clusters
        for d in sorted(dict_clusters, reverse=True):
            current_color = c_map_cluster(norm_colors(d))
            if title:
                list_legend.append(mpatch.Patch(color=current_color, label=str(round(d, 2)), ))
            for n in dict_clusters[d]:
                plt.hlines(n[2], xmin=n[0], linewidth=5,
                           xmax=n[1], color=current_color, alpha=1)
        # Dirty trick to scale legend
        while len(list_legend) > 10:
            list_legend = list_legend[:len(list_legend) - 1:2] + [list_legend[-1]]

        #  Plot links
        if links:
            for k, lp in zip(self.links, self.link_presence):
                id1 = self.nodes.index(k[0])
                id2 = self.nodes.index(k[1])
                idmax = max(id1, id2)
                idmin = min(id1, id2)
                plt.vlines(lp[::2], ymin=idmin, ymax=idmax,
                           linewidth=1.5, alpha=0.15)
                eps = random.choice([1, -1]) * (random.random() / 5)
                plt.hlines([(idmax + idmin) / 2 + eps] * (len(lp) // 2), xmin=lp[::2], xmax=lp[1::2],
                           colors='k', linewidth=1.7, alpha=0.5)

                # ARRIVALS
                plt.plot([lp[::2]], [idmin], color='#004d00', marker='^', alpha=1, markersize=7)
                plt.plot([lp[::2]], [idmax], color='#004d00', marker='v', alpha=1, markersize=7)

                # DEPARTURE
                # plt.plot([lp[1::2]],[idmin,],color='#800000',marker='v',alpha =1,markersize=7)
                # plt.plot([lp[1::2]],[idmax],color='#800000',marker='^',alpha =1,markersize=7)
        if label_nodes_by_letters == True and lnodes <= 26:
            plt.yticks(range(lnodes), ascii_uppercase[0:lnodes], fontname='Ubuntu',
                       fontsize=fontsize, color='#666699')
        elif self.node_to_label:
            plt.yticks(range(lnodes), [self.node_to_label[i] for i in range(lnodes)], fontname='Ubuntu',
                       fontsize=fontsize,
                       color='#666699')
        else:
            plt.yticks(range(lnodes), range(lnodes), fontname='Ubuntu', fontsize=fontsize,
                       color='#666699')

        if xtimestamp == True:
            timestamps = np.linspace(int(self.times[0]), int(self.times[1]), 11)
            dates = [dt.datetime.fromtimestamp(ts) for ts in timestamps]
            xfmt = md.DateFormatter('%Y-%m-%d %H:%M:%S')
            ax = plt.gca()
            ax.xaxis.set_major_formatter(xfmt)
            plt.xticks(ticks=np.linspace(int(self.times[0]), int(self.times[1]), 11),
                       labels=dates,
                       rotation=30, fontname='Ubuntu', fontsize=12, color='#476b6b')
        else:
            plt.xticks(ticks=np.linspace(int(self.times[0]), int(self.times[1]), 11),
                       fontname='Ubuntu', fontsize=fontsize, color='#476b6b')
        plt.ylabel("Nodes", fontname='Ubuntu', fontsize=fontsize, color='#666699')
        plt.xlabel("t", fontname='Ubuntu', fontsize=fontsize, color='#476b6b')
        plt.legend(handles=list_legend, loc='upper right')
        # plt.colorbar(c_map_cluster)
        if title:
            plt.title(title, fontname='Ubuntu', fontsize=fontsize)
        for place, spine in plt.gca().spines.items():
            if place != 'bottom':
                spine.set_visible(False)
            else:
                spine.set_bounds(self.times[0], self.times[1])
                spine.set_color('#476b6b')

        plt.tick_params(top=False, bottom=True, right=False, left=True, labelleft=True, labelbottom=True)
        plt.tight_layout()
        if saving_path:
            fig.savefig(saving_path + "." + format, format=format)
        return fig

    def plot_3d_stream(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        lnodes = len(self.nodes)
        c_map = get_cmap(lnodes)
        T = self.times
        ax.set_xlim(0, lnodes)
        ax.set_ylim(T[0], T[1])
        z_max = max([max([j - i for i, j in zip(t[::2], t[1::2])]) for t in self.link_presence])
        ax.set_zlim(0, z_max)
        y = [T[0], T[1]]

        for p, t in zip(range(lnodes), self.node_presence):
            clr = c_map(p)
            # dotted line all the time corresponding to the existence of the node,
            # not his presence
            if len(self.nodes) < 350:
                ax.plot(xs=[p, p], ys=y, zs=0, zdir='z', linestyle='--', linewidth=0.7,
                        color=clr, alpha=0.1)
            for i, j in zip(t[::2], t[1::2]):
                # Line corresponding to the node presence
                ax.plot(xs=[p, p], ys=[i, j], zs=0, zdir='z', linewidth=1.1,
                        color=clr, alpha=0.9)
        for k, t in zip(self.links, self.link_presence):
            id1 = self.nodes.index(k[0])
            id2 = self.nodes.index(k[1])
            idmax = max(id1, id2)
            idmin = min(id1, id2)
            for i, j in zip(t[::2], t[1::2]):

                # eps = random.choice([1, -1]) * (random.random() / 5)
                id_tmp = (idmax + idmin) / 2  # + eps
                # Line corresponding to the duration of the link
                ax.plot(xs=[id_tmp, id_tmp], ys=[i, i], zs=[0, j - i], zdir='z',
                        color='#152686', linewidth=7, alpha=0.25)
                if len(self.links) < 150:
                    # Line corresponding to the link between two nodes
                    ax.plot(xs=[idmin, idmax], ys=[i, i], zs=0, zdir='z',
                            linewidth=1.5, alpha=0.15, color='#BDBDBD')
                    # Lines corresponding of both ending of the link to each node
                    ax.plot(xs=[idmax], ys=[i], zs=0, zdir='z', color='k', marker='o', alpha=0.2)
                    ax.plot(xs=[idmin], ys=[i], zs=0, zdir='z', color='k', marker='o', alpha=0.2)
        ax.set_title('3d representation of a Stream Graph')
        ax.set_xticklabels(self.nodes)
        ax.set_xlabel('Nodes')
        ax.set_ylabel('Time')
        ax.set_zlabel('Link\'s duration ')
        ax.view_init(elev=20., azim=-35)
        return fig

        # plt.show()

    def continous_to_discrete(self):
        """
        Discretize the intervals in node_presence and links
        :return: Nothing
        """
        self.times_d = discretize_time_interval(self.times_c, self.bins)
        self.node_presence_d = {}
        # Seee dict comprehension for efficiency
        for i in self.node_presence_c:
            tmp = []
            for t in self.node_presence_c[i]:
                tmp.append(discretize_time_interval(t, self.bins))
            self.node_presence_d[i] = np.concatenate(tmp)
        self.links_d = {}
        # Maybe make a funcion, it's the same ^v
        for i in self.links_c:
            tmp = []
            for t in self.links_c[i]:
                tmp.append(discretize_time_interval(t, self.bins))
            self.links_d[i] = np.concatenate(tmp)

    # Time must be divided by the same number as the other intervals
    # ( in node_presence and in link)
    def discrete_to_continuous(self):
        """
        Transform discrete intervals in node_presence and links into continous intervals
        :return:
        """
        self.times_c = [min(self.times_d), max(self.times_c)]
        self.node_presence_c = {}
        for i in self.node_presence_d:
            self.node_presence_c[i] = \
                get_bounds_for_evenly_spaced_intervals(self.node_presence_d[i],
                                                       self.bins)
        self.links_c = {}
        for i in self.links_d:
            self.links_c[i] = \
                get_bounds_for_evenly_spaced_intervals(self.links_d[i],
                                                       self.bins)

    def plot_aggregated_graph(self):
        fig = plt.figure()
        ax = plt.axes()
        G_glob = self.aggregated_graph(to_networkx=True)
        pos = nx.circular_layout(G_glob)
        l = len(self.nodes)
        c_map = get_cmap(l)
        dict_node_colors = {n: c_map(p) for n, p in zip(G_glob.nodes(), range(l))}

        nx.draw_networkx_labels(G_glob, pos, font_size=28, ax=ax)
        nodes = nx.draw_networkx_nodes(G_glob, pos, node_size=900, node_color='w', alpha=0, ax=ax)
        G = G_glob
        node_c = [dict_node_colors[n] for n in G.nodes()]
        nodes = nx.draw_networkx_nodes(G, pos, node_size=900,
                                       node_color=node_c, alpha=0.5, ax=ax)
        edges = nx.draw_networkx_edges(G, pos, edge_color='#2d5986',
                                       alpha=0.3, width=6, ax=ax)
        for place, spine in plt.gca().spines.items():
            spine.set_visible(False)
        plt.tick_params(top=False, bottom=False, right=False, left=False, labelbottom=False, labelleft=False)
        return fig

    def plot_induced_graph(self, time=0):

        fig = plt.figure()
        ax = plt.axes()
        G_glob = self.aggregated_graph(to_networkx=True)
        pos = nx.circular_layout(G_glob)
        l = len(self.nodes)
        c_map = get_cmap(l)
        dict_node_colors = {n: c_map(p) for n, p in zip(G_glob.nodes(), range(l))}

        nx.draw_networkx_labels(G_glob, pos, font_size=28, ax=ax)
        nodes = nx.draw_networkx_nodes(G_glob, pos, node_size=900, node_color='w', alpha=0, ax=ax)

        G = nx.from_dict_of_lists(self.instant_graph(time))
        node_c = [dict_node_colors[n] for n in G.nodes()]
        nodes = nx.draw_networkx_nodes(G, pos, node_size=900,
                                       node_color=node_c, alpha=0.5, ax=ax)
        edges = nx.draw_networkx_edges(G, pos, edge_color='#2d5986',
                                       alpha=0.3, width=6, ax=ax)
        ax.set_xlabel("Time : {0}".format(time), fontsize=32)
        for place, spine in plt.gca().spines.items():
            spine.set_visible(False)
        plt.tick_params(top=False, bottom=False, right=False, left=False, labelbottom=False, labelleft=False)
        return fig

    def plot_induced_graphs(self, repeat=False):
        """
        Draw and present an animation of the stream graph from the induced graphs
        we can save this animation
        :return:
        """
        fig = plt.figure()
        ax = plt.axes()
        G_glob = self.aggregated_graph(to_networkx=True)
        pos = nx.circular_layout(G_glob)
        l = len(self.nodes)
        c_map = get_cmap(l)
        dict_node_colors = {n: c_map(p) for n, p in zip(G_glob.nodes(), range(l))}

        def init():
            nx.draw_networkx_labels(G_glob, pos, font_size=28, ax=ax)
            nodes = nx.draw_networkx_nodes(G_glob, pos, node_size=900, node_color='w', alpha=0, ax=ax)
            return nodes,  # labels

        def update(t):
            G = self.instant_graph(t, to_networkx=True)
            node_c = [dict_node_colors[n] for n in G.nodes()]
            nodes = nx.draw_networkx_nodes(G, pos, node_size=900,
                                           node_color=node_c, alpha=0.5, ax=ax)
            edges = nx.draw_networkx_edges(G, pos, edge_color='#2d5986',
                                           alpha=0.3, width=6, ax=ax)
            ax.set_xlabel("Time : {0}".format(t), fontsize=32)
            if not nodes:
                nodes = nx.draw_networkx_nodes(G_glob, pos, node_size=900, node_color='w', alpha=0, ax=ax)
            if not G.edges():
                return nodes,
            return nodes, edges,  # labels

        interactions_times = set()
        for k, v in self.get_interaction_times().items():
            for i in v:
                interactions_times.add(i)

        for place, spine in plt.gca().spines.items():
            spine.set_visible(False)

        plt.tick_params(top=False, bottom=False, right=False, left=False, labelbottom=False, labelleft=False)

        anim = animation.FuncAnimation(fig, update, init_func=init,
                                       frames=list(sorted(interactions_times)), interval=1000,
                                       blit=True, repeat=repeat)
        # anim = animation.FuncAnimation(fig, update, init_func=init,
        #                                frames=np.linspace(self.times[0], self.times[1], 100), interval=1000,
        #                                blit=True, repeat=repeat)
        # anim.save("induce_graphs_anim.mp4",
        #           writer=Writer(bitrate=-1,fps =1),dpi = 120)
        return anim

    def set_card_W(self):
        self.card_W = sum([t1 - t0 for i in self.node_presence
                           for t0, t1 in zip(i[::2], i[1::2])])

    def set_card_E(self):
        self.card_E = sum([t1 - t0 for i in self.link_presence
                           for t0, t1 in zip(i[::2], i[1::2])])

    def duration(self):
        return self.times[1] - self.times[0]

    def surface(self):
        return sum([t1 - t0 for v in self.link_presence for t0, t1 in zip(v[::2], v[1::2])])

    def coverage(self):
        T = self.duration()
        cov = self.card_W / (len(self.nodes) * T)
        return cov

    def nb_nodes(self):
        """
        :return: The number nodes according to their presence in the stream graph
        """
        T = self.times[1] - self.times[0]
        nb_nodes = self.card_W / T
        # for i in self.node_presence:
        #     nb_nodes+=sum(map((lambda ta,tb: tb-ta),i[::2],i[1::2]))/T
        # nb_nodes=0
        # for i in self.node_presence:
        #     nb_nodes += reduce((lambda a, b: b - a), i)
        return nb_nodes

    def nb_links(self):
        """
        :return: the number of links according to their presence in the stream graph
        """
        T = self.times[1] - self.times[0]
        nb_links = self.card_E / T
        return nb_links

    def nodes_over_time(self):
        '''
        Compute the number of nodes over time in the Stream Graph.
        :return: A dictionnary event time to the corresponding number of nodes
        '''
        t_to_nb_nodes = {}
        eps = 10 ** -7
        cnt_nodes = 0
        t_last_event = None
        for e in self.ordered_events():
            c = e[0]
            t = e[1]
            if c == 2:
                _, t0, t1, n = e
                if t == t_last_event:
                    cnt_nodes += 1
                    t_to_nb_nodes[t] = cnt_nodes
                else:
                    t_to_nb_nodes[t - eps] = cnt_nodes
                    cnt_nodes += 1
                    t_to_nb_nodes[t] = cnt_nodes
                    t_last_event = t
            elif c == -2:
                _, t1, n = e
                if t == t_last_event:
                    cnt_nodes -= 1
                    t_to_nb_nodes[t] = cnt_nodes
                else:
                    t_to_nb_nodes[t - eps] = cnt_nodes
                    cnt_nodes -= 1
                    t_to_nb_nodes[t] = cnt_nodes
                    t_last_event = t
        return t_to_nb_nodes

    def mean_degree_over_time(self):
        '''
        Compute the mean degree over time in the Stream Graph.
        :return: A dictionnary event time to the corresponding number of nodes
        '''
        t_to_mean_degree = {}
        eps = 10 ** -7
        cnt_nodes = 0
        cnt_global_degree = 0
        t_last_event = None
        for e in self.ordered_events():
            c = e[0]
            t = e[1]
            if c == 1:
                _, t0, t1, u, v = e
                if t == t_last_event:
                    cnt_global_degree += 2
                    t_to_mean_degree[t] = cnt_global_degree / cnt_nodes
                else:
                    t_to_mean_degree[t - eps] = cnt_global_degree / cnt_nodes
                    cnt_global_degree += 2
                    t_to_mean_degree[t] = cnt_global_degree / cnt_nodes
                    t_last_event = t
            elif c == -1:
                _, t1, u, v = e
                if t == t_last_event:
                    cnt_global_degree -= 2
                    t_to_mean_degree[t] = cnt_global_degree / cnt_nodes
                else:
                    t_to_mean_degree[t - eps] = cnt_global_degree / cnt_nodes
                    cnt_global_degree -= 2
                    t_to_mean_degree[t] = cnt_global_degree / cnt_nodes
                    t_last_event = t

            elif c == 2:
                _, t0, t1, n = e
                if t == t_last_event:
                    cnt_nodes += 1
                    if cnt_nodes > 0:
                        t_to_mean_degree[t] = cnt_global_degree / cnt_nodes
                    else:
                        t_to_mean_degree[t] = 0
                else:
                    if cnt_nodes > 0:
                        t_to_mean_degree[t - eps] = cnt_global_degree / cnt_nodes
                    else:
                        t_to_mean_degree[t - eps] = 0
                    cnt_nodes += 1
                    t_to_mean_degree[t] = cnt_global_degree / cnt_nodes
                    t_last_event = t
            elif c == -2:
                _, t1, n = e
                if t == t_last_event:
                    cnt_nodes -= 1
                    if cnt_nodes > 0:
                        t_to_mean_degree[t] = cnt_global_degree / cnt_nodes
                    else:
                        t_to_mean_degree[t] = 0
                else:
                    t_to_mean_degree[t - eps] = cnt_global_degree / cnt_nodes
                    cnt_nodes -= 1
                    if cnt_nodes > 0:
                        t_to_mean_degree[t] = cnt_global_degree / cnt_nodes
                    else:
                        t_to_mean_degree[t] = 0
                    t_last_event = t
        return t_to_mean_degree

    def node_weight_at_t(self, t):
        V = len(self.nodes)
        # print("V :",V)
        Vt = sum([1 for i in self.node_presence for ta, tb in zip(i[::2], i[1::2]) if ta <= t and t <= tb])
        # print("V_"+str(t)+" :",Vt)
        return Vt / V

    def link_weight_at_t(self, t):
        l = len(self.nodes)
        E = (l * (l - 1)) / 2
        # print("E :",E)
        Et = sum([1 for i in self.link_presence for ta, tb in zip(i[::2], i[1::2]) if ta <= t and t <= tb])
        # print("E_"+str(t)+" :",Et)
        return Et / E

    def plot_link_weight(self):
        fig = plt.figure()
        T = self.event_times()
        # add epsilon to capture information just before an event and just after
        epsilon = 10 ** -9
        resampled_T = []
        for t in T:
            resampled_T.append(t - epsilon)
            resampled_T.append(t)
            resampled_T.append(t + epsilon)
        lt = [self.link_weight_at_t(t) for t in resampled_T]
        plt.plot(resampled_T, lt)
        plt.xlabel("Time")
        plt.ylabel("Link weight")
        plt.title("Link weight through time")
        plt.tight_layout()
        plt.show()
        return

    def node_weight_on_I(self, I):
        V = len(self.nodes)
        # print("V :",V)
        Vt = [sum([1 for i in self.node_presence for ta, tb in zip(i[::2], i[1::2]) if ta <= t <= tb]) / V
              for t in I]
        # print("V_"+str(t)+" :",Vt)
        return Vt

    def link_weight_on_I(self, I):
        l = len(self.nodes)
        E = (l * (l - 1)) / 2
        # print("E :",E)
        Et = [sum([1 for i in self.link_presence for ta, tb in zip(i[::2], i[1::2]) if ta <= t <= tb]) / E
              for t in I]
        # print("E_"+str(t)+" :",Et)
        return Et

    def plot_node_weight(self):
        fig = plt.figure()
        T = self.event_times()
        # add delta to capture information just before an event and just after
        epsilon = 10 ** -9
        resampled_T = []
        for t in T:
            resampled_T.append(t - epsilon)
            resampled_T.append(t)
            resampled_T.append(t + epsilon)
        kt = [self.node_weight_at_t(t) for t in resampled_T]
        # kt = self.node_weight_on_I(T)
        plt.plot(resampled_T, kt)
        plt.xlabel("Time")
        plt.ylabel("Node weight")
        plt.title(" Node weight through time")
        plt.tight_layout()
        plt.show()
        return

    def node_duration(self):
        nd = self.card_W / len(self.nodes)
        return nd

    def link_duration(self):
        l = len(self.nodes)
        possible_pairs = (l * (l - 1)) / 2
        return self.card_E / possible_pairs

    def get_sum_intersection(self):
        l = len(self.nodes)
        sum_node_intersection = sum([min(a1, b1) - max(a0, b0)
                                     for i1 in range(l)
                                     for i2 in range(i1 + 1, l)
                                     for a0, a1 in zip(self.node_presence[i1][::2], self.node_presence[i1][1::2])
                                     for b0, b1 in zip(self.node_presence[i2][::2], self.node_presence[i2][1::2])
                                     if (a0 <= b1 and b0 <= a1)
                                     ])
        return sum_node_intersection

    def get_sum_union(self):
        l = len(self.nodes)
        sum_node_union = sum([max(a1, b1) - min(a0, b0)
                              if (a0 <= b1 and b0 <= a1)
                              else (b1 - b0) + (a1 - a0)
                              for i1 in range(l)
                              for i2 in range(i1 + 1, l)
                              for a0, a1 in zip(self.node_presence[i1][::2], self.node_presence[i1][1::2])
                              for b0, b1 in zip(self.node_presence[i2][::2], self.node_presence[i2][1::2])
                              ])
        return sum_node_union

    def uniformity(self):
        sum_node_intersection = self.get_sum_intersection()
        sum_node_union = self.get_sum_union()
        return sum_node_intersection / sum_node_union

    def density(self):
        """
        :return: the density of the stream graph (aka the probability
        if we randomly chose two nodes that they are linked together)
        """
        l = len(self.nodes)
        sum_links = self.card_E
        sum_node_intersection = self.get_sum_intersection()
        return sum_links / sum_node_intersection

    def link_densities(self):
        densities = Counter()
        for l, lp in zip(self.links, self.link_presence):
            n1 = l[0]
            n2 = l[1]
            Tuv = sum([t1 - t0 for t0, t1 in zip(lp[::2], lp[1::2])])
            #  Pour que ça marche il faut que les noeuds soient représentés par leurs indices dans la liste
            intersection_Tu_Tv = sum([min(a1, b1) - max(a0, b0)
                                      for a0, a1 in zip(self.node_presence[n1][::2], self.node_presence[n1][1::2])
                                      for b0, b1 in zip(self.node_presence[n2][::2], self.node_presence[n2][1::2])
                                      if (a0 <= b1 and b0 <= a1)
                                      ])
            densities[l] = Tuv / intersection_Tu_Tv
        return densities

    def node_densities(self):
        densities = Counter()
        for n, np in zip(self.nodes, self.node_presence):
            sum_links = sum([t1 - t0
                             for l, lp in zip(self.links, self.link_presence)
                             if (l[0] == n or l[1] == n)
                             for t0, t1 in zip(lp[::2], lp[1::2])
                             ])

            sum_node_intersection = sum([min(a1, b1) - max(a0, b0)
                                         for n2, np2 in zip(self.nodes, self.node_presence)
                                         if (n2 != n)
                                         for a0, a1 in zip(np[::2], np[1::2])
                                         for b0, b1 in zip(np2[::2], np2[1::2])
                                         if (a0 <= b1 and b0 <= a1)
                                         ])
            densities[n] = sum_links / sum_node_intersection
        return densities

    def neighborhood(self):
        """
        :return: the neighborhood of each node globally in the form of a dictionary
        each node --> their neighbors with the associated interval of time (of the link)
        """
        N = {n: {} for n in self.nodes}
        for l, lp in zip(self.links, self.link_presence):
            n1 = l[0]
            n2 = l[1]
            N[n2][n1] = lp
            N[n1][n2] = lp
        return N

    def neighborhood_from_node(self, node):
        """
        :return: the neighborhood of a specific node in the form of a dictionary
        node <-> their neighbor with the associated interval of time
        """
        N = {}
        for l, lp in zip(self.links, self.link_presence):
            n1 = l[0]
            n2 = l[1]
            if n1 == node:
                N[n2] = lp
            if n2 == node:
                N[n1] = lp
        return N

    def degrees(self):
        """
        :return: A dictionary with to each nodes their degree
        """
        T = self.times[1] - self.times[0]
        degrees = Counter()
        for l, lp in zip(self.links, self.link_presence):
            n1 = l[0]
            n2 = l[1]
            s_d = sum([tb - ta for ta, tb in zip(lp[::2], lp[1::2])]) / T
            # For our dataset, due to splitting links between windows,
            # check that the degree is not < to something not coherent
            if self.bins and self.bins / T > s_d:
                degrees[n1] += self.bins / T
                degrees[n2] += self.bins / T
            else:
                degrees[n1] += s_d
                degrees[n2] += s_d
        return degrees

    def nb_neighbors(self):
        nb_neighbors = Counter()
        for l in self.links:
            nb_neighbors[l[0]] += 1
            nb_neighbors[l[1]] += 1
        return nb_neighbors

    def sum_degree_neighborhood(self, degrees=None):
        if not degrees:
            degrees = self.degrees()
        sum_degree_neighborhood = Counter()
        for l, lp in zip(self.links, self.link_presence):
            n1 = l[0]
            n2 = l[1]
            sum_degree_neighborhood[n1] += degrees[n2]
            sum_degree_neighborhood[n2] += degrees[n1]
        return sum_degree_neighborhood

    def mean_degree_neighborhood(self, degrees=None):
        S = self.sum_degree_neighborhood(degrees)
        N = self.nb_neighbors()
        mean_degree_neighborhood = Counter()
        # std_degree_neighborhood = Counter()
        for n in self.nodes:
            mean_degree_neighborhood[n] = S[n] / N[n]
        return mean_degree_neighborhood  # ,std_degree_neighborhood

    def expected_node_degrees(self, d=None):
        T = self.times[1] - self.times[0]
        expected_degrees = Counter()
        if not d:
            d = self.degrees()
        for n, np in zip(self.nodes, self.node_presence):
            expected_degrees[n] = (d[n] * T) / sum([t1 - t0 for t0, t1 in zip(np[::2], np[1::2])])
        return expected_degrees

    def average_degree(self, degrees=None, nb_nodes=None):
        """
        :return: the average degree of a stream graph
        """
        if not degrees:
            degrees = self.degrees()
        if not nb_nodes:
            nb_nodes = self.nb_nodes()
        T = self.times[1] - self.times[0]
        d_bar = sum([degrees[n] * sum([tb - ta for ta, tb in zip(np[::2], np[1::2])])
                     for n, np in enumerate(self.node_presence)
                     ]) / (nb_nodes * T)
        return d_bar

    def stream_graph_degree(self):
        T = self.times[1] - self.times[0]
        return self.card_E / (len(self.nodes) * T)

    def expected_stream_graph_degree(self):
        return 2 * self.card_E / self.card_W

    def clustering_coefficient(self):
        """
        :return: A dictionary with for each node his clustering coefficient
        """
        cc = Counter()
        N = self.neighborhood()
        for current_node in N:
            current_neighbors = N[current_node]
            # print("current node :", current_node)
            # print("current neighbors :", current_neighbors)
            l = len(current_neighbors)
            list_neighbors = list(current_neighbors)
            # print("l :", l)
            # print("list_nodes :", list_neighbors)
            sum_link_intersection = sum([min(a1, b1) - max(a0, b0)
                                         for i1 in range(l)
                                         for i2 in range(i1 + 1, l)
                                         for a0, a1 in zip(current_neighbors[list_neighbors[i1]][::2],
                                                           current_neighbors[list_neighbors[i1]][1::2])
                                         for b0, b1 in zip(current_neighbors[list_neighbors[i2]][::2],
                                                           current_neighbors[list_neighbors[i2]][1::2])
                                         if (b1 >= a0 and a1 >= b0)
                                         ])
            # print("sum_link_intersection :", sum_link_intersection)
            sum_triplet_intersection = sum([
                min(a1, b1, c1) - max(a0, b0, c0)
                for i1 in range(l)
                for i2 in range(i1 + 1, l)
                for a0, a1 in zip(current_neighbors[list_neighbors[i1]][::2],
                                  current_neighbors[list_neighbors[i1]][1::2])
                for b0, b1 in zip(current_neighbors[list_neighbors[i2]][::2],
                                  current_neighbors[list_neighbors[i2]][1::2])
                if list_neighbors[i2] in N[list_neighbors[i1]]
                for c0, c1 in zip(N[list_neighbors[i1]][list_neighbors[i2]][::2],
                                  N[list_neighbors[i1]][list_neighbors[i2]][1::2])
                if (a1 >= b0 and a1 >= c0 and b1 >= a0 and b1 >= c0 and c1 >= b0 and c1 >= a0)
            ])
            # Last index on neighbors can be inverted if u is a neighbor of v,
            #  v is a neighbor of u
            # print("sum_triplet_intersection :", sum_triplet_intersection)
            if sum_link_intersection != 0:
                cc[current_node] = sum_triplet_intersection / sum_link_intersection
            else:
                cc[current_node] = 0
        return cc

    def average_clustering(self, cc=None):
        """
        :return: The avaerage clustering coefficient of a stream graph
        """
        if not cc:
            cc = self.clustering_coefficient()
        T = self.times[1] - self.times[0]
        cc_bar = sum([coef * sum([tb - ta for ta, tb in zip(i[::2], i[1::2])])
                      for i, coef in zip(self.node_presence, cc.values())]) / (self.nb_nodes() * T)
        return cc_bar

    def induced_line_stream(self):
        """
        :return: The induced line stream (which is a stream graph too) corresponding to the stream graph
        """
        induced_nodes = self.links
        induced_links = []
        link_presence = []
        for n1 in induced_nodes:
            for n2 in induced_nodes:
                if n1 == n2:
                    continue
                if n1[0] in n2 or n1[1] in n2:
                    if (n1, n2) in induced_links or (n2, n1) in induced_links:
                        continue
                    id1 = self.links.index(n1)
                    id2 = self.links.index(n2)
                    intersection_time = []
                    for t1a, t1b in zip(self.link_presence[id1][::2], self.link_presence[id1][1::2]):
                        for t2a, t2b in zip(self.link_presence[id2][::2], self.link_presence[id2][1::2]):
                            if (min(t1b, t2b) - max(t1a, t2a)) > 0:
                                intersection_time += [max(t1a, t2a)]
                                intersection_time += [min(t1b, t2b)]
                    # A=np.minimum(self.link_presence[id1][1::2],self.link_presence[id2][1::2])\
                    #   -np.maximum(self.link_presence[id1][::2],self.link_presence[id2][::2])
                    if intersection_time:
                        # print("n1 :", n1, "n2 :", n2)
                        induced_links.append((n1, n2))
                        link_presence.append(intersection_time)
                        # print("induced links :", induced_links)
        LS = stream_graph(times=self.times,
                          nodes=induced_nodes,
                          node_presence=self.link_presence,
                          links=induced_links,
                          link_presence=link_presence
                          )
        # print("LS.nodes :", LS.nodes)
        # print("LS.node_presence :", LS.node_presence)
        # print("LS links :", LS.links)
        # print("LS.link_presence :", LS.link_presence)
        return LS

    def check_node_presence(self, A, t):
        """
        :param A: A node
        :param t: An instant
        :return: True if the node is present False otherwise
        """
        for i in self.node_presence_c[A]:
            if t >= i[0] and t <= i[1]:
                return True
        return False

    def check_link_presence(self, L, t):
        """
        :param L: A link
        :param t: An instant
        :return: True if the link is present False otherwise
        """
        for i in self.links_c[L]:
            if t >= i[0] and t <= i[1]:
                return True
        return False

    def get_neighborhood_time_t(self, A, t):
        """
        :param A: A node
        :param t: An instant
        :return: a dictionary : the node <-> his neighbors at time t
        """
        N = {A: []}
        if self.check_node_presence(A, t):
            for l in self.links_c:
                n1 = l[0]
                n2 = l[1]
                if A != n1 and A != n2:
                    continue
                for i in self.links_c[l]:
                    if i[0] <= t and i[1] >= t:
                        if A == n1:
                            N[A].append(n2)
                        if A == n2:
                            N[A].append(n1)
            return N
        return None

    def number_of_link_per_node(self, dataframe=False):
        d = {n: 0 for n in self.nodes}
        for l, lp in zip(self.links, self.link_presence):
            s = len(lp) / 2
            d[l[0]] += s
            d[l[1]] += s
        print("Number of links per node DONE !")
        # print(d)
        if dataframe:
            D = pd.DataFrame.from_dict(d, orient='index', dtype=int)
            # D.plot.bar()
            # plt.show()
            return D
        else:
            return d

    def number_of_pair_per_link(self, dataframe=False):
        d = {n: 0 for n in self.links}
        for l, lp in zip(self.links, self.link_presence):
            s = len(lp) / 2
            d[l] += s
        print("Number of pair per link DONE !")
        # print(d)
        if dataframe:
            D = pd.DataFrame.from_dict(d, orient='index', dtype=int)
            # D.plot.bar()
            # plt.show()
            return D
        else:
            return d

    def transform_links_label_to_int(self, index=False):
        dict_node_label_2_int = {}
        dict_int_2_node_label = {}
        # print("nodes 1:",self.nodes)
        # print("links 1:",self.links)
        for i, n in enumerate(self.nodes):
            dict_node_label_2_int[n] = i
            dict_int_2_node_label[i] = n
            self.nodes[i] = i
        for l, j in zip(self.links, range(len(self.links))):
            self.links[j] = (dict_node_label_2_int[l[0]], dict_node_label_2_int[l[1]])
        # print("nodes 2:",self.nodes)
        # print("links 2:",self.links)
        self.node_to_label = dict_int_2_node_label
        print("NODE LABEL :", self.node_to_label)

        if index:
            return dict_int_2_node_label, dict_node_label_2_int

    #####################################################
    #           Methods to suppress nodes or links      #
    #####################################################

    def remove_node(self, n, time_presence):
        """
        Remove the (node,time_presence) from the Stream Graph
        :param n:
        :return:
        """
        self.nodes.remove(n)
        del self.node_presence[n]

    def remove_link(self, l):
        """
        Remove the (link,time_presence) from the Stream Graph
        :param l:
        :return:
        """
        del self.link_presence[self.links.index(l)]
        self.links.remove(l)

    #########################################
    #       Methods connected components    #
    #########################################

    def segmented_neighborhood(self):
        """
        On veut des segments de noeuds comme clefs ainsi que dans les voisins (b,e,v) : [(t0,t1,(nt0,nt1,u))] :)
        :param links:
        :param node_presence:
        :return:
        """
        # See 'neighborhood_with_node_presence' in stream.py
        N = defaultdict(SortedSet)
        for l, lp in zip(self.links, self.link_presence):
            for t0, t1 in zip(lp[::2], lp[1::2]):
                nu, nv = self.get_according_node_presence_from_link(l, t0, t1)
                if nu and nv:
                    N[nu].add((t0, t1, nv))
                    N[nv].add((t0, t1, nu))
        return N

    def get_node_presence_from_interval(self, n, b, e):
        """
        Return the maximal temporal node corresponding to (b,e,n)
        :param n: node
        :param b: beginning of the interval (time)
        :param e: ending of the interval (time)
        :return: Maximal temporal node presence corresponding to (b,e,n) : (t0,t1)
        """
        for t0, t1 in zip(self.node_presence[n][::2], self.node_presence[n][1::2]):
            if b <= t0 and t1 <= e:
                return (t0, t1)
        return None

    def get_according_node_presence_from_link(self, l, t0, t1):
        """
        Get both extremities (which are maximal temporal nodes) of the link.
        :param l:
        :param t0: Beginning of the link
        :param t1: Ending of the link
        :return:
        """
        u, v = l[0], l[1]
        nu, nv = None, None
        for a, b in zip(self.node_presence[u][::2], self.node_presence[u][1::2]):
            if a <= t0 and t1 <= b:
                nu = (a, b, u)
                break
        for a, b in zip(self.node_presence[v][::2], self.node_presence[v][1::2]):
            if a <= t0 and t1 <= b:
                nv = (a, b, v)
                break
        return nu, nv

    def neighborhood_with_node_presence(self):
        """
        :return: the neighborhood of each node globally in the form of a dictionary
        each node --> their neighbors with their presences
        """
        N_presence = {}
        node_to_segmented_node = {}
        for e in self.ordered_events():
            c = e[0]
            if c == 2:
                _, t0, t1, n = e
                node_to_segmented_node[n] = (t0, t1, n)
                N_presence[(t0, t1, n)] = []
            elif c == 1:
                _, t0, t1, u, v = e
                u, v = node_to_segmented_node[u], node_to_segmented_node[v]
                N_presence[u].append(v)
                N_presence[v].append(u)

        # N_presence = {(t0, t1, n): [] for n, np in zip(self.nodes, self.node_presence)
        #               for t0, t1 in zip(np[::2], np[1::2])}
        # for l, lp in zip(self.links, self.link_presence):
        #     for t0, t1 in zip(lp[::2], lp[1::2]):
        #         np1, np2 = self.get_according_node_presence_from_link(l, t0, t1)
        #         if np1 and np2:
        #             N_presence[np1].append(np2)
        #             N_presence[np2].append(np1)
        return N_presence

    ####################################################################
    #                       Sort Functions                             #
    ####################################################################

    def node_position_by_arrival(self):
        """
        Sort nodes by their first time of arrival
        :return:
        """
        node_to_pos = {}
        time_of_arrival = {}
        for n, np in zip(self.nodes, self.node_presence):
            time_of_arrival[n] = min(np)
        cnt_pos = 0
        for n in time_of_arrival:
            node_to_pos[self.node_to_id[n]] = cnt_pos
            cnt_pos += 1
        return node_to_pos

    def node_position_by_increasing_degree(self, degrees_partition=None):
        """
        Sort nodes by their degree ; TODO : see S.nodes_degrees()
        :param degrees_partition:
        :return:
        """
        if degrees_partition is None:
            degrees_partition = self.degrees_partition()
            isolated_nodes = self.get_isolated_nodes()
            degrees_partition[0] = isolated_nodes
        node_to_pos = {}
        cnt_pos = 0
        for d in sorted(degrees_partition):
            sum_interval_degrees = Counter()
            for e in degrees_partition[d]:
                (t0, t1, n) = e
                sum_interval_degrees[n] += t1 - t0
            for n in sorted(sum_interval_degrees, key=sum_interval_degrees.get):
                node_to_pos[self.node_to_id[n]] = cnt_pos
                cnt_pos += 1
        return node_to_pos

    def node_position_by_peak_degree_arrival(self, degrees_partition=None):
        """
        Sort nodes increasingly by their maximal value of degree arrival time
        :return:
        """
        if degrees_partition is None:
            degrees_partition = self.degrees_partition()
            isolated_nodes = self.get_isolated_nodes()
            degrees_partition[0] = isolated_nodes
        node_to_pos = {}
        cnt_pos = len(self.nodes) - 1
        visited = set()
        peak_degree_arrival = Counter()
        for d in sorted(degrees_partition, reverse=True):
            # print("Degree :",d)
            for e in degrees_partition[d]:
                (t0, t1, n) = e
                if n not in visited:
                    if n in peak_degree_arrival:
                        peak_degree_arrival[n] = min(peak_degree_arrival[n], t0)
                    else:
                        peak_degree_arrival[n] = t0
            for n in peak_degree_arrival:
                visited.add(n)
        # print("Peak degree arrival :",peak_degree_arrival)
        for n in sorted(peak_degree_arrival, key=peak_degree_arrival.get):
            # print("n :",n,"cnt_pos :",cnt_pos)
            node_to_pos[self.node_to_id[n]] = cnt_pos
            cnt_pos -= 1
        return node_to_pos

    def link_position_by_duration(self):
        """
        sort links by their duration
        :return:
        """
        link_to_pos = {}
        link_duration = Counter()
        for l, lp in zip(self.links, self.link_presence):
            link_duration[l] = sum([t0 - t1 for t0, t1 in zip(lp[::2], lp[1::2])])
        cnt_pos = 0
        for l in sorted(link_duration, key=link_duration.get):
            link_to_pos[(self.node_to_id[l[0]], self.node_to_id[l[1]])] = cnt_pos
            cnt_pos += 1
        return link_to_pos

    def link_position_by_arrival(self):
        """
        sort links by their arrival
        :return:
        """
        link_to_pos = {}
        time_of_arrival = {}
        for l, lp in zip(self.links, self.link_presence):
            time_of_arrival[l] = min(lp)
        cnt_pos = 0
        for l in time_of_arrival:
            link_to_pos[(self.node_to_id[l[0]], self.node_to_id[l[1]])] = cnt_pos
            cnt_pos += 1
        return link_to_pos

    ####################################################################
    #       Induced Graphs and Substreams                              #
    ####################################################################

    def aggregated_graph(self, label=True, to_networkx=False, to_networkit=False):
        """
        Return the aggregated induced graph.

        :param label: True if we want the node's label in the adjacency list.
        :return: An adjacency list (networkx compatible)
        """
        if label:
            adjacency_list = {self.node_to_label[n]: [] for n in self.nodes}
            for l in self.links:
                n1 = self.node_to_label[l[0]]
                n2 = self.node_to_label[l[1]]
                adjacency_list[n1].append(n2)
                adjacency_list[n2].append(n1)
        else:
            adjacency_list = {n: [] for n in self.nodes}
            for l in self.links:
                n1 = l[0]
                n2 = l[1]
                adjacency_list[n1].append(n2)
                adjacency_list[n2].append(n1)
        if to_networkx:
            return nx.from_dict_of_lists(adjacency_list)
        if to_networkit:
            import networkit as nk
            G = nx.from_dict_of_lists(adjacency_list)
            return nk.nxadapter.nx2nk(G)
        return adjacency_list

    def instant_graph(self, t, label=True, to_networkx=False, to_networkit=False):
        """
        Return an adjacency list corresponding to the induced graph
        from the stream graph at a specified time *t*.

        :param t: Time instant
        :param label: True if we want the node's label in the adjacency list.
        :return: An adjacency list (networkx compatible)
        """
        adjacency_list = defaultdict(list)
        for i, l in enumerate(self.links):
            u = l[0]
            v = l[1]
            for t0, t1 in zip(self.link_presence[i][::2], self.link_presence[i][1::2]):
                if t0 <= t <= t1:
                    if label:
                        adjacency_list[self.node_to_label[u]].append(self.node_to_label[v])
                        adjacency_list[self.node_to_label[v]].append(self.node_to_label[u])
                    else:
                        adjacency_list[u].append(v)
                        adjacency_list[v].append(u)
        # Add isolated nodes:
        for n, np in zip(self.nodes, self.node_presence):
            for t0, t1 in zip(np[::2], np[1::2]):
                if t0 <= t <= t1:
                    if label:
                        adjacency_list[self.node_to_label[n]]  # Defaultdict functionality
                    else:
                        adjacency_list[n]  # Defaultdict functionality
        if to_networkx:
            return nx.from_dict_of_lists(adjacency_list)
        if to_networkit:
            G = nx.from_dict_of_lists(adjacency_list)
            return nk.nxadapter.nx2nk(G)
        return adjacency_list

    def induced_substream_by_nodes(self, nodes_list):
        """
        Return the sub stream induced by *nodes_list*.
        IMPORTANT: Only include links with both extremities in *nodes_list*. (see filter_links otherwise)

        :param nodes_list: A list of nodes (Ids or Labels)
        :return: A Stream Graph
        """

        # Check if element of nodes_list are labels or nodes id:
        if type(nodes_list[0]) == str:
            label_to_node = {v: k for k, v in self.node_to_label.items()}
            nodes_list = [label_to_node[n] for n in nodes_list]

        for n in nodes_list:
            if n not in self.nodes:
                raise ("Trying to filter Stream Graph by nodes that does not exist in Stream Graph")

        new_nodes = []
        new_node_presence = [[] for _ in nodes_list]
        nodes_to_new_nodes = defaultdict(lambda: len(nodes_to_new_nodes))
        new_nodes_to_label = {}
        new_nodes_to_id = {}
        new_links = []
        new_link_presence = []

        for n, np in zip(self.nodes, self.node_presence):
            if n in nodes_list:
                new_n = nodes_to_new_nodes[n]
                if self.node_to_label:
                    new_nodes_to_label[new_n] = self.node_to_label[n]
                if self.node_to_id:
                    new_nodes_to_id[new_n] = self.node_to_id[n]
                new_nodes.append(new_n)
                new_node_presence[new_n] = np
        new_nodes.sort()  #  To corresponds to emplacement in new_node_presence

        for l, lp in zip(self.links, self.link_presence):
            u, v = l
            if u in nodes_list and v in nodes_list:
                new_u = nodes_to_new_nodes[u]
                new_v = nodes_to_new_nodes[v]
                new_links.append((new_u, new_v))
                new_link_presence.append(lp)

        return stream_graph(times=self.times,
                            nodes=new_nodes,
                            node_to_label=new_nodes_to_label,
                            node_to_id=new_nodes_to_id,
                            node_presence=new_node_presence,
                            links=new_links,
                            link_presence=new_link_presence)

    def induced_substream_by_links(self, links_list):
        """
        Return the sub stream induced by links in *links_list* and nodes implicated in these links.

        :param links_list: list of links ((id1,id2) or (label1,label2))
        :return: A Stream Graph
        """

        # Check if element of links_list are labels or nodes id:
        if type(links_list[0][0]) == str:
            label_to_node = {v: k for k, v in self.node_to_label.items()}
            links_list = [(label_to_node[l[0]], label_to_node[l[1]]) for l in links_list]

        new_nodes = []
        new_node_presence = []
        new_nodes_to_label = {}
        new_nodes_to_id = {}
        nodes_to_new_nodes = defaultdict(lambda: len(nodes_to_new_nodes))
        new_links = []
        new_link_presence = []
        for l, lp in zip(self.links, self.link_presence):
            u, v = l
            if (u, v) in links_list or (v, u) in links_list:
                if u not in nodes_to_new_nodes:
                    new_u = nodes_to_new_nodes[u]
                    new_nodes_to_label[new_u] = self.node_to_label[u]
                    if self.node_to_id:
                        new_nodes_to_id[new_u] = self.node_to_id[u]
                    new_nodes.append(new_u)
                    new_node_presence.append(self.node_presence[u])
                else:
                    new_u = nodes_to_new_nodes[u]
                if v not in nodes_to_new_nodes:
                    new_v = nodes_to_new_nodes[v]
                    new_nodes_to_label[new_v] = self.node_to_label[v]
                    if self.node_to_id:
                        new_nodes_to_id[new_v] = self.node_to_id[v]
                    new_nodes.append(new_v)
                    new_node_presence.append(self.node_presence[v])
                else:
                    new_v = nodes_to_new_nodes[v]
                new_links.append((new_u, new_v))
                new_link_presence.append(lp)

        return stream_graph(times=self.times,
                            nodes=new_nodes,
                            node_to_label=new_nodes_to_label,
                            node_to_id=new_nodes_to_id,
                            node_presence=new_node_presence,
                            links=new_links,
                            link_presence=new_link_presence)

    def induced_substream_by_time_window(self, time_window):
        """
        Return the sub stream induced by T = [a,b].

        :param time_window: The time_window that delimit the substream
        :return: A Stream Graph
        """
        a, b = time_window
        new_nodes = []
        new_node_presence = []
        new_nodes_to_label = {}
        new_nodes_to_id = {}
        nodes_to_new_nodes = defaultdict(lambda: len(nodes_to_new_nodes))

        new_links = []
        new_link_presence = []

        for n, np in zip(self.nodes, self.node_presence):
            new_np = []
            for t0, t1 in zip(np[::2], np[1::2]):
                if t0 <= b and a <= t1:  # Intersection
                    new_np += [max(a, t0), min(b, t1)]
            if new_np:
                new_n = nodes_to_new_nodes[n]
                if self.node_to_label:
                    new_nodes_to_label[new_n] = self.node_to_label[n]
                if self.node_to_id:
                    new_nodes_to_id[new_n] = self.node_to_id[n]
                new_nodes.append(new_n)
                new_node_presence.append(new_np)

        for l, lp in zip(self.links, self.link_presence):
            new_lp = []
            for t0, t1 in zip(lp[::2], lp[1::2]):
                if t0 <= b and a <= t1:
                    new_lp += [max(a, t0), min(b, t1)]
            if new_lp:
                u, v = l
                new_u = nodes_to_new_nodes[u]
                new_v = nodes_to_new_nodes[v]
                new_links.append((new_u, new_v))
                new_link_presence.append(new_lp)

        return stream_graph(times=time_window,
                            nodes=new_nodes,
                            node_to_label=new_nodes_to_label,
                            node_to_id=new_nodes_to_id,
                            node_presence=new_node_presence,
                            links=new_links,
                            link_presence=new_link_presence)

    def substream(self, cluster, substream_id=None, return_new_node_label=False):
        """
        Return the substream corresponding to the cluster:
        A Stream Graph containing the nodes of the cluster along with the segmented links that exist between them.

        :param cluster: A sequence of temporal nodes: [(t0,t1,u),...,(t0,t1,v)]
        :param substream_id: (optional parameter) assign an id to the substream (usefull when dealing with components)
        :return:
        """

        if type(cluster[0][2]) == str:
            label_to_node = {v: k for k, v in self.node_to_label.items()}
            cluster = [(t0, t1, label_to_node[n]) for t0, t1, n in cluster]

        new_nodes = []
        new_node_presence = []
        nodes_to_new_nodes = defaultdict(lambda: len(nodes_to_new_nodes))
        new_nodes_to_label = {}
        new_nodes_to_id = {}

        for e in cluster:
            t0, t1, n = e
            if n in nodes_to_new_nodes:
                new_n = nodes_to_new_nodes[n]

                new_node_presence[new_n] += [t0, t1]
            else:
                new_n = nodes_to_new_nodes[n]
                if self.node_to_label:
                    new_nodes_to_label[new_n] = self.node_to_label[n]
                if self.node_to_id:
                    new_nodes_to_id[new_n] = self.node_to_id[n]
                new_nodes.append(new_n)
                new_node_presence.append([t0, t1])

        new_links = []
        new_link_presence = []
        for l, lp in zip(self.links, self.link_presence):
            u, v = l
            if u in nodes_to_new_nodes and v in nodes_to_new_nodes:
                new_u = nodes_to_new_nodes[u]
                new_v = nodes_to_new_nodes[v]
                lu, lv = [], []
                for t0, t1 in zip(lp[::2], lp[1::2]):
                    tu, tv = None, None
                    for a, b in zip(new_node_presence[new_u][::2], new_node_presence[new_u][1::2]):
                        if a <= t1 and t0 <= b:  # Intersection
                            tu = (max(a, t0), min(t1, b))
                            break
                    for a, b in zip(new_node_presence[new_v][::2], new_node_presence[new_v][1::2]):
                        if a <= t1 and t0 <= b:  # Intersection
                            tv = (max(a, t0), min(t1, b))
                            break
                    if tu and tv:
                        lu.append(tu)
                        lv.append(tv)
                if lu and lv:
                    new_links.append((new_u, new_v))
                    uv_presence = []
                    for tu, tv in zip(lu, lv):
                        uv_presence += [max(tu[0], tv[0]), min(tv[1], tv[1])]
                    new_link_presence.append(uv_presence)

        S = stream_graph(id=substream_id,
                         times=self.times,
                         nodes=new_nodes,
                         node_to_label=new_nodes_to_label,
                         node_to_id=new_nodes_to_id,
                         node_presence=new_node_presence,
                         links=new_links,
                         link_presence=new_link_presence)
        if return_new_node_label:
            return S, nodes_to_new_nodes
        return S

    ############################################################################################
    #          FILTERS( dot not use without reading doc, intended for StraphViz)               #
    ############################################################################################

    def filter_by_time_window(self, a, b):
        """
        Return the sub stream induced by T = [a,b].
        :param a:
        :param b:
        :return:
        """
        new_node_presence = []
        new_link_presence = []
        for n, np in zip(self.nodes, self.node_presence):
            new_np = []
            for t0, t1 in zip(np[::2], np[1::2]):
                if t0 <= b and a <= t1:  # Intersection
                    new_np += [max(a, t0), min(b, t1)]
            # if new_np:
            new_node_presence.append(new_np)

        for l, lp in zip(self.links, self.link_presence):
            new_lp = []
            for t0, t1 in zip(lp[::2], lp[1::2]):
                if t0 <= b and a <= t1:
                    new_lp += [max(a, t0), min(b, t1)]
            # if new_lp:
            new_link_presence.append(new_lp)

        return stream_graph(times=[a, b], nodes=self.nodes,
                            node_to_label=self.node_to_label,
                            node_to_id=self.node_to_id,
                            node_presence=new_node_presence,
                            links=self.links,
                            link_presence=new_link_presence)

    def filter_by_links(self, links_list):
        """
        Return the sub stream induced by links in *links_list* and nodes implicated in these links.
        :param links_list: list of links to filter by
        :return:
        """
        new_nodes = []
        new_node_presence = []
        new_nodes_to_label = {}
        new_nodes_to_id = {}
        nodes_to_new_nodes = defaultdict(lambda: len(nodes_to_new_nodes))
        links_target = [(js['source'], js['target']) for js in links_list]
        new_links = []
        new_link_presence = []
        for l, lp in zip(self.links, self.link_presence):
            u, v = l
            if (u, v) in links_target or (v, u) in links_target:
                if u not in nodes_to_new_nodes:
                    new_u = nodes_to_new_nodes[u]
                    new_nodes_to_label[new_u] = self.node_to_label[u]
                    new_nodes_to_id[new_u] = self.node_to_id[u]
                    new_nodes.append(new_u)
                    new_node_presence.append(self.node_presence[u])
                else:
                    new_u = nodes_to_new_nodes[u]
                if v not in nodes_to_new_nodes:
                    new_v = nodes_to_new_nodes[v]
                    new_nodes_to_label[new_v] = self.node_to_label[v]
                    new_nodes_to_id[new_v] = self.node_to_id[v]
                    new_nodes.append(new_v)
                    new_node_presence.append(self.node_presence[v])
                else:
                    new_v = nodes_to_new_nodes[v]
                new_links.append((new_u, new_v))
                new_link_presence.append(lp)
        return stream_graph(times=self.times,
                            nodes=new_nodes,
                            node_to_label=new_nodes_to_label,
                            node_to_id=new_nodes_to_id,
                            node_presence=new_node_presence,
                            links=new_links,
                            link_presence=new_link_presence)

    def filter_by_links_id(self, links_list):
        """
        Return the sub stream induced by links in *links_list* and nodes implicated in these links.
        :param links_list: list of links to filter by
        :return:
        """
        new_nodes = []
        new_node_presence = []
        new_nodes_to_label = {}
        new_nodes_to_id = {}
        nodes_to_new_nodes = defaultdict(lambda: len(nodes_to_new_nodes))
        links_target = [(int(js['source']), int(js['target'])) for js in links_list]
        new_links = []
        new_link_presence = []
        for l, lp in zip(self.links, self.link_presence):
            u, v = l
            if (self.node_to_id[u], self.node_to_id[v]) in links_target \
                    or (self.node_to_id[v], self.node_to_id[u]) in links_target:
                if u not in nodes_to_new_nodes:
                    new_u = nodes_to_new_nodes[u]
                    new_nodes_to_label[new_u] = self.node_to_label[u]
                    new_nodes_to_id[new_u] = self.node_to_id[u]
                    new_nodes.append(new_u)
                    new_node_presence.append(self.node_presence[u])
                else:
                    new_u = nodes_to_new_nodes[u]
                if v not in nodes_to_new_nodes:
                    new_v = nodes_to_new_nodes[v]
                    new_nodes_to_label[new_v] = self.node_to_label[v]
                    new_nodes_to_id[new_v] = self.node_to_id[v]
                    new_nodes.append(new_v)
                    new_node_presence.append(self.node_presence[v])
                else:
                    new_v = nodes_to_new_nodes[v]
                new_links.append((new_u, new_v))
                new_link_presence.append(lp)
        return stream_graph(times=self.times,
                            nodes=new_nodes,
                            node_to_label=new_nodes_to_label,
                            node_to_id=new_nodes_to_id,
                            node_presence=new_node_presence,
                            links=new_links,
                            link_presence=new_link_presence)

    def filter_by_links_label(self, links_list):
        """
        Return the sub stream induced by links in *links_list* and nodes implicated in these links.
        :param links_list: list of links to filter by
        :return:
        """
        new_nodes = []
        new_node_presence = []
        new_nodes_to_label = {}
        new_nodes_to_id = {}
        nodes_to_new_nodes = defaultdict(lambda: len(nodes_to_new_nodes))
        new_links = []
        links_target = [(js['source'], js['target']) for js in links_list]
        new_link_presence = []
        for l, lp in zip(self.links, self.link_presence):
            u, v = l
            if (self.node_to_label[u], self.node_to_label[v]) in links_target or \
                    (self.node_to_label[v], self.node_to_label[u]) in links_target:
                if u not in nodes_to_new_nodes:
                    new_u = nodes_to_new_nodes[u]
                    new_nodes_to_label[new_u] = self.node_to_label[u]
                    new_nodes_to_id[new_u] = self.node_to_id[u]
                    new_nodes.append(new_u)
                    new_node_presence.append(self.node_presence[u])
                else:
                    new_u = nodes_to_new_nodes[u]
                if v not in nodes_to_new_nodes:
                    new_v = nodes_to_new_nodes[v]
                    new_nodes_to_label[new_v] = self.node_to_label[v]
                    new_nodes_to_id[new_v] = self.node_to_id[v]
                    new_nodes.append(new_v)
                    new_node_presence.append(self.node_presence[v])
                else:
                    new_v = nodes_to_new_nodes[v]
                new_links.append((new_u, new_v))
                new_link_presence.append(lp)
        return stream_graph(times=self.times,
                            nodes=new_nodes,
                            node_to_label=new_nodes_to_label,
                            node_to_id=new_nodes_to_id,
                            node_presence=new_node_presence,
                            links=new_links,
                            link_presence=new_link_presence)

    def filter_by_nodes(self, nodes_list):
        """
        Return the sub stream graph induced by nodes in *nodes_target*.
        If an extremity of a link is in *nodes_list* we add the other node.
        Different from induced stream graph

        :param nodes_list:
        :return:
        """
        new_nodes = []
        new_node_presence = []
        new_nodes_to_label = {}
        new_nodes_to_id = {}
        nodes_to_new_nodes = defaultdict(lambda: len(nodes_to_new_nodes))
        for n, np in zip(self.nodes, self.node_presence):
            if n in nodes_list:
                new_n = nodes_to_new_nodes[n]
                new_nodes_to_label[new_n] = self.node_to_label[n]
                if self.node_to_id:
                    new_nodes_to_id[new_n] = self.node_to_id[n]
                new_nodes.append(new_n)
                new_node_presence.append(np)

        new_links = []
        new_link_presence = []
        for l, lp in zip(self.links, self.link_presence):
            u, v = l
            if u in nodes_list or v in nodes_list:
                if u not in nodes_to_new_nodes:
                    new_u = nodes_to_new_nodes[u]
                    new_nodes_to_label[new_u] = self.node_to_label[u]
                    if self.node_to_id:
                        new_nodes_to_id[new_u] = self.node_to_id[u]
                else:
                    new_u = nodes_to_new_nodes[u]

                if v not in nodes_to_new_nodes:
                    new_v = nodes_to_new_nodes[v]
                    new_nodes_to_label[new_v] = self.node_to_label[v]
                    if self.node_to_id:
                        new_nodes_to_id[new_v] = self.node_to_id[v]

                else:
                    new_v = nodes_to_new_nodes[v]

                new_links.append((new_u, new_v))
                new_link_presence.append(lp)

        return stream_graph(times=self.times,
                            nodes=new_nodes,
                            node_to_label=new_nodes_to_label,
                            node_to_id=new_nodes_to_id,
                            node_presence=new_node_presence,
                            links=new_links,
                            link_presence=new_link_presence)

    def filter_by_node_id(self, nodes_target):
        """
        Return the sub stream graph induced by nodes in *nodes_target*.
        Different from induced stream graph (that include nodes in links containing a node in nodes_list)
        :param nodes_target:
        :return:
        """
        nodes_target = [int(n) for n in nodes_target]
        new_nodes = []
        new_node_presence = []
        new_nodes_to_label = {}
        new_nodes_to_id = {}
        nodes_to_new_nodes = defaultdict(lambda: len(nodes_to_new_nodes))
        for n, np in zip(self.nodes, self.node_presence):
            if self.node_to_id[n] in nodes_target:
                new_n = nodes_to_new_nodes[n]
                new_nodes_to_label[new_n] = self.node_to_label[n]
                new_nodes_to_id[new_n] = self.node_to_id[n]
                new_nodes.append(new_n)
                new_node_presence.append(np)

        new_links = []
        new_link_presence = []
        for l, lp in zip(self.links, self.link_presence):
            u, v = l
            if self.node_to_id[u] in nodes_target or self.node_to_id[v] in nodes_target:
                if u not in nodes_to_new_nodes:
                    new_u = nodes_to_new_nodes[u]
                    new_nodes_to_label[new_u] = self.node_to_label[u]
                    new_nodes_to_id[new_u] = self.node_to_id[u]
                else:
                    new_u = nodes_to_new_nodes[u]

                if v not in nodes_to_new_nodes:
                    new_v = nodes_to_new_nodes[v]
                    new_nodes_to_label[new_v] = self.node_to_label[v]
                    new_nodes_to_id[new_v] = self.node_to_id[v]

                else:
                    new_v = nodes_to_new_nodes[v]

                new_links.append((new_u, new_v))
                new_link_presence.append(lp)

        return stream_graph(times=self.times,
                            nodes=new_nodes,
                            node_to_label=new_nodes_to_label,
                            node_to_id=new_nodes_to_id,
                            node_presence=new_node_presence,
                            links=new_links,
                            link_presence=new_link_presence)

    def filter_by_node_label(self, nodes_target):
        """
        Return the sub stream graph induced by nodes in *nodes_target*.
        :param nodes_target:
        :return:
        """
        new_nodes = []
        new_node_presence = []
        new_nodes_to_label = {}
        new_nodes_to_id = {}
        nodes_to_new_nodes = defaultdict(lambda: len(nodes_to_new_nodes))
        for n, np in zip(self.nodes, self.node_presence):
            if self.node_to_label[n] in nodes_target:
                new_n = nodes_to_new_nodes[n]
                new_nodes_to_label[new_n] = self.node_to_label[n]
                new_nodes_to_id[new_n] = self.node_to_id[n]
                new_nodes.append(new_n)
                new_node_presence.append(np)

        new_links = []
        new_link_presence = []
        for l, lp in zip(self.links, self.link_presence):
            u, v = l
            if self.node_to_label[u] in nodes_target or self.node_to_label[v] in nodes_target:
                if u not in nodes_to_new_nodes:
                    new_u = nodes_to_new_nodes[u]
                    new_nodes_to_label[new_u] = self.node_to_label[u]
                    new_nodes_to_id[new_u] = self.node_to_id[u]
                    # new_nodes.append(new_u)
                    # new_node_presence.append([])
                else:
                    new_u = nodes_to_new_nodes[u]

                if v not in nodes_to_new_nodes:
                    new_v = nodes_to_new_nodes[v]
                    new_nodes_to_label[new_v] = self.node_to_label[v]
                    new_nodes_to_id[new_v] = self.node_to_id[v]
                    # new_nodes.append(new_v)
                    # new_node_presence.append([])
                else:
                    new_v = nodes_to_new_nodes[v]

                new_links.append((new_u, new_v))
                new_link_presence.append(lp)

        return stream_graph(times=self.times,
                            nodes=new_nodes,
                            node_to_label=new_nodes_to_label,
                            node_to_id=new_nodes_to_id,
                            node_presence=new_node_presence,
                            links=new_links,
                            link_presence=new_link_presence)

    ########################################################
    #           Weakly Connected Components : DFS          #
    ########################################################

    def DFS_iterative(self, v, Neighborhood):
        """
        Performs a Depth First search from the node *v*
        :param v: Source node for the DFS
        :param Neighborhood: Neighborhood of the Stream Graph : (t0,t1,u) -> [(t0,t1,v),....]
        :return: The weakly connected component of *v* along with visited nodes (other nodes in the wcc)
        """
        visited = set()
        Q = [v]
        component = []
        while len(Q) > 0:
            v = Q.pop()
            if v not in visited:
                visited.add(v)
                component.append(v)
                for w in Neighborhood[v]:
                    Q.append(w)
        return component, visited

    def weakly_connected_components(self, streaming=False, reformat=False, free_memory=False):
        """
        Compute the weakly connected components of the Stream Graph.
        Time Complexity in $O(M+N)$.
        :return: Weakly Connected Components of the Stream Graph (a cluster : [(t0,t1,u),...,(t0,t1,v)]
        """

        if streaming:
            return self.wcc_streaming(reformat=reformat, free_memory=free_memory)

        components = []
        Neighborhood = self.neighborhood_with_node_presence()
        if free_memory:
            self.nodes = []
            self.node_presence = []
            self.links = []
            self.link_presence = []

        unvisited = set(Neighborhood.keys())
        while len(unvisited) > 0:
            v = unvisited.pop()
            current_component, visited = self.DFS_iterative(v, Neighborhood)
            unvisited -= visited
            components.append(current_component)

        return components

    def number_of_weakly_connected_component(self):
        """
        Return the number of weakly connected components in the stream graph.
        :return:
        """
        return len(self.weakly_connected_components())

    def wcc_as_substreams(self):
        """
        Return the weakly connected components of the stream graph
        in the form of substreams (stream graph induced by the component/cluster)
        :return:
        """
        list_WCC = self.weakly_connected_components()
        node_to_id_wcc = {v: [] for v in self.nodes}
        list_substreams = []
        if len(list_WCC) > 1:
            for i, c in enumerate(list_WCC):
                substream, node_to_new_node = self.substream(c, substream_id=i, return_new_node_label=True)
                list_substreams.append(substream)
                for v in node_to_new_node:
                    node_to_id_wcc[v].append((i, node_to_new_node[v]))  # (id_wcc, new id)
        else:
            self.id = 0  # Only one WCC
            list_substreams = [self]
            for v in self.nodes:
                node_to_id_wcc[v].append((0, v))
        return list_substreams, node_to_id_wcc

    #############################################################
    #           Weakly Connected Components : Streaming         #
    #############################################################

    def find_node(self, u, dict_components):
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

    def link_components(self, u, v, dict_components, rank):
        # Choose the "biggest component" to append the node
        if rank[u] > rank[v]:
            dict_components[v] = u
        else:
            dict_components[u] = v
            if rank[u] == rank[v] and u != v:
                rank[v] += 1
        return dict_components, rank

    def wcc_streaming(self, reformat=False, free_memory=False):
        ''''
        Compute Weakly Connected Components in the lightest way.
        '''
        node_to_wcc = {}
        rank = {}
        node_to_segmented_node = {}
        E = self.ordered_arrivals(free_memory=free_memory)  # self.ordered_events(free_memory = free_memory)
        print("Ordered Arrivals DONE")
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

                root_u = self.find_node(u, node_to_wcc)
                root_v = self.find_node(v, node_to_wcc)
                self.link_components(root_u, root_v,
                                     node_to_wcc,
                                     rank)

        if reformat:
            def reformat_components(dict_components):
                k = 0
                dict_roots = {}
                for u in dict_components:
                    if dict_components[u] == u:
                        dict_roots[u] = k
                        k += 1
                components = [[] for i in range(k)]
                for u in dict_components:
                    p = u
                    while dict_components[p] != p:
                        p = dict_components[p]
                    components[dict_roots[p]].append(u)
                return components

            wcc = reformat_components(node_to_wcc)

            return wcc
        return node_to_wcc

    #################################################################
    #           Alternative, weakly connected with union find       #
    #################################################################

    def union_of_nodes(self, u, v, dict_components, rank):
        self.link_components(
            self.find_node(u, dict_components),
            self.find_node(v, dict_components),
            dict_components,
            rank
        )
        return

    def weakly_connected_components_union_find(self):
        dict_components = {(n, t0, t1): (n, t0, t1)
                           for n, np in zip(self.nodes, self.node_presence)
                           for t0, t1 in zip(np[::2], np[1::2])
                           }
        rank = {(n, t0, t1): 0
                for n, np in zip(self.nodes, self.node_presence)
                for t0, t1 in zip(np[::2], np[1::2])}
        # From now on, we only consider nodes as couples (n,In) In : interval of presence of n
        #  Moreover the value in each key will be considered as components
        for l, lp in zip(self.links, self.link_presence):
            for t0, t1 in zip(lp[::2], lp[1::2]):
                # Récupérer les temps de présence de l[0] et l[1] correspondants au lien (l[0],l[1],t0,t1) courant
                u, v = self.get_according_node_presence_from_link(l, t0, t1)
                self.union_of_nodes(u, v, dict_components, rank)

        def reformat_components(dict_components):
            k = 0
            dict_roots = {}
            for u in dict_components:
                if dict_components[u] == u:
                    dict_roots[u] = k
                    k += 1
            components = [[] for i in range(k)]
            for u in dict_components:
                p = u
                while dict_components[p] != p:
                    p = dict_components[p]
                components[dict_roots[p]].append(u)
            return components

        components = reformat_components(dict_components)
        return components

    ###########################################################
    #               Weakly Bounded connected components       #
    ###########################################################

    def get_interaction_times(self):
        """
        O(M)
        :return:
        """
        interaction_times = defaultdict(list)
        for l, lp in zip(self.links, self.link_presence):
            interaction_times[l[0]] += lp
            interaction_times[l[1]] += lp
        return interaction_times

    def get_interaction_segmented_nodes(self):
        interactions_times = defaultdict(list)
        for l, lp in zip(self.links, self.link_presence):
            for lt0, lt1 in zip(lp[::2], lp[1::2]):
                u, v = self.get_according_node_presence_from_link(l, lt0, lt1)
                interactions_times[u] += [lt0, lt1]
                interactions_times[v] += [lt0, lt1]
        return interactions_times

    def get_node_presence_according_to_node_partition_from_link(self, l, lt0, lt1, node_partition):
        n1, n2 = l[0], l[1]
        np1, np2 = None, None
        for n1t0, n1t1 in zip(node_partition[n1][::2], node_partition[n1][1::2]):
            if n1t0 <= lt0 and lt1 <= n1t1:
                np1 = (n1, n1t0, n1t1)
                break
        for n2t0, n2t1 in zip(node_partition[n2][::2], node_partition[n2][1::2]):
            if n2t0 <= lt0 and lt1 <= n2t1:
                np2 = (n2, n2t0, n2t1)
                break
        return np1, np2

    def get_neighborhood_from_node_partition(self, node_partition):
        N = {(n, t0, t1): [] for n in self.nodes
             for t0, t1 in zip(node_partition[n][::2], node_partition[n][1::2])}
        for l, lp in zip(self.links, self.link_presence):
            for t0, t1 in zip(lp[::2], lp[1::2]):
                np1, np2 = self.get_node_presence_according_to_node_partition_from_link(l, t0, t1, node_partition)
                if np1 and np2:
                    N[np1].append(np2)
                    N[np2].append(np1)
        return N

    def get_ordered_augmented_neighborhood_from_node_partition(self, node_partition):
        N = {(n, t0, t1): [] for n in self.nodes
             for t0, t1 in zip(node_partition[n][::2], node_partition[n][1::2])}
        for l, lp in zip(self.links, self.link_presence):
            for t0, t1 in zip(lp[::2], lp[1::2]):
                np1, np2 = self.get_node_presence_according_to_node_partition_from_link(l, t0, t1, node_partition)
                if np1 and np2:
                    N[np1].append((np2, (t0, t1)))
                    N[np2].append((np1, (t0, t1)))
        print("old neighborhood:")
        print(N)
        for k in N:
            N[k] = [n for _, n in sorted(zip([v[1][0] for v in N[k]], N[k]), reverse=True)]  # popopo
        print("ordered neighborhood:")
        print(N)
        return N

    def weakly_bounded_connected_components(self, Neighborhood=None, refactor=False):
        """
        At leat O(2M)
        :return:
        """
        if not Neighborhood:
            bounds = self.get_links_bounds_per_nodes()
            Neighborhood = self.get_neighborhood_from_node_partition(bounds)
        wb_components = []
        unvisited = set(Neighborhood.keys())
        # print("Neighborhood :\n",Neighborhood)
        while len(unvisited) > 0:
            v = unvisited.pop()
            current_component, visited = self.DFS_iterative(v, Neighborhood)
            unvisited -= visited
            wb_components.append(current_component)
        if refactor:
            L = []
            for wcc in wb_components:
                L.append([(i[1], i[2], i[0]) for i in wcc])
            return L
        else:
            return wb_components

    ##########################################################################
    #               Snapshots                                                #
    ##########################################################################

    def to_snapshots(self):
        G_list = []
        E = self.ordered_events()
        current_a_l = {}
        t_old = self.times[0]
        last_c = E[0][0]
        for ev in E:
            c, t = ev[0], ev[1]
            if t != t_old:
                G_list.append((t_old, t, {k: list(v) for k, v in current_a_l.items()}))
            if c == 1:
                (b, e, u, v) = ev[1:]
                current_a_l[u].add(v)
                current_a_l[v].add(u)
            elif c == -1:
                if t == t_old and last_c > 0:
                    G_list.append((t_old, t, {k: list(v) for k, v in current_a_l.items()}))
                (e, u, v) = ev[1:]
                current_a_l[u].remove(v)
                current_a_l[v].remove(u)
            elif c == 2:
                (b, e, u) = ev[1:]
                current_a_l[u] = set()
            elif c == -2:
                if t == t_old and last_c > 0:
                    G_list.append((t_old, t, {k: list(v) for k, v in current_a_l.items()}))
                (e, u) = ev[1:]
                current_a_l.pop(u)
            t_old = t
            last_c = c
        return G_list

    def nb_cc_from_adjacency_list(self, a_l):
        nb_components = 0
        unvisited = set(a_l.keys())
        while len(unvisited) > 0:
            v = unvisited.pop()
            current_component, visited = self.DFS_iterative(v, a_l)
            unvisited -= visited
            nb_components += 1
        return nb_components

    def number_of_connected_components_in_snapshots(self):
        cnt_cc = 0
        E = self.ordered_events()
        current_a_l = {}
        t_old = self.times[0]
        last_c = E[0][0]
        for ev in E:
            c, t = ev[0], ev[1]
            if t != t_old:
                cnt_cc += self.nb_cc_from_adjacency_list(current_a_l)
            if c == 1:
                (b, e, u, v) = ev[1:]
                current_a_l[u].add(v)
                current_a_l[v].add(u)
            elif c == -1:
                if t == t_old and last_c > 0:
                    cnt_cc += self.nb_cc_from_adjacency_list(current_a_l)
                (e, u, v) = ev[1:]
                current_a_l[u].remove(v)
                current_a_l[v].remove(u)
            elif c == 2:
                (b, e, u) = ev[1:]
                current_a_l[u] = set()
            elif c == -2:
                if t == t_old and last_c > 0:
                    cnt_cc += self.nb_cc_from_adjacency_list(current_a_l)
                (e, u) = ev[1:]
                current_a_l.pop(u)
            t_old = t
            last_c = c
        return cnt_cc

    ##########################################################################
    #               Stable Connected Components                              #
    ##########################################################################

    def stable_connected_components(self, format="object_with_links"):
        """

        :param format: "cluster" to obtain StCC under the form of [(t0,t1,u),...]
                        "object" to obtain StCC under the form of objects of class "connected_components"
        :return:
        """
        return comp.compute_stable_connected_components(self, format=format)

    def number_of_stable_connected_component(self, format="cluster"):
        """
        Return the number of strongly connected components in the stream graph.

        :return:
        """
        return len(self.stable_connected_components(format=format))

    ##########################################################################
    #               Strongly Connected Components                            #
    ##########################################################################

    def strongly_connected_components(self, format="object", streaming_output=None,
                                      method="Direct", isolated_nodes=True, free_memory=False):
        """
        Compute the strongly connected components of the Stream Graph

        :param format: "cluster" to obtain SCC under the form of [(t0,t1,u),...]
                        "object" to obtain SCC under the form of objects of class "strongly_connected_components"
        :return:
        """
        if method == "Direct":
            return comp.compute_strongly_connected_components(self, format=format,
                                                              streaming_output=streaming_output,
                                                              isolated_nodes=isolated_nodes,
                                                              free_memory=free_memory)
        if method == "UF":
            return comp.strongly_connected_components_UF(self, format=format)

        if method == "FD":
            return etf.strongly_connected_components_ETF(self, format=format,
                                                         streaming_output=streaming_output,
                                                         isolated_nodes=isolated_nodes,
                                                         free_memory=free_memory)

    def number_of_strongly_connected_component(self, format="cluster"):
        """
        Return the number of strongly connected components in the stream graph.

        :return:
        """
        return len(self.strongly_connected_components(format=format))

    ############################################################
    #                     Condensation DAG                     #
    ############################################################

    def condensation_dag(self, format="object_with_links", free_memory=False, isolated_nodes=True):
        """
        Compute the Condensation Graph of the Stream Graph.

        :return: Condensation Graph
        """

        _, dag = comp.compute_strongly_connected_components(self, condensation_dag=True, format=format,
                                                            free_memory=free_memory, isolated_nodes=isolated_nodes)
        dag.set_id_scc_to_scc()

        if format == "cluster":
            dag.times = (min([c[0][0] for c in dag.c_nodes]), max([c[0][1] for c in dag.c_nodes]))
            dag.compute_links_inplace()  # Add Links

        elif format == "object" or format == "object_with_links":
            dag.times = (min([c.times[0] for c in dag.c_nodes]), max([c.times[1] for c in dag.c_nodes]))
            dag.set_index_node_to_id_scc()
            dag.compute_links_inplace()  # Add Links
            if format == "object_with_links":
                dag.set_index_segmented_node_to_id_scc([(t0, t1, n)
                                                        for n, np in zip(self.nodes, self.node_presence)
                                                        for t0, t1 in zip(np[::2], np[1::2])])
        return dag

    def stable_dag(self, format="object_with_links", free_memory=False, isolated_nodes=True):
        """
        Compute the Stable Graph of the Stream Graph.
        # TODO : Comparer kcores parallel sur Satble Dag et Kcores parallel sur Condensation DAG

        :return: Stable Graph
        """

        _, dag = comp.compute_stable_connected_components(self, stable_dag=True, format=format,
                                                          free_memory=free_memory, isolated_nodes=isolated_nodes)
        dag.set_id_scc_to_scc()

        if format == "cluster":
            dag.times = (min([c[0][0] for c in dag.c_nodes]), max([c[0][1] for c in dag.c_nodes]))
            dag.compute_links_inplace()  # Add Links

        elif format == "object" or format == "object_with_links":
            dag.times = (min([c.times[0] for c in dag.c_nodes]), max([c.times[1] for c in dag.c_nodes]))
            dag.set_index_node_to_id_scc()
            dag.compute_links_inplace()  # Add Links
            if format == "object_with_links":
                dag.set_index_segmented_node_to_id_scc([(t0, t1, n)
                                                        for n, np in zip(self.nodes, self.node_presence)
                                                        for t0, t1 in zip(np[::2], np[1::2])])
        return dag

    # def old_stable_dag(self, condensation_dag=None, format="object_with_links", free_memory=False):
    #     """
    #     Not very efficient (i think)
    #     :return:
    #     """
    #     if condensation_dag is None:
    #         _, condensation_dag = comp.compute_strongly_connected_components(self, condensation_dag=True, format=format,
    #                                                                          free_memory=free_memory)
    #     stable_dag = condensation_dag.get_stable_dag()
    #     stable_dag.set_id_scc_to_scc()
    #
    #     if format == "cluster":
    #         stable_dag.times = (min([c[0][0] for c in stable_dag.c_nodes]), max([c[0][1] for c in stable_dag.c_nodes]))
    #         stable_dag.compute_links_inplace()  # Add Links
    #
    #     elif format == "object" or format == "object_with_links":
    #         stable_dag.times = (
    #         min([c.times[0] for c in stable_dag.c_nodes]), max([c.times[1] for c in stable_dag.c_nodes]))
    #         stable_dag.set_index_node_to_id_scc()
    #         stable_dag.compute_links_inplace()  # Add Links
    #
    #         if format == "object_with_links":
    #             stable_dag.set_index_segmented_node_to_id_scc([(t0, t1, n)
    #                                                            for n, np in zip(self.nodes, self.node_presence)
    #                                                            for t0, t1 in zip(np[::2], np[1::2])])
    #
    #     return stable_dag

    ############################################################
    #               Kcores                                     #
    ############################################################
    def k_core(self, k, condensation_dag=None, n_jobs=-1):
        """
        Compute the *k*-core cluster of the Stream Graph.

        :param k: *k*-core
        :param condensation_dag: condensation_dag: Condensation DAG (optional)
        :return:
        """
        if condensation_dag is None:
            _, condensation_dag = comp.compute_strongly_connected_components(self, condensation_dag=True)
        k_core = condensation_dag.k_core(k, n_jobs=n_jobs)
        return k_core

    def core_number(self, stable_dag=None, n_jobs=-1):
        """
        Compute the core number of each temporal node.

        :param condensation_dag: Condensation DAG (optional)
        :return:
        """
        if stable_dag is None:
            stable_dag = self.stable_dag()
        core_number = stable_dag.core_number(n_jobs=n_jobs)
        isolated_nodes = self.get_isolated_nodes()
        core_number[0] = isolated_nodes
        return core_number

    def average_core_size(self, n_jobs=-1):
        kcores = self.core_number(n_jobs=n_jobs)
        T = self.times[1] - self.times[0]
        average_core_size = sum([k * sum([n[1] - n[0] for n in v])
                                 for k, v in kcores.items()]) / (self.nb_nodes() * T)
        return average_core_size

    ############################################################
    #               Kcliques                                   #
    ############################################################

    def k_clique(self, k, stable_dag=None, n_jobs=-1):
        """
        Compute the *k*-cliques of the Stream Graph
        :param k:
        :param stable_dag:
        :return:
        """
        if stable_dag is None:
            stable_dag = self.stable_dag()
        k_clique = stable_dag.k_clique(k, n_jobs=n_jobs)
        return k_clique

    def all_cliques(self, stable_dag=None, n_jobs=-1):
        """
        return all the cliques of the Stream Graph
        :param k:
        :param stable_dag:
        :return:
        """
        if stable_dag is None:
            stable_dag = self.stable_dag()
        all_cliques = stable_dag.all_cliques(n_jobs=n_jobs)
        return all_cliques

    def max_clique_number(self, stable_dag=None, n_jobs=-1):
        """
        Return the size of the maximal clique for each node (in the form of temporal windows : a cluster)
        :return:
        """
        print("NOT IMPLEMENTED YET !!")
        if stable_dag is None:
            stable_dag = self.stable_dag()
        all_cliques = stable_dag.all_cliques(n_jobs=n_jobs)
        max_clique = {v: [] for v in self.nodes}
        for k, cluster in all_cliques.items():
            for t0, t1, n in cluster:
                max_clique[n] = [t0, t1, n]
        # TODO: To finish...
        return

    def average_clique_size(self, n_jobs=-1):
        kcliques = self.all_cliques(n_jobs=n_jobs)
        T = self.times[1] - self.times[0]
        average_clique_size = sum([k * sum([n[1] - n[0] for i in v for n in i])
                                   for k, v in kcliques.items()]) / (self.nb_nodes() * T)
        return average_clique_size

    ###########################################################
    #                   Degree Partition                      #
    ###########################################################

    def neighborhood_node_2_link_presence(self):
        N = {(t0, t1, n): []
             for n, np in zip(self.nodes, self.node_presence)
             for t0, t1 in zip(np[::2], np[1::2])}
        for l, lp in zip(self.links, self.link_presence):
            for t0, t1 in zip(lp[::2], lp[1::2]):
                np1, np2 = self.get_according_node_presence_from_link(l, t0, t1)
                N[np1] += [t0, t1]
                N[np2] += [t0, t1]
        return N

    def count_arrival_departure_degree(self, N):
        # We count links arriving vs links departing
        arrival = Counter(N[::2])
        departure = Counter(N[1::2])
        k = 0
        degree = defaultdict(list)
        N = list(set(N))
        N.sort()
        for t in N:
            # Some specificity if multiple things happen at the same time
            if t in arrival and t in departure:
                if departure[t] == arrival[t]:
                    continue
                # More departure than arrival
                elif departure[t] > arrival[t]:
                    for i in range(departure[t] - arrival[t]):
                        try:
                            degree[k].append(t)
                        except:
                            degree[k] = [t]
                        try:
                            degree[k - 1].append(t)
                            k -= 1
                        except:
                            k = 0
                else:
                    for i in range(arrival[t] - departure[t]):
                        k += 1
                        try:
                            degree[k].append(t)
                        except:
                            degree[k] = [t]
                        try:
                            degree[k - 1].append(t)
                        except:
                            pass
            elif t in arrival:
                for i in range(arrival[t]):
                    k += 1
                    try:
                        degree[k].append(t)
                    except:
                        degree[k] = [t]
                    try:
                        degree[k - 1].append(t)
                    except:
                        pass
            elif t in departure:
                for i in range(departure[t]):
                    try:
                        degree[k].append(t)
                    except:
                        degree[k] = [t]
                    try:
                        degree[k - 1].append(t)
                        k -= 1
                    except:
                        k = 0
        return degree

    def degrees_with_presence(self, N=None):
        if not N:
            N = self.neighborhood_node_2_link_presence()
        N_degree = {n: {} for n in N}
        for n in N:
            if not N[n]:
                warn("No link for :" + str(n) + " probably an isolated node")
                continue
            N_degree[n] = self.count_arrival_departure_degree(N[n])
        return N_degree

    # def degrees_partition_v1(self, N=None):
    #     if not N:
    #         N = self.neighborhood_node_2_link_presence()
    #     N_ensemble = defaultdict(list)
    #     for n in N:
    #         if not N[n]:
    #             warn("No link for :" + str(n) + " probably an isolated node")
    #             continue
    #         d = self.count_arrival_departure_degree(N[n])
    #         # print("d :",d)
    #         d_transform = {k: [(t0, t1, n[0]) for t0, t1 in zip(d[k][::2], d[k][1::2])] for k in d}
    #         # print("d transform :", d_transform)
    #         for k in d_transform:
    #             N_ensemble[k] += d_transform[k]
    #     return N_ensemble

    def degrees_partition(self, interaction_times=None):
        # Probably the best function to get degrees
        # note : some split node are not merged together
        # because a discontinuity (de mesure nulle) has occured
        # and we choose to keep it this way (never know)
        if not interaction_times:
            interaction_times = self.get_interaction_times()
        degrees = defaultdict(set)
        isolated_nodes = self.get_isolated_nodes()
        degrees[0] = isolated_nodes
        for n in interaction_times:
            t_n = interaction_times[n]
            if t_n:
                weights = [0] * len(t_n)
                l = len(t_n) // 2
                #  Warning !! Seems strange but this a trick to consider arrival before departure.
                weights[::2] = [-1] * l
                weights[1::2] = [1] * l
                weights = [-x for _, x in sorted(zip(t_n, weights))]
                t_n = sorted(t_n)
                cum_weights = list(itertools.accumulate(weights))
                for i in range(len(weights)):
                    if cum_weights[i] == 0:
                        continue
                    degrees[cum_weights[i]].add((t_n[i], t_n[i + 1], n))
        return degrees

    def nodes_partition_on_degrees(self, interaction_times=None):
        if not interaction_times:
            interaction_times = self.get_interaction_times()
        nodes_partition = {n: [] for n in self.nodes}
        for n in interaction_times:
            t_n = interaction_times[n]
            if t_n:
                weights = [0] * len(t_n)
                l = len(t_n) // 2
                weights[::2] = [1] * l
                weights[1::2] = [-1] * l
                weights = [x for _, x in sorted(zip(t_n, weights))]
                t_n = sorted(t_n)
                cum_weights = list(itertools.accumulate(weights))
                for i in range(len(weights)):
                    if cum_weights[i] == 0:
                        continue
                    if t_n[i] != t_n[i + 1]:
                        nodes_partition[n] += [t_n[i], t_n[i + 1]]
        return nodes_partition

    def nodes_degree(self, interaction_times=None):
        if not interaction_times:
            interaction_times = self.get_interaction_times()
        nodes_degree = {n: [] for n in self.nodes}
        for n in interaction_times:
            t_n = interaction_times[n]
            if t_n:
                weights = [0] * len(t_n)
                l = len(t_n) // 2
                weights[::2] = [1] * l
                weights[1::2] = [-1] * l
                weights = [x for _, x in sorted(zip(t_n, weights))]
                t_n = sorted(t_n)
                cum_weights = list(itertools.accumulate(weights))
                for i in range(len(weights)):
                    if cum_weights[i] == 0:
                        continue
                    if t_n[i] != t_n[i + 1]:
                        nodes_degree[n] += [(t_n[i], t_n[i + 1], cum_weights[i])]
        return nodes_degree

    ###################################
    #             PATHS               #
    ###################################

    def times_to_reach(self, source=None, destination=None, start_time=None, enumerate=False, E=None):
        """
        Compute temporal paths in the Stream Graph.
        If both *source* and *destination* are specified: return the time to reach of the foremost path
        If only *source* is specified: return time to reach all other nodes (single source)
        If *source* and *destination* aren't specified: return times to reach pairwise.

        :param source: Can be a source node or a temporal source node(t_0,t_1,u)
        :param destination: Can be a destination node or a temporal destination node(t_0,t_1,v)
        :param enumerate: if True we enumerate the number of path
        :return: Times to reach
        """
        if source is not None:
            return ap.FoP(self, source, destination, start_time, E=E)
        else:
            return ap.FoP_pw(self)

    def times_to_reach_and_lengths(self, source=None, destination=None, start_time=None, enumerate=False, E=None):
        """
        Compute temporal paths in the Stream Graph.
        If both *source* and *destination* are specified: return the ttr and distance of the shortest foremost path
        If only *source* is specified: return ttr and distances to all other nodes (single source)
        If *source* and *destination* aren't specified: return ttr and distances pairwise.

        :param source: Can be a source node or a temporal source node(t_0,t_1,u)
        :param destination: Can be a destination node or a temporal destination node(t_0,t_1,v)
        :param enumerate: if True we enumerate the number of path
        :return: Times to reach and distances
        """
        if source is not None:
            return ap.SFoP(self, source, destination, start_time, E=E)
        else:
            return ap.SFoP_pw(self)

    def latencies(self, source=None, destination=None, enumerate=False, E=None):
        """
        Compute temporal paths in the Stream Graph.
        If both *source* and *destination* are specified: return the latency of the fastest path
        If only *source* is specified: return latencies to all other nodes (single source)
        If *source* and *destination* aren't specified: return latencies pairwise.

        :param source: Can be a source node or a temporal source node(t_0,t_1,u)
        :param destination: Can be a destination node or a temporal destination node(t_0,t_1,v)
        :param enumerate: if True we enumerate the number of path
        :return: Latencies
        """
        if source is not None:
            return ap.FP(self, source, destination, E=E)
        else:
            return ap.FP_pw(self)

    def latencies_and_lengths(self, source=None, destination=None, enumerate=False, E=None):
        """
        Compute temporal paths in the Stream Graph.
        If both *source* and *destination* are specified: return the latency_and_length of the shortest fastest path
        If only *source* is specified: return latencies and lengths to all other nodes (single source)
        If *source* and *destination* aren't specified: return latencies and lengths pairwise.

        :param source: Can be a source node or a temporal source node(t_0,t_1,u)
        :param destination: Can be a destination node or a temporal destination node(t_0,t_1,v)
        :param enumerate: if True we enumerate the number of path
        :return: Latencies and lengths
        """
        if source is not None:
            return ap. \
                SFP(self, source, destination, E=E)
        else:
            return ap.SFP_pw(self)

    def distances(self, source=None, destination=None, enumerate=False, E=None):
        """
        Compute temporal paths in the Stream Graph.
        If both *source* and *destination* are specified: return the distance of the shortest path
        If only *source* is specified: return distances to all other nodes (single source)
        If *source* and *destination* aren't specified: return distances pairwise.

        :param source: Can be a source node or a temporal source node(t_0,t_1,u)
        :param destination: Can be a destination node or a temporal destination node(t_0,t_1,v)
        :param enumerate: if True we enumerate the number of path
        :return: Distances
        """
        if source is not None:
            return ap.SP(self, source, destination, E=E)
        else:
            return ap.SP_pw(self)

    def distances_and_durations(self, source=None, destination=None, enumerate=False, E=None):
        """
        Compute temporal paths in the Stream Graph.
        If both *source* and *destination* are specified: return the distance and duration of the fastest shortest path
        If only *source* is specified: return distances and durations to all other nodes (single source)
        If *source* and *destination* aren't specified: return distances and durations pairwise.

        :param source: Can be a source node or a temporal source node(t_0,t_1,u)
        :param destination: Can be a destination node or a temporal destination node(t_0,t_1,v)
        :param enumerate: if True we enumerate the number of path
        :return: Distances and durations
        """
        if source is not None:
            return ap.FSP(self, source, destination, E=E)
        else:
            return ap.FSP_pw(self)

    ##############################
    #    Optimised Paths
    #############################

    def times_to_reach2(self, source=None, destination=None, start_time=None, enumerate=False):
        """
        Compute temporal paths in the Stream Graph.
        If both *source* and *destination* are specified: return the time to reach of the foremost path
        If only *source* is specified: return time to reach all other nodes (single source)
        If *source* and *destination* aren't specified: return times to reach pairwise.

        :param source: Can be a source node or a temporal source node(t_0,t_1,u)
        :param destination: Can be a destination node or a temporal destination node(t_0,t_1,v)
        :param enumerate: if True we enumerate the number of path
        :return: Times to reach
        """
        if source is not None:
            return ap2.FoP(self, source, destination, start_time)
        else:
            return ap2.FoP_pw(self)

    def times_to_reach_and_lengths2(self, source=None, destination=None, start_time=None, enumerate=False):
        """
        Compute temporal paths in the Stream Graph.
        If both *source* and *destination* are specified: return the ttr and distance of the shortest foremost path
        If only *source* is specified: return ttr and distances to all other nodes (single source)
        If *source* and *destination* aren't specified: return ttr and distances pairwise.

        :param source: Can be a source node or a temporal source node(t_0,t_1,u)
        :param destination: Can be a destination node or a temporal destination node(t_0,t_1,v)
        :param enumerate: if True we enumerate the number of path
        :return: Times to reach and distances
        """
        if source is not None:
            return ap2.SFoP(self, source, destination, start_time)
        else:
            return ap2.SFoP_pw(self)

    def latencies2(self, source=None, destination=None, enumerate=False):
        """
        Compute temporal paths in the Stream Graph.
        If both *source* and *destination* are specified: return the latency of the fastest path
        If only *source* is specified: return latencies to all other nodes (single source)
        If *source* and *destination* aren't specified: return latencies pairwise.

        :param source: Can be a source node or a temporal source node(t_0,t_1,u)
        :param destination: Can be a destination node or a temporal destination node(t_0,t_1,v)
        :param enumerate: if True we enumerate the number of path
        :return: Latencies
        """
        if source is not None:
            return ap2.FP(self, source, destination)
        else:
            return ap2.FP_pw(self)

    def latencies_and_lengths2(self, source=None, destination=None, enumerate=False):
        """
        Compute temporal paths in the Stream Graph.
        If both *source* and *destination* are specified: return the latency_and_length of the shortest fastest path
        If only *source* is specified: return latencies and lengths to all other nodes (single source)
        If *source* and *destination* aren't specified: return latencies and lengths pairwise.

        :param source: Can be a source node or a temporal source node(t_0,t_1,u)
        :param destination: Can be a destination node or a temporal destination node(t_0,t_1,v)
        :param enumerate: if True we enumerate the number of path
        :return: Latencies and lengths
        """
        if source is not None:
            return ap2.SFP(self, source, destination)
        else:
            return ap2.SFP_pw(self)

    def distances2(self, source=None, destination=None, enumerate=False):
        """
        Compute temporal paths in the Stream Graph.
        If both *source* and *destination* are specified: return the distance of the shortest path
        If only *source* is specified: return distances to all other nodes (single source)
        If *source* and *destination* aren't specified: return distances pairwise.

        :param source: Can be a source node or a temporal source node(t_0,t_1,u)
        :param destination: Can be a destination node or a temporal destination node(t_0,t_1,v)
        :param enumerate: if True we enumerate the number of path
        :return: Distances
        """
        if source is not None:
            return ap2.SP(self, source, destination)
        else:
            return ap2.SP_pw(self)

    def distances_and_durations2(self, source=None, destination=None, enumerate=False):
        """
        Compute temporal paths in the Stream Graph.
        If both *source* and *destination* are specified: return the distance and duration of the fastest shortest path
        If only *source* is specified: return distances and durations to all other nodes (single source)
        If *source* and *destination* aren't specified: return distances and durations pairwise.

        :param source: Can be a source node or a temporal source node(t_0,t_1,u)
        :param destination: Can be a destination node or a temporal destination node(t_0,t_1,v)
        :param enumerate: if True we enumerate the number of path
        :return: Distances and durations
        """
        if source is not None:
            return ap2.FSP(self, source, destination)
        else:
            return ap2.FSP_pw(self)

    #############################################################
    #               Properties
    #############################################################

    def links_properties(self, list_properties, to_pandas=False):
        # TODO : TO FINISH
        """

        :param list_properties:  properties of the Stream Graph
        :param to_pandas: If 'True' return a pandas DataFrame
        :return:
        """
        self.set_card_W()
        self.set_card_E()
        property_to_method = {'duration': self.duration,
                              'surface': self.card_E,
                              'stream_graph_degree': self.stream_graph_degree,
                              'expected_stream_graph_degree': self.expected_stream_graph_degree,
                              'average_degree': self.average_degree,
                              'nb_nodes': self.nb_nodes,
                              'nb_links': self.nb_links,
                              'node_duration': self.node_duration,
                              'link_duration': self.link_duration,
                              'coverage': self.coverage,
                              'density': self.density,
                              'average_clustering_coefficient': self.average_clustering,
                              'nb_of_wcc': self.number_of_weakly_connected_component,
                              'nb_of_scc': self.number_of_strongly_connected_component,
                              'average_core_size': self.average_core_size,
                              'average_clique_size': self.average_clique_size,
                              # 'uniformity':self.uniformity, Too costly !!
                              }
        default_properties = {'n': len(self.nodes),
                              'N': sum([len(item) / 2 for item in self.node_presence]),
                              'm': len(self.links),
                              'M': sum([len(item) / 2 for item in self.link_presence]),
                              }

        if list_properties == "all":
            dict_data = {p: [] for p in property_to_method}
            for p in property_to_method:
                chrono = time.time()
                dict_data[p].append(property_to_method[p]())
                print("Time to compute ", p, " : ", time.time() - chrono)
        elif list_properties != "default":
            dict_data = {p: [] for p in list_properties}
            for p in list_properties:
                dict_data[p].append(property_to_method[p]())

        for p in default_properties:
            dict_data[p] = [default_properties[p]]

        if to_pandas:
            if self.id:
                Index = pd.Index([self.id])
            else:
                Index = pd.Index(['Stream Graph'])
            D = pd.DataFrame(dict_data, index=Index)
            return D

        return dict_data

    def nodes_max_degree(self):
        max_degree = Counter()
        nodes_degree = self.nodes_degree()
        for n, d in nodes_degree.items():
            max_degree[n] = max(d)
        return max_degree

    def nodes_average_core_number(self, core_number=None):
        nodes_duration = self.nodes_duration()
        nodes_average_core_number = {n: 0 for n in self.nodes}
        if core_number is None:
            core_number = self.core_number()
        for c in core_number:
            for n in core_number[c]:
                t0, t1, u = n
                nodes_average_core_number[u] += (t1 - t0) * c
        for u in nodes_average_core_number:
            #  TODO : NODES DURATION
            nodes_average_core_number[u] = nodes_average_core_number[u] / nodes_duration[u]
        return

    def nodes_properties(self, list_properties, to_pandas=False):
        """
        # TODO : TO FINISH
        :param list_properties:  properties of the Stream Graph
        :param to_pandas: If 'True' return a pandas DataFrame
        :return:
        """
        self.set_card_W()
        self.set_card_E()
        property_to_method = {'degree': self.degrees,
                              'max_degree': self.nodes_max_degree,
                              'average_core_number': self.nodes_average_core_number,
                              'max_core_number': self.nodes_max_core_number,
                              'node_duration': self.nodes_duration,
                              'average_clique_number': self.nodes_average_clique_number,
                              'max_clique_number': self.nodes_max_clique_number,
                              'clustering_coefficient': self.clustering_coefficient,
                              'mean_degree_neighborhood': self.mean_degree_neighborhood,
                              'expected_node_degree': self.expected_node_degree,
                              'average_size_wcc': self.nodes_average_size_wcc,
                              'average_size_scc': self.nodes_average_size_scc,
                              }

        default_properties = {'n': len(self.nodes),
                              'N': sum([len(item) / 2 for item in self.node_presence]),
                              'm': len(self.links),
                              'M': sum([len(item) / 2 for item in self.link_presence]),
                              }

        if list_properties == "all":
            dict_data = {p: [] for p in property_to_method}
            for p in property_to_method:
                chrono = time.time()
                dict_data[p].append(property_to_method[p]())
                print("Time to compute ", p, " : ", time.time() - chrono)
        elif list_properties != "default":
            dict_data = {p: [] for p in list_properties}
            for p in list_properties:
                dict_data[p].append(property_to_method[p]())

        for p in default_properties:
            dict_data[p] = [default_properties[p]]

        if to_pandas:
            if self.id:
                Index = pd.Index([self.id])
            else:
                Index = pd.Index(['Stream Graph'])
            D = pd.DataFrame(dict_data, index=Index)
            return D

        return dict_data

    def properties(self, list_properties="default", to_pandas=False):
        """

        :param list_properties:  properties of the Stream Graph
        :param to_pandas: If 'True' return a pandas DataFrame
        :return:
        """
        self.set_card_W()
        self.set_card_E()
        property_to_method = {'duration': self.duration,
                              'stream_graph_degree': self.stream_graph_degree,
                              'expected_stream_graph_degree': self.expected_stream_graph_degree,
                              'average_degree': self.average_degree,
                              'nb_nodes': self.nb_nodes,
                              'nb_links': self.nb_links,
                              'node_duration': self.node_duration,
                              'link_duration': self.link_duration,
                              'coverage': self.coverage,
                              'density': self.density,
                              'average_clustering_coefficient': self.average_clustering,
                              'nb_of_wcc': self.number_of_weakly_connected_component,
                              'nb_of_scc': self.number_of_strongly_connected_component,
                              'average_core_size': self.average_core_size,
                              'average_clique_size': self.average_clique_size,
                              # 'uniformity':self.uniformity, Too costly !!
                              }
        default_properties = {'n': len(self.nodes),
                              'N': sum([len(item) / 2 for item in self.node_presence]),
                              'm': len(self.links),
                              'M': sum([len(item) / 2 for item in self.link_presence]),
                              }

        if list_properties == "all":
            dict_data = {p: [] for p in property_to_method}
            for p in property_to_method:
                chrono = time.time()
                dict_data[p].append(property_to_method[p]())
                print("Time to compute ", p, " : ", time.time() - chrono)
        elif list_properties != "default":
            dict_data = {p: [] for p in list_properties}
            for p in list_properties:
                dict_data[p].append(property_to_method[p]())

        for p in default_properties:
            dict_data[p] = [default_properties[p]]

        if to_pandas:
            if self.id:
                Index = pd.Index([self.id])
            else:
                Index = pd.Index(['Stream Graph'])
            D = pd.DataFrame(dict_data, index=Index)
            return D

        return dict_data

    ############################################################
    #               Extend graph                               #
    ############################################################

    def extend_on_node_presence(self):
        # Split on nodes...

        New_nodes = [(n, t0, t1)
                     for n in self.nodes
                     for t0, t1 in zip(self.node_presence[n][::2],
                                       self.node_presence[n][1::2])
                     ]

        dict_node_2_id = {n: i for n, i in zip(New_nodes, range(len(New_nodes)))}
        New_links = []
        pos_link = {}
        New_link_presence = []
        for l, lp in zip(self.links, self.link_presence):
            for lt0, lt1 in zip(lp[::2], lp[1::2]):
                np1, np2 = self.get_according_node_presence_from_link(l, lt0, lt1)
                n1, n2 = dict_node_2_id[np1], dict_node_2_id[np2]
                if (np1, np2) not in New_links and (np2, np1) not in New_links:
                    pos_link[(n1, n2)] = len(New_links)
                    pos_link[(n2, n1)] = len(New_links)
                    New_links.append((n1, n2))
                    New_link_presence.append([lt0, lt1])
                else:
                    New_link_presence[pos_link[(n1, n2)]] += [lt0, lt1]

        return stream_graph(nodes=[i for i in dict_node_2_id.values()],
                            node_presence=[[n[1], n[2]] for n in dict_node_2_id.keys()],
                            link_presence=New_link_presence,
                            links=New_links)

    def get_links_bounds_per_nodes(self, interaction_times=None):
        # Probably the best function to get bounds
        if not interaction_times:
            interaction_times = self.get_interaction_times()
        bounds = defaultdict(list)
        for n in interaction_times:
            t_n = interaction_times[n]
            if t_n:
                weights = [0] * len(t_n)
                l = len(t_n) // 2
                weights[::2] = [1] * l
                weights[1::2] = [-1] * l
                weights = [x for _, x in sorted(zip(t_n, weights))]
                t_n.sort()
                t0 = t_n[0]
                cnt = 0
                for i in range(len(weights)):
                    cnt_old = cnt
                    cnt += weights[i]
                    if cnt == 0:
                        bounds[n] += [t0, t_n[i]]
                    elif cnt_old == 0 and cnt == 1:
                        t0 = t_n[i]
        return bounds

    def extend_on_link_presence(self):
        # Split on links...
        new_nodes = []
        new_links = []
        new_node_presence = []
        new_link_presence = []
        interaction_times = {n: [] for n in self.nodes}
        for l, lp in zip(self.links, self.link_presence):
            interaction_times[l[0]] += lp
            interaction_times[l[1]] += lp
        nn = 0
        dict_nodes_2_new_nodes = {n: [] for n in self.nodes}
        for n in interaction_times:
            t_n = interaction_times[n]
            if not t_n:
                continue
            weights = np.empty(len(t_n))
            weights[::2] = 1
            weights[1::2] = -1
            ord = np.argsort(t_n)
            weights = weights[ord]
            weights = weights.cumsum()
            t_n.sort()
            t0 = t_n[0]

            for j in range(len(weights)):
                if weights[j] == 0:
                    dict_nodes_2_new_nodes[n].append(nn)
                    new_nodes.append(nn)
                    new_node_presence.append([t0, t_n[j]])
                    nn += 1
                elif weights[j - 1] == 0 and weights[j] == 1:
                    t0 = t_n[j]
        # for n in self.nodes:
        #     print(" Old node :",n," Time :",self.node_presence[n])
        #     for nn in dict_nodes_2_new_nodes[n]:
        #         print("\t New node :",nn," Time :",new_node_presence[nn])

        pos_link = {}
        for l, lp in zip(self.links, self.link_presence):
            # print("Old Link :",l,"Time :",lp)
            for lt0, lt1 in zip(lp[::2], lp[1::2]):
                for n1 in dict_nodes_2_new_nodes[l[0]]:
                    if new_node_presence[n1][0] <= lt0 and lt1 <= new_node_presence[n1][1]:
                        nl1 = n1
                        break
                for n2 in dict_nodes_2_new_nodes[l[1]]:
                    if new_node_presence[n2][0] <= lt0 and lt1 <= new_node_presence[n2][1]:
                        nl2 = n2
                        break
                if (nl1, nl2) not in new_links and (nl2, nl1) not in new_links:
                    pos_link[(nl1, nl2)] = len(new_links)
                    pos_link[(nl2, nl1)] = len(new_links)
                    new_links.append((nl1, nl2))
                    new_link_presence.append([lt0, lt1])
                    # print("\tNew link :",(nl1,nl2),"Time :",[lt0,lt1])
                else:

                    new_link_presence[pos_link[(nl1, nl2)]] += [lt0, lt1]
                    # print("\tNew link :",(nl1,nl2),"Time :",[lt0,lt1])
        print("After extension on links !")
        print("Nb nodes old :", len(self.nodes), "Nb new nodes :", len(new_nodes))
        print("Nb links old :", len(self.links), "Nb new links :", len(new_links))

        return stream_graph(nodes=new_nodes,
                            node_presence=new_node_presence,
                            links=new_links,
                            link_presence=new_link_presence)

    ###########################################################
    #                   Isolated Nodes                        #
    ###########################################################

    def get_isolated_nodes(self):
        """
        Get isolated temporal nodes (node with temporal degree == 0)
        :return:
        """
        E = self.augmented_ordered_links()
        isolated_nodes = []
        d = defaultdict(lambda: 0)
        last_activity = {}
        for l in E:
            c = l[0]
            if c == 1:
                _, t0, t1, u, v = l
                if u not in last_activity:
                    last_activity[u] = u[0]
                if v not in last_activity:
                    last_activity[v] = v[0]
                if d[u] == 0 and t0 != last_activity[u]:
                    isolated_nodes.append((last_activity[u], t0, u[2]))
                if d[v] == 0 and t0 != last_activity[v]:
                    isolated_nodes.append((last_activity[v], t0, v[2]))
                d[u] += 1
                d[v] += 1
            else:
                _, t1, u, v = l
                d[u] -= 1
                d[v] -= 1
                if d[u] == 0:
                    last_activity[u] = t1
                if d[v] == 0:
                    last_activity[v] = t1
        for n in last_activity:
            if n[1] != last_activity[n]:
                isolated_nodes.append((last_activity[n], n[1], n[2]))
        for n, np in zip(self.nodes, self.node_presence):
            for t0, t1 in zip(np[::2], np[1::2]):
                if (t0, t1, n) not in d:
                    isolated_nodes.append((t0, t1, n))
        return isolated_nodes

    ###########################################################
    #              Nodes partitionning, Neighborhood          #
    ###########################################################

    def split_nodes(self):
        """
        Must REMOVE ISOLATED NODE BEFORE (maybe)
        Divide nodes into minimum partition to characterize k_cores
        :return: Nodes partition
        """
        ### INIT :
        Neighborhood = [[] for n in self.nodes]
        interaction_times = [[] for n in self.nodes]  #
        for l, lp in zip(self.links, self.link_presence):
            n1 = l[0]
            n2 = l[1]
            Neighborhood[n1].append(n2)
            Neighborhood[n2].append(n1)
            interaction_times[n1] += lp  #
            interaction_times[n2] += lp  #
        ## ITER :
        depth = 1
        bounds = self.get_links_bounds_per_nodes(interaction_times)  # Get bounds (~where link exists)
        split_nodes = [set(item) for item in interaction_times]
        del interaction_times
        time_while = time.time()
        undone = set(self.nodes)
        # Iteritavely add iterate neighborhood interaction times
        while undone:
            split_old = split_nodes
            split_nodes = [set()
                           if n in undone
                           else split_old[n]
                           for n in self.nodes]
            done = set()
            for n in undone:
                for w in Neighborhood[n]:
                    split_nodes[n] |= split_old[w]
                if split_nodes[n] == split_old[n]:
                    done.add(n)
            undone -= done
            depth += 1
        del Neighborhood
        # print("DEPTH SPLIT : ",depth)
        # print("time while :",time.time()-time_while)
        ## FINAL
        time_final_step = time.time()
        split_nodes = [sorted(split_nodes[n]) for n in self.nodes]
        split_nodes_final = {n: [] for n in self.nodes}
        for n in self.nodes:
            t_n = split_nodes[n]
            for b0, b1 in zip(bounds[n][::2], bounds[n][1::2]):
                for t0, t1 in zip(t_n[:len(t_n) - 1], t_n[1:len(t_n)]):  # range(len(t_n) - 1):
                    if b0 <= t0 and t1 <= b1:
                        split_nodes_final[n] += [t0, t1]
        # print("Time final step :",time.time()-time_final_step)
        return split_nodes_final

    # DEPRECATED : get_neighborhood_from_ordered_nodes
    def get_according_node_from_link_and_node_presence(self, l, lt0, lt1, dict_nodes):
        """
        :param l:
        :param lt0:
        :param lt1:
        :param dict_nodes:
        :return:
        """
        n1, n2 = l[0], l[1]
        np1 = [(n1, n1t0, n1t1)
               for n1t0, n1t1 in zip(dict_nodes[n1][::2], dict_nodes[n1][1::2])
               if n1t0 >= lt0 and n1t1 <= lt1]
        np2 = [(n2, n2t0, n2t1)
               for n2t0, n2t1 in zip(dict_nodes[n2][::2], dict_nodes[n2][1::2])
               if n2t0 >= lt0 and n2t1 <= lt1]
        return np1, np2

    # DEPRECATED cf :get_neighborhood_from_ordered_nodes
    def get_neighborhood_from_nodes(self, dict_nodes):
        """
        :param dict_nodes:
        :return:
        """
        N = {(n, t0, t1): []
             for n in dict_nodes
             for t0, t1 in zip(dict_nodes[n][::2], dict_nodes[n][1::2])}
        for l, lp in zip(self.links, self.link_presence):
            for t0, t1 in zip(lp[::2], lp[1::2]):
                np1, np2 = self.get_according_node_from_link_and_node_presence(l, t0, t1, dict_nodes)
                for n1 in np1:
                    for n2 in np2:
                        if n2[1] < n1[2] and n1[1] < n2[2]:
                            N[n1].append(n2)
                            N[n2].append(n1)
        return N

    def get_neighborhood_from_ordered_nodes(self, split_nodes):
        N = {(n, t0, t1): []
             for n in split_nodes
             for t0, t1 in zip(split_nodes[n][::2], split_nodes[n][1::2])}
        for l, lp in zip(self.links, self.link_presence):
            n1 = l[0]
            n2 = l[1]
            # We pick the ones with the less time instant to evaluate
            if len(split_nodes[n1]) > len(split_nodes[n2]):
                t_current = split_nodes[n2]
            else:
                t_current = split_nodes[n1]
            # For each link we update neighbors
            for lt0, lt1 in zip(lp[::2], lp[1::2]):
                for t0, t1 in zip(t_current[::2], t_current[1::2]):
                    if lt0 <= t0 and t1 <= lt1:
                        N[(n1, t0, t1)].append((n2, t0, t1))
                        N[(n2, t0, t1)].append((n1, t0, t1))
        return N

    # Version with two passes on edges, DEPRECATED
    def algo_kcores_v1(self, components=None):
        """
        :return:
        """
        if components:
            split_nodes = self.split_nodes_components(components)
        else:
            split_nodes = self.split_nodes()
        a_l_tmp = self.get_neighborhood_from_nodes(split_nodes)
        degrees = {}
        a_l = {}
        for n in a_l_tmp:
            neighbors = [w[0] for w in a_l_tmp[n]]
            d = len(set(neighbors))
            if d > 0:
                degrees[n] = d
                a_l[n] = a_l_tmp[n]
        m_d = max(degrees.values())
        # t = time.time()
        D = [set([l for l in degrees if degrees[l] == i]) for i in range(m_d + 1)]
        # print("market t:",time.time()-t)
        k = 1
        L = defaultdict(list)
        visited = set()
        visited_tmp = set()
        for j in range(len(degrees)):
            for d in range(m_d + 1):
                if D[d]:
                    i = d
                    break
            if i > k:
                L[k] = visited_tmp
                visited |= visited_tmp
                k = i
                visited_tmp = set()
            v = D[i].pop()
            visited_tmp.add(v)
            for w in a_l[v]:
                if w not in visited and w not in visited_tmp:
                    D[degrees[w]].discard(w)
                    degrees[w] -= 1
                    D[degrees[w]].add(w)
        L[k] += visited_tmp
        return L

    def k_cores_partitionning(self):
        """

        :return:
        """
        split_nodes = self.split_nodes()
        # print("Nb new nodes :",sum([len(s) for s in split_nodes.values()]))
        a_l_tmp = self.get_neighborhood_from_ordered_nodes(split_nodes)
        degrees = {}
        a_l = {}
        for n in a_l_tmp:
            d = len(set([w[0] for w in a_l_tmp[n]]))
            if d > 0:
                degrees[n] = d
                a_l[n] = set(a_l_tmp[n])
        cores = algo_kcores_batagelj(a_l, degrees)
        L = defaultdict(list)
        for v in cores:
            L[cores[v]].append(v)
        return L

    #############################################################################
    #           K cores from (strongly) connected components                    #
    #############################################################################

    def get_links_components(self, components, bounds):
        dict_node_2_components = {n: k
                                  for k, scc in zip(range(len(components)), components)
                                  for n in scc}
        l_links = [[] for k in range(len(components))]
        l_link_presence = [defaultdict(list) for k in range(len(components))]
        for l, lp in zip(self.links, self.link_presence):
            for lt0, lt1 in zip(lp[::2], lp[1::2]):
                np1, np2 = self.get_node_presence_according_to_node_partition_from_link(l, lt0, lt1, bounds)
                l_links[dict_node_2_components[np1]].append((np1, np2))
                l_link_presence[dict_node_2_components[np1]][(np1, np2)] += [lt0, lt1]
        return l_links, l_link_presence

    def neighborhood_and_degrees_from_links(self, links, link_presence, scc):
        times_scc = sorted(set([t
                                for l in link_presence
                                for t in link_presence[l]]))  # Sort + remove duplicate
        time_2_pos = {t: i for t, i in zip(times_scc, range(len(times_scc)))}
        a_l = {(n[0], t0, t1): set()
               for n in scc
               for t0, t1 in zip(times_scc[time_2_pos[n[1]]:time_2_pos[n[2]]],
                                 times_scc[(time_2_pos[n[1]] + 1):(time_2_pos[n[2]] + 1)])}
        # Garder en mémoire les liens de chaque composantes fortement connexes....
        for l in links:
            n1 = l[0][0]
            n2 = l[1][0]
            for lt0, lt1 in zip(link_presence[l][::2], link_presence[l][1::2]):
                for t0, t1 in zip(times_scc[time_2_pos[lt0]:time_2_pos[lt1]],
                                  times_scc[(time_2_pos[lt0] + 1):(time_2_pos[lt1] + 1)]):
                    nn1 = (n1, t0, t1)
                    nn2 = (n2, t0, t1)
                    a_l[nn1].add(nn2)
                    a_l[nn2].add(nn1)
        degrees = {n: len(set([w[0] for w in a_l[n]])) for n in a_l}
        return a_l, degrees

    # OLD VERSION : DEPRECATED
    def k_cores_old(self):
        bounds = self.get_links_bounds_per_nodes()
        Neighborhood = self.get_neighborhood_from_node_partition(bounds)
        w_components = self.weakly_bounded_connected_components(Neighborhood)
        l_links, l_link_presence = self.get_links_components(w_components, bounds)

        # Decomposition of nodes according to each connected components
        L = defaultdict(list)
        # n_split_nodes = 0
        split_nodes = []
        for wcc, links, link_presence in zip(w_components, l_links, l_link_presence):
            a_l, degrees = self.neighborhood_and_degrees_from_links(links, link_presence, wcc)
            cores = algo_kcores_batagelj(a_l, degrees)  # Algo Batelgej
            for v in cores:
                L[cores[v]].append(v)
            # n_split_nodes+=len(a_l)
            split_nodes += [[n] for n in a_l]
            del cores, a_l, degrees
        return L, split_nodes  # ,n_split_nodes

    ##############################################################
    #       Kcores v3    (parallelized)                          #
    ##############################################################

    def neighborhood_and_degrees_from_links_v3(self, t0, t1, links):
        """

        :param t0:
        :param t1:
        :param links:
        :return:
        """
        a_l = defaultdict(set)
        for l in links:
            u = (l[0], t0, t1)
            v = (l[1], t0, t1)
            a_l[u].add(v)
            a_l[v].add(u)
        # degrees1 = {n: len(set([w[0] for w in a_l[n]])) for n in a_l}
        degrees = {n: len(a_l[n]) for n in a_l}
        return a_l, degrees

    def k_cores_v3(self):
        bounds = self.get_links_bounds_per_nodes()
        Neighborhood = self.get_neighborhood_from_node_partition(bounds)
        w_components = self.weakly_bounded_connected_components(Neighborhood)
        l_links, l_link_presence = self.get_links_components(w_components, bounds)
        del w_components
        L = defaultdict(list)
        split_nodes = []
        for links, link_presence in zip(l_links, l_link_presence):
            inter_times = sorted(set([t
                                      for l in link_presence
                                      for t in link_presence[l]]))
            time_2_pos = {t: i for t, i in zip(inter_times, range(len(inter_times)))}
            # Construisons les graphes un à un en suivant les fenêtres.
            inter_links = [[] for k in range(len(inter_times) - 1)]
            for l in links:
                for lp0, lp1 in zip(link_presence[l][::2], link_presence[l][1::2]):
                    for i in range(time_2_pos[lp0], time_2_pos[lp1]):
                        inter_links[i].append((l[0][0], l[1][0]))
            for i in range(len(inter_times) - 1):
                t0, t1 = inter_times[i], inter_times[i + 1]
                current_links = inter_links[time_2_pos[t0]]
                a_l, degrees = self.neighborhood_and_degrees_from_links_v3(t0, t1, current_links)
                cores = algo_kcores_batagelj(a_l, degrees)  # Algo Batelgej
                for v in cores:
                    L[cores[v]].append(v)
                # n_split_nodes+=len(a_l)
                split_nodes += [[n] for n in a_l]
        return L, split_nodes

    def get_links_from_scc(self, components):
        # Get bounds for each nodes, and sort them
        number_of_component = len(components)
        dict_node_2_components = {n: set() for n in self.nodes}
        for k, scc in zip(range(number_of_component), components):
            for n in scc[0]:
                dict_node_2_components[n].add(k)
        l_links = [[] for k in range(number_of_component)]
        l_link_presence = [defaultdict(list) for k in range(number_of_component)]
        for l, lp in zip(self.links, self.link_presence):
            potential_components = dict_node_2_components[l[0]] & dict_node_2_components[l[1]]
            for lt0, lt1 in zip(lp[::2], lp[1::2]):
                for comp_number in potential_components:
                    t0, t1 = components[comp_number][1], components[comp_number][2]
                    if lt0 <= t1 and t0 <= lt1:
                        u, v = l[0], l[1]
                        l_links[comp_number].append((u, v))
                        l_link_presence[comp_number][(u, v)] += [max(lt0, t0), min(lt1, t1)]
        return l_links, l_link_presence

    def k_cores_from_scc(self, sc_components=None):
        if not sc_components:
            sc_components = self.propagation(refactor=False)
        l_links, l_link_presence = self.get_links_from_scc(sc_components)
        del sc_components
        L = defaultdict(list)
        split_nodes = []
        for links, link_presence in zip(l_links, l_link_presence):
            inter_times = sorted(set([t
                                      for l in link_presence
                                      for t in link_presence[l]]))
            time_2_pos = {t: i for t, i in zip(inter_times, range(len(inter_times)))}
            # Si il n'y a que deux noeuds dans la composantes : heuristiquement ils sont dans le 1core:
            if len(links) == 1:
                u, v = links[0][0], links[0][1]
                t0, t1 = link_presence[links[0]][0], link_presence[links[0]][1]
                L[1].append((u, t0, t1))
                L[1].append((v, t0, t1))
                split_nodes.append([(u, t0, t1)])
                split_nodes.append([(v, t0, t1)])
            else:
                # Construisons les graphes un à un en suivant les fenêtres.
                inter_links = [[] for k in range(len(inter_times) - 1)]
                for l in links:
                    for lp0, lp1 in zip(link_presence[l][::2], link_presence[l][1::2]):
                        for i in range(time_2_pos[lp0], time_2_pos[lp1]):
                            inter_links[i].append((l[0], l[1]))
                for i in range(len(inter_times) - 1):
                    t0, t1 = inter_times[i], inter_times[i + 1]
                    current_links = inter_links[time_2_pos[t0]]
                    a_l, degrees = self.neighborhood_and_degrees_from_links_v3(t0, t1, current_links)
                    cores = algo_kcores_batagelj(a_l, degrees)  # Algo Batelgej
                    for v in cores:
                        L[cores[v]].append(v)
                    # n_split_nodes+=len(a_l)
                    split_nodes += [[n] for n in a_l]
        return L, split_nodes

    def save_plot(self, path, format='pdf'):
        wc_components = self.weakly_connected_components()
        wbc_components = self.weakly_bounded_connected_components()
        SCC_stream = self.compute_scc_stream()
        SCC_stream.compute_links()
        Kcores = SCC_stream.get_kcores()
        SCC = SCC_stream.refactor()
        degrees_partition = self.degrees_partition()
        random.shuffle(SCC)

        self.plot(saving_path=path + "sg", format=format)
        self.plot(clusters=wc_components, title="SG with WC Components", saving_path=path + "sg_with_wcc",
                  format=format)  # links= False)
        self.plot(clusters=wbc_components, title="SG with WBC Components", saving_path=path + "sg_with_wbcc",
                  format=format)  # links= False)
        color_SCC = self.plot(clusters=SCC, title="SG with SC Components", saving_path=path + "sg_with_scc",
                              format=format)
        SCC_stream.plot(color_SCC, saving_path=path + "scc_stream", format=format)
        self.plot_dict_clusters(dict_clusters=degrees_partition, title="-Degree", saving_path=path + "sg_with_degrees",
                                format=format)
        self.plot_dict_clusters(Kcores, title="Core", saving_path=path + "sg_with_kcores", format=format)

    #############################################################################
    #                       END                                                 #
    #############################################################################


def dump_stream_graph(S, path):
    with open(path, 'wb') as output_file:
        pickle.dump(S, output_file, -1)


def load_stream_graph(path):
    with open(path, 'rb') as input_file:
        S = pickle.load(input_file)
    return S


def postprocess_kcliques_into_clusters(Kcliques):
    """
    We only keep the biggest 'K' for each segmented temporal nodes
    :param Kcliques:
    :return:
    """
    K = defaultdict(list)
    seen = set()
    for k, v in sorted(Kcliques.items(), reverse=True):
        # print("\n k : ", k)
        for clique in v:
            for e in clique:
                if e not in seen:
                    seen.add(e)
                    K[k].append(e)
                    # print("e : ", e)
    return K


def adjacency_list_to_json(a_l):
    js = {'nodes': [],
          'links': []}
    for u in a_l:
        js['nodes'].append({'id': u})
        for v in a_l[u]:
            js['links'].append({'source': u,
                                'target': v})
    return js


def write_adjacency_list_to_json(a_l, path_json):
    js = {'nodes': [],
          'links': []}
    for u in a_l:
        js['nodes'].append({'id': u})
        for v in a_l[u]:
            js['links'].append({'source': u,
                                'target': v})
    with open(path_json, 'w') as file_output:
        json.dump(js, file_output)


def get_index_node_to_id_wcc(list_WCC):
    """
    For source node paths
    :param list_WCC:
    :return:
    """
    index_node_to_wcc = defaultdict(set)
    for SG in list_WCC:
        for n in SG.nodes:
            index_node_to_wcc[n].add(SG.id)
    return index_node_to_wcc


def get_index_segmented_node_to_id_wcc(list_WCC):
    """
    for temporal source node paths
    :param list_WCC:
    :return:
    """
    index_segmented_nodes_to_wcc = {}
    for SG in list_WCC:
        for n, np in zip(SG.nodes, SG.node_presence):
            for t0, t1 in zip(np[::2], np[1::2]):
                index_segmented_nodes_to_wcc[(t0, t1, n)] = [SG.id]
    return index_segmented_nodes_to_wcc


def condensation_dag_from_wcc(list_WCC):
    """
    Get the corresponding condensation DAG from a list of Weakly Connected Components (Stream Graphs Objects)
    :param list_WCC:
    :return:
    """
    n_wcc = len(list_WCC)

    def para_cdag(SG):
        return SG.condensation_dag()

    if n_wcc >= 4:
        # We parallelize only if the number of WCC is high enough :)
        list_CDAG = Parallel(n_jobs=-1, mmap_mode='r+')(
            delayed(para_cdag)(SG) for SG in list_WCC)
    else:
        list_CDAG = [SG.condensation_dag() for SG in list_WCC]
    return list_CDAG


