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

import copy
import math
import matplotlib.patches as mpatch
import matplotlib.pyplot as plt
import random
from collections import defaultdict, deque

from straph import components as comp
from straph import paths as pt
from straph import stream as sg
from straph.dags.dag import Dag
from straph.dags.stable_dag import StableDag


def compute_all_foremost_paths(dict_id_wcc_to_dag, index_node_to_scc, source, destination,
                               duration_threshold=None, start_comp=None):
    """
    Compute a foremost path between 'source' and 'destination" starting at 'begin time'.

    :param dict_id_wcc_to_dag:
    :param index_node_to_scc:
    :param source:
    :param destination:
    :param duration_threshold:
    :param start_comp:
    :return:
    """
    start = source[0]
    source = source[1]
    # print("Foremost path from ", source, " to ", destination, " starting at time ", start)

    # Get starting component
    if start_comp:
        id_wcc, start_comp = start_comp
    else:
        id_wcc, start_comp = None, None
        for (i, id_scc) in index_node_to_scc[source]:
            id_wcc = i
            c = dict_id_wcc_to_dag[id_wcc].id_comp_to_comp[id_scc]
            if c.times[0] <= start <= c.times[1]:
                if destination in c.nodes:
                    return [[id_scc]], 0, id_wcc
                start_comp = id_scc
                break

    if start_comp is None:
        print("Node : " + str(source) + " does not exist in the Stream Graph at time " + str(start) + " !")
    adj_list = defaultdict(list)
    end_time_comp = []
    G = dict_id_wcc_to_dag[id_wcc]
    if G.c_links:
        for l in G.c_links:
            # ONLY IF the destination is accesible from begin time
            c = G.id_comp_to_comp[l[1]]
            if c.times[0] >= start:
                adj_list[l[0]].append(l[1])
                end_time_comp.append(c.times[1])

    if duration_threshold is None:
        if not end_time_comp:
            return None, None, None
        duration_threshold = max(end_time_comp) + 1 - start

    # print("duration threshold :",duration_threshold)
    # print(" adj_list : ", adj_list)

    # Custom BFS on DAG
    def bfs_scc(a_l, start_cmp, dst):
        path_queue = [(start_cmp, [start_cmp])]
        foremost_duration = duration_threshold  # This variable, once assigned, is used as a threshold (yeah)
        while path_queue:
            (v, path) = path_queue.pop(0)
            # print(" v : ",v)
            if len(path_queue) > 0 and len(path_queue) % 10000 == 0:
                print(" len path queue :", len(path_queue))
            if v in a_l and a_l[v]:
                for cc in a_l[v]:
                    cmp = dict_id_wcc_to_dag[id_wcc].id_comp_to_comp[cc]
                    print(" comp id :", cc)
                    print(" comp times :", cmp.times)
                    if dst in cmp.nodes:
                        # print("PATH Found :",path)
                        foremost_duration = min(foremost_duration, cmp.times[0] - start)
                        yield path + [cc]
                        # return path + [c]
                    elif cmp.times[0] <= start + foremost_duration:
                        path_queue.append((cc, path + [cc]))

    # print(" Start BFS")
    foremost_paths = list(bfs_scc(adj_list, start_comp, destination))

    if not foremost_paths:
        return None, None, None

    path_times = []
    for p in foremost_paths:
        t = []
        for c in p:
            t.append(dict_id_wcc_to_dag[id_wcc].id_comp_to_comp[c].times[0])
        path_times.append(t)
    # print("Path times :", path_times)
    min_t = min([p[-1] for p in path_times])
    # print("Foremost time : ", min_t)
    fm_paths = []
    for p, t in zip(foremost_paths, path_times):
        if t[-1] == min_t:
            fm_paths.append(p)
    return fm_paths, min_t - start, id_wcc


class CondensationDag(Dag):
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
        A basic constructor for a ``CondensationDag``

        :param id:
        :param times:
        :param c_nodes : A list of nodes id (each node represents a ``StronglyConnectedComponent`` object:
        a set of nodes, a begin time, an end time)
        :param c_links : A list of directed link (each directed links connects two adjacent
        ``StronglyConnectedComponent``)
        :param id_comp_to_comp:
        :param node_to_id_comp:
        :param segmented_node_to_id_comp:
        :param adj_list:
        """
        super().__init__(id, times, c_nodes, c_links, id_comp_to_comp, node_to_id_comp,
                         segmented_node_to_id_comp, adj_list)

    def get_stable_dag(self):
        # Add stables parts as stable connected components
        stable_DAG = StableDag()
        stable_DAG.set_id(self.id)
        stable_DAG.times = self.times
        cnt_c_nodes = 0
        for cmp in self.c_nodes:
            stable_comps = cmp.get_stable_components(format="object")
            for c in stable_comps:
                c.id = cnt_c_nodes
                cnt_c_nodes += 1
            stable_DAG.add_nodes(stable_comps)
        return stable_DAG

    ################################
    #       FORMAT                 #
    ################################

    def cluster_to_object(self):
        new_cnodes = []
        for id_cc, cc in self.id_comp_to_comp.items():
            assert type(cc) == list
            new_cnodes = comp.StronglyConnectedComponent(id=id_cc, times=(cc[0][0], cc[0][1]),
                                                         nodes=set([c[2] for c in cc]))
        self.c_nodes = new_cnodes
        self.id_comp_to_comp = {cc.id: cc for cc in new_cnodes}

    ###############################
    #       Paths Methods         #
    ###############################

    def path_induced_substream(self, path, node_to_label=None,
                               path_bounds=None):
        """
        Transform a path in the condensation dag into a substream

        :param path: Sequence of identifiers of ``StronglyConnectedComponent`` objects in the current \
        ``CondensationDag``
        :param node_to_label:
        :param path_bounds:
        :return:
        """
        if type(path[0]) is int:
            path = [self.id_comp_to_comp[id_scc] for id_scc in path]

        new_nodes = {}
        nodes_to_new_nodes = defaultdict(lambda: len(nodes_to_new_nodes))
        new_node_to_label = {}
        new_node_to_id = {}
        new_links = {}
        t0_min, t1_max = math.inf, -math.inf
        for c in path:
            t0, t1 = c.times  # Initial Bounds
            if path_bounds is not None:
                t0, t1 = max(path_bounds[0], t0), min(path_bounds[1], t1)  # Bounds
            t0_min, t1_max = min(t0, t0_min), max(t1, t1_max)

            for n in c.nodes:
                if n in nodes_to_new_nodes:
                    new_n = nodes_to_new_nodes[n]
                    if t0 <= new_nodes[new_n][-1]:
                        new_nodes[new_n][-1] = t1
                    else:
                        new_nodes[new_n] += [t0, t1]
                else:
                    new_n = nodes_to_new_nodes[n]
                    new_node_to_id[new_n] = n
                    if node_to_label:
                        new_node_to_label[new_n] = node_to_label[n]
                    new_nodes[new_n] = [t0, t1]

            if c.links is not None:
                for l in c.links:
                    lt0, lt1 = l[0], l[1]
                    if lt1 < t0 or lt0 > t1:
                        # We do not consider links that end or begin outside the bounds
                        continue
                    lt0, lt1 = max(t0, lt0), min(t1, lt1)
                    u, v = l[2], l[3]
                    new_u, new_v = nodes_to_new_nodes[u], nodes_to_new_nodes[v]
                    if (new_u, new_v) in new_links:
                        if lt0 <= new_links[(new_u, new_v)][-1]:
                            new_links[(new_u, new_v)][-1] = lt1
                        else:
                            new_links[(new_u, new_v)] += [lt0, lt1]
                    else:
                        new_links[(new_u, new_v)] = [lt0, lt1]
                    # print("l :",(node_to_label[new_u],node_to_label[new_v]))

        F = sg.StreamGraph(times=[t0_min, t1_max],
                           nodes=list(new_nodes.keys()),
                           node_presence=list(new_nodes.values()),
                           node_to_label=new_node_to_label,
                           node_to_id=new_node_to_id,
                           links=list(new_links.keys()),
                           link_presence=list(new_links.values()))
        return F

    ####################################################
    #   1. Source-Destination Time to reach/ Latencies #
    ####################################################
    # TODO: Need to update functions in 1. with below functions in 2.

    def temporal_node_to_scc(self, node):
        """
        Return the ``StronglyConnectedComponent`` containing the temporal source *node*.

        :param node:
        :return:
        """
        n = node[2]
        t0, t1 = node[0], node[1]
        for c in self.c_nodes:
            if n in c.nodes and t0 <= c.times[0] <= c.times[1] <= t1:
                return c
        print("Node : " + str(n) + " does not exist in the Stream Graph at time " + str(t0) + " !")

    def node_to_scc(self, n):
        list_scc = []
        for c in self.c_nodes:
            if n in c.nodes:
                list_scc.append(c)
        return list_scc

    def time_to_reach(self, source, destination):

        if type(self.c_nodes[0]) == list:
            self.cluster_to_object()
        if type(source) is int:
            return self._time_to_reach(source, destination)
        else:
            return self._time_to_reach_temporal_nodes(source, destination)

    def _time_to_reach(self, source, destination):
        # TODO : to finish
        return

    def _time_to_reach_temporal_nodes(self, source, destination):
        """
        Return the time to reach the *destination* from the temporal source node *source* in the SG.

        :param source:
        :param destination:
        :return:
        """
        ttr = math.inf
        a_l = self.adjacency_list()
        start_comp = self.temporal_node_to_scc(source)
        if start_comp is None:
            return ttr
        st = source[0]
        print("Start comp:", start_comp)
        queue = deque([start_comp])  # comp
        visited = {start_comp.id}
        while queue:
            cmp = queue.popleft()
            if destination[2] in cmp.nodes and \
                    destination[0] <= cmp.times[0] <= cmp.times[1] <= destination[1]:
                ttr = max(min(ttr, cmp.times[0] - st), 0)
            if cmp.id in a_l and a_l[cmp.id]:
                for c_id in a_l[cmp.id]:
                    if c_id not in visited:
                        c = self.id_comp_to_comp[c_id]
                        if c.times[0] <= st + ttr:
                            queue.append(c)
                        visited.add(c_id)
        return ttr

    def latency(self, source, destination):

        if type(self.c_nodes[0]) == list:
            self.cluster_to_object()

        if type(source) is int:
            raise ValueError("Source node is not yet supported as input for latency computation in CondensationDag.")
        else:
            return self._latency_temporal_nodes(source, destination)

    def _latency_temporal_nodes(self, source, destination):
        """
        Return the latency between the temporal node *source* and the temporal node *destination* in the SG.

        :param source:
        :param destination:
        :return:
        """
        latency = math.inf
        unvisited = set()
        a_l = self.adjacency_list()
        for c in self.c_nodes:  # On itere sur les comp contenant source
            if source[2] in c.nodes and source[0] <= c.times[0] <= c.times[1] <= source[1]:
                unvisited.add(c.id)
                if destination[2] in c.nodes and destination[0] <= c.times[0] <= c.times[1] <= destination[1]:
                    latency = 0
                    return latency

        while len(unvisited) != 0:
            start_comp = self.id_comp_to_comp[unvisited.pop()]
            st = start_comp.times[1]  # starting time
            queue = deque([(start_comp, st)])  # comp, starting time
            visited = {start_comp.id}
            while queue:
                cmp, st = queue.popleft()
                # We can reset the starting time
                if source in cmp.nodes and source[0] <= cmp.times[0] <= cmp.times[1] <= source[1]:
                    st = cmp.times[1]
                    unvisited.discard(cmp.id)
                # Update latencies
                if destination[2] in cmp.nodes and destination[0] <= cmp.times[0] <= cmp.times[1] <= destination[1]:
                    latency = max(min(latency, cmp.times[0] - st), 0)

                if cmp.id in a_l and a_l[cmp.id]:
                    for c_id in a_l[cmp.id]:
                        if c_id not in visited:
                            c = self.id_comp_to_comp[c_id]
                            if c.times[0] <= st + latency:
                                queue.append((c, st))
                            visited.add(c.id)
        return latency

    #####################################################
    #   2. Single-Source Time to reach and Latencies    #
    #####################################################

    def times_to_reach_ss(self, source):

        if type(self.c_nodes[0]) == list:
            self.cluster_to_object()

        if type(source) is int:
            return self._times_to_reach_ss(source)
        else:
            return self._times_to_reach_temporal_nodes_ss(source)

    def postprocess_ttr(self, source, ttr_comp):
        ttr = {}
        for n in self.node_to_id_comp:
            potential_ttr = [ttr_comp[c] for c in self.node_to_id_comp[n] if c in ttr_comp]
            if potential_ttr:
                ttr[n] = min(potential_ttr)
        ttr[source] = 0
        return ttr

    def _times_to_reach_ss(self, source):
        id_start_comp = self.node_to_id_comp[source][0]  # The first SCC where the source appears
        start_comp = self.id_comp_to_comp[id_start_comp]
        ttr_comp = self._times_to_reach_comp_ss(id_start_comp, start_time=start_comp.times[0])
        ttr = self.postprocess_ttr(source, ttr_comp)
        return ttr

    def _times_to_reach_temporal_nodes_ss(self, source):
        potential_start_comp = self.segmented_node_to_id_comp[source]
        id_start_comp = None
        for id_comp in potential_start_comp:
            comp_t0, comp_t1 = self.id_comp_to_comp[id_comp].times
            if comp_t0 <= source[0] <= comp_t1:
                id_start_comp = id_comp
                break
        ttr_comp = self._times_to_reach_comp_ss(id_start_comp, start_time=source[0])
        ttr = self.postprocess_ttr(source[2], ttr_comp)
        return ttr

    def _times_to_reach_comp_ss(self, id_start_comp, start_time):
        """
        Return the times to reach from the temporal source node :*source* to every other node in the SG.

        :param id_start_comp:
        :param start_time:
        :return:
        """
        a_l = self.adjacency_list()
        times_to_reach_comp = {id_start_comp: 0}
        queue = deque([id_start_comp])
        visited = {id_start_comp}
        while queue:
            id_comp = queue.popleft()
            cmp = self.id_comp_to_comp[id_comp]
            times_to_reach_comp[id_comp] = cmp.times[0] - start_time
            if id_comp in a_l and a_l[id_comp]:
                for c_id in a_l[id_comp]:
                    if c_id not in visited:
                        queue.append(c_id)
                        visited.add(c_id)
        return times_to_reach_comp

    def latencies_ss(self, source):
        if type(self.c_nodes[0]) == list:
            self.cluster_to_object()

        if type(source) is int:
            return self._latencies_ss(source)
        else:
            return self._latencies_temporal_nodes_ss(source)

    def postprocess_latencies(self, latencies_comp, source):
        latencies = {}
        for n in self.node_to_id_comp:
            potential_latencies = [latencies_comp[c] for c in self.node_to_id_comp[n] if c in latencies_comp]
            if potential_latencies:
                latencies[n] = min(potential_latencies)
        latencies[source] = 0
        return latencies

    def _latencies_ss(self, source):
        ids_start_comp = self.node_to_id_comp[source]  # SCCs where the source appears
        latencies_comp = self._latencies_comp_ss(ids_start_comp)
        latencies = self.postprocess_latencies(latencies_comp, source)
        return latencies

    def _latencies_temporal_nodes_ss(self, source):
        ids_start_comp = self.segmented_node_to_id_comp[source]  # SCC where the source appears
        latencies_comp = self._latencies_comp_ss(ids_start_comp)
        latencies = self.postprocess_latencies(latencies_comp, source[2])
        return latencies

    def _latencies_comp_ss(self, ids_start_comp):
        """
        Return the latencies from the temporal source node: *source* to every other node in the SG.

        :param ids_start_comp:
        :return:
        """
        a_l = self.adjacency_list()
        latencies = {i: 0 for i in ids_start_comp}  # clefs : destination nodes ; valeurs : latency
        set_start_comps = set(ids_start_comp)
        unvisited = set(ids_start_comp)
        visited_to_st = {}
        # TODO : unvisited and set_start _comps necessaire ?
        while len(unvisited) != 0:
            id_start_comp = unvisited.pop()
            st = self.id_comp_to_comp[id_start_comp].times[1]  # starting time
            visited_to_st[id_start_comp] = st
            queue = deque([(id_start_comp, st)])  # comp, starting time
            while queue:
                id_comp, st = queue.popleft()
                cmp = self.id_comp_to_comp[id_comp]
                if id_comp in set_start_comps:
                    # We can reset the starting time
                    st = cmp.times[1]
                    unvisited.discard(id_comp)
                else:
                    if id_comp in latencies:
                        latencies[id_comp] = min(cmp.times[0] - st, latencies[id_comp])
                    else:
                        latencies[id_comp] = cmp.times[0] - st

                if id_comp in a_l and a_l[id_comp]:
                    for c_id in a_l[id_comp]:

                        if c_id not in visited_to_st or st > visited_to_st[c_id]:
                            # TODO : Verifier que la condition s'applique bien

                            # On ne doit pas avoir c_id plusieurs fois dans la queue !!
                            # Ou tester en sortie de pile !
                            # We leave later :)
                            queue.append((c_id, st))
                            visited_to_st[c_id] = st
        return latencies

    ##############################################
    #   3. Arbitrary Foremost and fastest Paths  #
    ##############################################

    def foremost_path(self, source, destination):
        """
        Compute a foremost path between 'source' and 'destination" starting at 'begin time'.

        :param source:
        :param destination:
        :return:
        """
        st = source[0]
        start_comp = self.temporal_node_to_scc(source)
        adj_list = self.adjacency_list()

        # Custom BFS on DAG
        def bfs_scc(a_l, start_cmp, dst):
            path_queue = deque([(start_cmp, [start_cmp.id])])
            ttr = math.inf  # This variable, once assigned, is used as a threshold (yeah)
            visited = {start_cmp.id}
            while path_queue:
                (cmp, path) = path_queue.popleft()

                if cmp.id in a_l and a_l[cmp.id]:
                    for v_id in a_l[cmp.id]:
                        if v_id not in visited:
                            v = self.id_comp_to_comp[v_id]
                            visited.add(v_id)

                            if dst[2] in v.nodes and dst[0] <= v.times[0] <= v.times[1] <= dst[1]:
                                ttr = max(min(ttr, v.times[0] - st), 0)
                                yield path + [v_id]
                            elif v.times[0] <= st + ttr:
                                path_queue.append((v, path + [v_id]))

        # print(" Start BFS")
        foremost_paths = list(bfs_scc(adj_list, start_comp, destination))
        if not foremost_paths:
            return None, None

        path_times = []
        for p in foremost_paths:
            t = []
            for c_id in p:
                c = self.id_comp_to_comp[c_id]
                t.append(c.times[0])
            path_times.append(t)
        # print("Path times :", path_times)
        min_t = min([p[-1] for p in path_times])
        # print("Foremost time : ", min_t)
        fm_paths = []
        for p, t in zip(foremost_paths, path_times):
            if t[-1] == min_t:
                fm_paths.append(p)
        return fm_paths, min_t - st

    def fastest_path(self, source, destination):
        """
        Compute the fastest path between 'source' and 'destination'

        :param source:
        :param destination:
        :return:
        """
        a_l = self.adjacency_list()
        ids_start_comp = set(self.segmented_node_to_id_comp[source])
        unvisited = copy.copy(ids_start_comp)  # SCC where the source appears
        visited_to_st = {}
        latency = math.inf
        fastest_paths = None
        while len(unvisited) != 0:
            id_start_comp = unvisited.pop()
            start_comp = self.id_comp_to_comp[id_start_comp]
            st = start_comp.times[1]  # starting time
            visited_to_st[id_start_comp] = st
            path_queue = deque([((id_start_comp, st), [id_start_comp])])  #  ((id comp, start time),path)
            while path_queue:
                e, path = path_queue.popleft()
                id_comp, st = e
                cmp = self.id_comp_to_comp[id_comp]

                if id_comp in ids_start_comp:
                    # We can reset the starting time
                    st = cmp.times[1]
                    path = [id_comp]
                    unvisited.discard(id_comp)
                if destination[2] in cmp.nodes and destination[0] <= cmp.times[0] <= cmp.times[1] <= destination[1]:
                    new_latency = max(cmp.times[0] - st, 0)
                    if new_latency < latency:
                        fastest_paths = path
                        latency = new_latency

                if id_comp in a_l and a_l[id_comp]:
                    for c_id in a_l[id_comp]:
                        if c_id not in visited_to_st:
                            visited_to_st[c_id] = st
                            path_queue.append(((c_id, st), path + [c_id]))
                        elif st > visited_to_st[c_id]:
                            # We leave later :)
                            path_queue.append(((c_id, st), path + [c_id]))
                            visited_to_st[c_id] = st
        return fastest_paths, latency

    ############################################
    #   4. All foremost path and fastest paths #
    ############################################

    def all_fastest_paths(self, source, destination):
        """
        Compute the fastest path between 'source' and 'destination'

        :param source:
        :param destination:
        :return:
        """
        a_l = self.adjacency_list()
        ids_start_comp = set(self.segmented_node_to_id_comp[source])
        unvisited = copy.copy(ids_start_comp)  # SCC where the source appears
        # We store visited component along with the starting time corresponding with the path that reached them
        visited_to_st = {}
        latency = math.inf
        fastest_paths = None
        while len(unvisited) != 0:
            id_start_comp = unvisited.pop()
            start_comp = self.id_comp_to_comp[id_start_comp]
            st = start_comp.times[1]  # starting time
            visited_to_st[id_start_comp] = st
            path_queue = deque([((id_start_comp, st), [id_start_comp])])  #  ((id comp, start time),path)
            while path_queue:
                e, path = path_queue.popleft()
                id_comp, st = e
                cmp = self.id_comp_to_comp[id_comp]

                if id_comp in ids_start_comp:
                    # We can reset the starting time
                    st = cmp.times[1]
                    path = [id_comp]
                    unvisited.discard(id_comp)
                if destination[2] in cmp.nodes and destination[0] <= cmp.times[0] <= cmp.times[1] <= destination[1]:
                    new_latency = max(cmp.times[0] - st, 0)
                    if new_latency < latency:
                        fastest_paths = [path]
                        latency = new_latency
                    elif new_latency == latency:
                        fastest_paths.append(path)
                if id_comp in a_l and a_l[id_comp]:
                    for c_id in a_l[id_comp]:
                        if c_id not in visited_to_st:
                            # We haven't seen the comp :)
                            visited_to_st[c_id] = st
                            path_queue.append(((c_id, st), path + [c_id]))
                        elif st > visited_to_st[c_id]:
                            # We leave later :)
                            path_queue.append(((c_id, st), path + [c_id]))
                            visited_to_st[c_id] = st
        return fastest_paths, latency

    def all_fastest_paths_ss(self, source):
        """
        Compute the fastest path between 'source' and all other nodes (in condensation DAG)

        :param source:
        :return:
        """
        a_l = self.adjacency_list()
        ids_start_comp = set(self.segmented_node_to_id_comp[source])
        unvisited = copy.copy(ids_start_comp)  # SCC where the source appears
        latencies_comp = {i: 0 for i in ids_start_comp}  # keys : destination nodes ;values : latency
        fastest_paths = {i: {i} for i in ids_start_comp}
        while len(unvisited) != 0:
            id_start_comp = unvisited.pop()
            start_comp = self.id_comp_to_comp[id_start_comp]
            st = start_comp.times[1]  # starting time
            path_queue = deque([((id_start_comp, st), [id_start_comp])])  #  ((id comp,start_time),path)
            while path_queue:
                e, path = path_queue.popleft()
                id_comp, st = e
                cmp = self.id_comp_to_comp[id_comp]
                new_latency = cmp.times[0] - st
                if id_comp in ids_start_comp:
                    # We can reset the starting time
                    st = cmp.times[1]
                    path = [id_comp]
                    unvisited.discard(id_comp)
                else:
                    if id_comp not in latencies_comp or latencies_comp[id_comp] > new_latency:
                        latencies_comp[id_comp] = new_latency
                        fastest_paths[id_comp] = {tuple(path)}
                    elif latencies_comp[id_comp] == new_latency:
                        fastest_paths[id_comp].add(tuple(path))
                    else:
                        continue

                if id_comp in a_l and a_l[id_comp]:
                    for c_id in a_l[id_comp]:
                        # We haven't seen the comp or its worth it to continue :)
                        cmp = self.id_comp_to_comp[c_id]
                        if c_id not in latencies_comp or latencies_comp[c_id] >= cmp.times[0] - st:
                            path_queue.append(((c_id, st), path + [c_id]))

        fastest_paths_nodes, latencies_nodes = self.postprocess_fastest_paths(fastest_paths, latencies_comp)
        del fastest_paths_nodes[source[2]]
        return fastest_paths_nodes, latencies_nodes

    def postprocess_fastest_paths(self, fastest_paths_comp, latencies_comp):
        # FP and Latencies in DAG to FP and Lat in Stream :
        latencies_nodes = defaultdict(lambda: math.inf)
        fastest_paths_nodes = {}
        for n in self.node_to_id_comp:
            lat = math.inf
            for c in self.node_to_id_comp[n]:
                if c in latencies_comp:
                    if latencies_comp[c] < lat:
                        fastest_paths_nodes[n] = fastest_paths_comp[c].copy()
                        latencies_nodes[n] = latencies_comp[c]
                        lat = latencies_comp[c]
                    elif latencies_comp[c] == lat:
                        fastest_paths_nodes[n] |= fastest_paths_comp[c]

        return fastest_paths_nodes, latencies_nodes

    ###############################################
    #   5. Shortest Fastest Path                  #
    #    (Hybrids Methods)                        #
    ###############################################

    def shortest_fastest_path_ss(self, source, node_to_label=None):
        """
        Compute Shortest Fastest Paths in a single source manner in a ``StreamGraph``
        using the current ``CondensationDag``.

        :param source:
        :param node_to_label:
        :return:
        """
        lengths = defaultdict(lambda: math.inf)
        fastest_paths_nodes, latencies_nodes = self.all_fastest_paths_ss(source)
        lengths[source[2]] = 0

        for destination in fastest_paths_nodes:
            list_fp = fastest_paths_nodes[destination]
            for fp in list_fp:
                if type(fp) == int:
                    fp = [fp]

                # print("\nfp :",fp)
                # print("destination :",destination)
                # assert source[2] in self.id_comp_to_comp[fp[0]].nodes
                # assert destination in self.id_comp_to_comp[fp[-1]].nodes
                # If length of the cdag path > 1 : last time of the first comp; first time of the last comp
                # Because they are SCC and all stream nodes are reachables at any instant.
                path_bounds = False
                if len(fp) > 1:
                    path_bounds = (self.id_comp_to_comp[fp[0]].times[1], self.id_comp_to_comp[fp[-1]].times[0])
                    # print("Path_bounds :",path_bounds)
                    # print("latencies destination :",latencies_nodes[destination])
                    # assert path_bounds[1] - path_bounds[0] == latencies_nodes[destination]

                # else:
                #     assert latencies_nodes[destination] == 0
                #     print("Path_bounds :",path_bounds)
                substream = self.path_induced_substream(fp, path_bounds=path_bounds,
                                                        node_to_label=node_to_label)  #  Get Substream for DAG paths
                stream_node_to_substream_node = {v: k for k, v in substream.node_to_id.items()}
                substream_source = stream_node_to_substream_node[source[2]]
                substream_destination = stream_node_to_substream_node[destination]
                # assert destination in stream_node_to_substream_node
                # assert source[2] in stream_node_to_substream_node
                #  Get source in Substream
                b, e = None, None
                npres = substream.node_presence[substream_source]
                for nt0, nt1 in zip(npres[::2], npres[1::2]):
                    if source[0] <= nt0 <= nt1 <= source[1]:
                        b, e = nt0, nt1
                        break
                # print("Source in stream :",source)
                # print("Source in substream :",(b,e,substream_source))
                # print("Destination in substream :",substream_destination)
                # substream.plot()
                # plt.show()
                #  Compute distances in Substream
                # print("Latency in Stream :",latencies_nodes[destination])
                # L,D = substream.latencies_and_lengths(((b, e, substream_source)))
                D = substream.distances((b, e, stream_node_to_substream_node[source[2]]))
                # if substream_destination not in D:
                #     print("nb wcc substream :",substream.number_of_weakly_connected_component())
                #     wcc = substream.weakly_connected_components()
                #     adj_list = substream.instant_graph(substream.times[0], label=False)
                #     plot_adjacency_list(substream, adj_list, label = False)
                #     substream.plot()
                #     plt.show()

                # if fp == [8, 25, 12, 13] and destination == 3 and source[2] == 7:
                #     # adj_list = substream.instant_graph(substream.times[0], label=False)
                #     # plot_adjacency_list(substream, adj_list, label = False)
                #     self.plot()
                #     substream.plot()
                #     plt.show()

                # print("Latency in substream :",L[substream_destination])
                # print("length in substream :",D[substream_destination])
                # assert L[substream_destination] == latencies_nodes[destination]
                lengths[destination] = min(lengths[destination], D[substream_destination])
            # print("Length SFP HYBRID :",lengths[destination])
        return latencies_nodes, lengths

    # END PATHS #

    #########################################
    #           Reachability Queries        #
    #########################################

    def is_reachable(self, source, target):
        if self.latency(source, target):
            return True
        else:
            return False

    #########################################
    #       Plot Functions                  #
    #########################################

    def plot_foremost_path(self, path, ttr, source, destination,
                           col="#8f246b"):
        """
        Draw a path on the current ``CondensationDag``

        :param path: A Stream Graph
        :param ttr:
        :param source:
        :param destination:
        :param col:
        :return:
        """

        if type(path[0]) is int:
            path = [self.id_comp_to_comp[id_scc] for id_scc in path]

        S = self.path_induced_substream(path)
        S.plot()
        rectangles = []
        for c in path:
            t0, t1 = c.times[0], c.times[1]
            for n in c.nodes:
                if n == source[2]:
                    color = '#4d79ff'
                elif n == destination[2]:
                    color = '#00cc00'
                else:
                    color = col
                rectangles.append(mpatch.Rectangle((t0, n - 0.15),
                                                   width=t1 - t0,
                                                   height=0.3,
                                                   edgecolor=color,
                                                   facecolor=color,
                                                   alpha=1
                                                   ))
            # Plot a single rectangle for COMP
            rectangles.append(mpatch.Rectangle((t0, min(c.nodes) + 0.15),
                                               width=t1 - t0,
                                               height=max(c.nodes) - min(c.nodes) - 0.3,
                                               edgecolor=col,
                                               facecolor=col,
                                               alpha=0.3
                                               ))

        ax = plt.gca()
        for r in rectangles:
            ax.add_artist(r)
        # Plot temporal source and temporal destination
        plt.plot([source[0]], [source[2]], color='#000099',
                 marker='^', alpha=1, markersize=13)
        plt.plot([source[0] + ttr], [destination[2]], color='#004d00',
                 marker='v', alpha=1, markersize=13)
        plt.title("Foremost Path from " + str(source) + " to " + str(destination))
        # plt.tight_layout()

    def plot_fastest_path(self, path, latency, source, destination,
                          col="#8f246b"):
        """
        Draw a path on the current ``CondensationDag``

        :param path:
        :param latency:
        :param source:
        :param destination:
        :param col:
        :return:
        """
        if type(path[0]) is int:
            path = [self.id_comp_to_comp[id_scc] for id_scc in path]
        S = self.path_induced_substream(path)
        S.plot()
        rectangles = []
        for c in path:
            t0, t1 = c.times[0], c.times[1]

            for n in c.nodes:
                if n == source[2]:
                    color = '#4d79ff'
                elif n == destination[2]:
                    color = '#00cc00'
                else:
                    color = col

                rectangles.append(mpatch.Rectangle((t0, n - 0.15),
                                                   width=t1 - t0,
                                                   height=0.3,
                                                   edgecolor=color,
                                                   facecolor=color,
                                                   alpha=1
                                                   ))
            # Plot a single rectangle for COMP
            rectangles.append(mpatch.Rectangle((t0, min(c.nodes) + 0.15),
                                               width=t1 - t0,
                                               height=max(c.nodes) - min(c.nodes) - 0.3,
                                               edgecolor=col,
                                               facecolor=col,
                                               alpha=0.3
                                               ))
        source_comp = path[0]

        ax = plt.gca()
        for r in rectangles:
            ax.add_artist(r)
        # Plot temporal source and temporal destination
        plt.plot([source_comp.times[1]], [source[2]], color='#000099',
                 marker='o', alpha=1, markersize=13)
        plt.plot([source_comp.times[1] + latency], [destination[2]], color='#004d00',
                 marker='d', alpha=1, markersize=13)
        plt.title("Fastest Path from " + str(source) + " to " + str(destination))
        # plt.tight_layout()

    def refactor_path(self, path_comp, path_times, source, destination):
        print("Path comp :", path_comp)
        print("Path Time :", path_times)
        P = pt.Path([], [])
        cur_node = source
        for i in range(len(path_comp) - 1):
            print("Current comp:", path_comp[i])
            print("Next comp:", path_comp[i + 1])
            cur_comp = self.c_nodes[path_comp[i]]
            next_comp = self.c_nodes[path_comp[i + 1]]
            intersec_nodes = cur_comp.nodes & next_comp.nodes
            if cur_node in intersec_nodes:
                # No need to jump
                continue
            else:
                next_node = random.choice(list(intersec_nodes))
            print("Intersec_nodes :", intersec_nodes)
            print("current node :", cur_node, " next node:", next_node)
            path_inside_comp = cur_comp.random_path(cur_node, next_node)
            print("path_inside_comp", path_inside_comp)
            cur_node = next_node
            for j in range(len(path_inside_comp[1]) - 1):
                P.add_link((path_inside_comp[1][j], path_inside_comp[1][j + 1]), path_times[i])
            print()
        path_inside_comp = self.c_nodes[path_comp[-1]].random_path(cur_node, destination)
        for j in range(len(path_inside_comp[1]) - 1):
            P.add_link((path_inside_comp[1][j], path_inside_comp[1][j + 1]), path_inside_comp[0])
        print("path times :", P.times)
        print("path nodes :", P.links)
        return P
