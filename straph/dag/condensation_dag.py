import msgpack, math, copy
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import pathlib, re, os, time
import matplotlib.collections as mcol
import matplotlib.patches as mpatch

import random

from joblib import Parallel, delayed
from collections import defaultdict, deque, Counter
from string import ascii_uppercase
from straph import paths as pt
from straph import stream as sg
from straph import components as comp
from straph.dag import stable_dag as sdag
import gc
from straph.utils import profile_shit,plot_adjacency_list

plt.rcParams["figure.figsize"] = (15, 10)


# TODO: Virer Kcores,Kcliques et compagnie (Garder les méthodes de chemins, sauf els foireuses genre SP).

def load_scc_dag(path_nodes, path_links=None):
    S = condensation_dag()
    with open(path_nodes, 'rb') as input:
        unpacker = msgpack.Unpacker(input)
        for i in unpacker:
            S.add_node(i)
    if path_links:
        with open(path_links, 'rb') as input:
            unpacker = msgpack.Unpacker(input)
            for i in unpacker:
                S.add_link(i)
    return S


def load_scc_dag_scf(path_scc_dag):
    '''
    Load a SCC DAG
    :param path_scc_dag:
    :return:
    '''
    dict_dag = {}
    with open(path_scc_dag, 'rb') as input:
        unpacker = msgpack.Unpacker(input, use_list=False)
        for i in unpacker:
            # if len(i[1]) > 1:  # We are only interested with DAG with at least 2 nodes
            S = condensation_dag()
            for n in i[1]:
                S.add_node(n)
            if len(i) == 3:
                for l in i[2]:
                    S.add_link(l)
            S.set_id(i[0])
            dict_dag[i[0]] = S
    return dict_dag


def compute_dict_offset(path_scc_scf, storage_path=False):
    scc_2_offset = {}
    with open(path_scc_scf, 'rb') as file_input:
        U = msgpack.Unpacker(file_input, use_list=False)
        offset = U.tell()
        for i in U:  # SCC components are ordered by their end time
            id_wcc, id_scc = i[0], i[1]
            scc_2_offset[(id_wcc, id_scc)] = offset
            offset = U.tell()
    if storage_path:
        with open(storage_path, 'wb') as output:
            msgpack.dump(scc_2_offset, output)
    return scc_2_offset


def load_dag_from_ids(id_wcc, id_scc, path_scc_scf, dict_offset=None):
    if not dict_offset:
        dict_offset = compute_dict_offset(path_scc_scf)

    with open(path_scc_scf, 'rb') as file_input:
        file_input.seek(dict_offset[(id_wcc, id_scc)])
        i = msgpack.Unpacker(file_input, use_list=False).unpack()
        print(i)
    return condensation_dag(id=i[0], c_nodes=i[1], c_links=i[2])


def read_global_dag(dag_path):
    G = condensation_dag()
    with open(dag_path, 'rb') as ipt:
        U = msgpack.Unpacker(ipt, use_list=False)
        for i in U:
            for n in i[1]:
                G.add_node(n)
            if len(i) > 2:
                for l in i[2]:
                    G.add_link(l)
    return G


def merge_scc_dags(l_scc):
    '''
    Merge every other streams into the first one, we assume there's no connection between them
    and that links aren't constructed yet
    :param l_scc
    :return:
    '''
    scc1 = l_scc[0]
    if scc1.times:
        min_t, max_t = scc1.times[0], scc1.times[1]
    for scc in l_scc[1:]:
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

    line:       ((a_x,_y),(b_x,b_y))
    ax:         matplotlib axes
    direction:  'left' or 'right'
    size:       size of the arrow in fontsize points
    color:      color of arrow (should be coherent with line color, or not )
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


class condensation_dag:
    def __init__(self,
                 id=None,
                 times=None,
                 c_nodes=None,
                 c_links=None,
                 id_scc_to_scc=None,
                 node_to_id_scc=None,
                 segmented_node_to_id_scc=None,
                 a_l=None,
                 ):
        '''
        A basic constructor for the condensation DAG
        :param c_nodes : A list of SCC nodes (each c node represent a SCC : a set of nodes, a begin time, an end time)
        :param c_links : A list of directed link (each link represent two connected SCC)
        '''
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
        self.node_to_id_scc = node_to_id_scc
        self.segmented_node_to_id_scc = segmented_node_to_id_scc
        self.id_scc_to_scc = id_scc_to_scc
        self.a_l = a_l

    def __repr__(self):
        rep = "\nId wcc :" + str(self.id)
        rep += "\nNodes :" + repr(self.c_nodes)
        rep += "\nLinks :" + str(self.c_links)
        return rep


    def describe(self):
        '''
        :return:
        '''
        print("DAG Id:",self.id)
        print("Nb of Nodes (comp) : ", len(self.c_nodes))
        print("Nb of stream nodes : ", len(set([n for c in self.c_nodes for n in c.nodes])))
        print("Nb of links (betw. comp): ", len(self.c_links))
        stream_links = set()
        degree_in = {n: 0 for n in self.id_scc_to_scc}
        degree_out = {n: 0 for n in self.id_scc_to_scc}
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



    def set_index_node_to_id_scc(self):
        node_to_id_scc = defaultdict(list)
        # order comp by their starting time
        for c in sorted(self.c_nodes, key=lambda x: x.times[0]):
            for n in c.nodes:
                node_to_id_scc[n].append(c.id)
        self.node_to_id_scc = node_to_id_scc

    def set_index_segmented_node_to_id_scc(self, segmented_nodes):
        if not self.node_to_id_scc:
            self.set_index_node_to_id_scc()
        segmented_node_to_id_scc = defaultdict(list)
        for n in segmented_nodes:
            for id_scc in self.node_to_id_scc[n[2]]:
                c = self.id_scc_to_scc[id_scc]
                if n[0] <= c.times[0] <= n[1]:
                    segmented_node_to_id_scc[n].append(id_scc)
        self.segmented_node_to_id_scc = segmented_node_to_id_scc

    def set_id_scc_to_scc(self):
        id_scc_to_scc = {}
        if type(self.c_nodes[0]) == list:
            cnt_id_scc = 0
            for c in self.c_nodes:
                id_scc_to_scc[cnt_id_scc] = c
                cnt_id_scc += 1
        else:
            for c in self.c_nodes:
                id_scc_to_scc[c.id] = c
        self.id_scc_to_scc = id_scc_to_scc

    def set_id(self, id):
        if id is None:
            self.id = 0
        else:
            self.id = id

    def set_index_id_scc_to_scc(self, index):
        self.id_scc_to_scc = index

    def add_node(self, n):
        self.c_nodes.append(n)

    def add_nodes(self, l):
        self.c_nodes += l

    def add_link(self, l):
        self.c_links.append(l)

    def store(self, output_file):
        packer = msgpack.Packer(use_bin_type=True)
        output_file.write(packer.pack((self.id, self.c_nodes, self.c_links)))

    def store_links(self, path, id):
        packer = msgpack.Packer(use_bin_type=True)
        if self.c_links:
            with open(path + "scc_dag_" + str(id) + "_links.mspk", 'wb') as output:
                output.write(packer.pack(self.c_links))

    def refactor(self):
        return [[(cc.times[0], cc.times[1], n) for n in cc.nodes] for cc in self.c_nodes]

    def plot_as_nx(self, label=True):
        fig = plt.figure()
        g_adjacency_list = {i: [] for i, _ in enumerate(self.c_nodes)}
        print("a_l : ", g_adjacency_list)
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

    def plot_custom(self, label=True, arrow=True, fontsize=30, path=[]):
        fig = plt.figure()
        ax = plt.gca()
        pos = defaultdict(list)
        min_x = math.inf
        max_x = -math.inf
        max_y = 0
        min_y = 0
        for i, n in enumerate(self.c_nodes):
            t0, t1 = n.times
            y = np.mean(list(n.nodes))
            x = (t0 + t1) / 2
            min_x = min(min_x, t0)
            max_x = max(max_x, t1)
            max_y = max(max_y, y)
            pos[i] = [x, y]

        g_adjacency_list = self.adjacency_list()
        # Plot nodes and links : plot nodes in increasing order
        # adjust 'y' depending on the number of neighbors
        visited = set()

        E = list()

        if path:
            E_path = list()
            # If we choose to plot a path on the DAG, we need a second edge_collections
            for n, _ in sorted(pos.items(), key=lambda x: x[1][0]):
                for u in g_adjacency_list[n]:
                    if (n, u) in path or (u, n) in path:
                        E_path.append((pos[n], pos[u]))
                    else:
                        E.append((pos[n], pos[u]))

        else:
            for n, _ in sorted(pos.items(), key=lambda x: x[1][0]):
                for u in g_adjacency_list[n]:
                    E.append((pos[n], pos[u]))
        xy = np.asarray([pos[v] for v in pos])
        ax.scatter(xy[:, 0], xy[:, 1],
                   c="#339966",
                   s=120,
                   marker='o',
                   alpha=0.7
                   )
        if label:
            for n in pos:
                ax.annotate(n, pos[n], fontsize=22)

        edge_collections = mcol.LineCollection(E, colors=['#2d5986'], linewidths=2.5, alpha=0.5)
        if arrow == True:
            for e in E:
                add_arrow(e, ax, direction='right', color='#2d5986')
        ax.add_collection(edge_collections)

        if path:
            edge_collections_path = mcol.LineCollection(E_path, colors=['#8f246b'], linewidths=3, alpha=0.7)
            if arrow == True:
                for e in E_path:
                    add_arrow(e, ax, direction='right', color='#8f246b')
            ax.add_collection(edge_collections_path)


        ax.set_ylim((min_y - 3, max_y + 3))
        ax.set_xlim((min_x, max_x))

        ax.set_xticks(np.linspace(int(min_x), int(max_x), 11))
        ax.set_xlabel("t", fontname='Ubuntu', fontsize=fontsize, color='#476b6b')
        for place,spine in plt.gca().spines.items():
            if place != 'bottom':
                spine.set_visible(False)
            else:
                spine.set_bounds(self.times[0],self.times[1])
                spine.set_color('#476b6b')
        ax.tick_params(right=False, left=False, labelleft=False,
                       labelbottom=True,
                       labelsize=fontsize,
                       colors='#476b6b')
        plt.tight_layout()

        return fig

    def size(self):
        return len(self.c_nodes)

    def core_number(self, storage_path=None, n_jobs=-1):
        L = defaultdict(list)

        def para_cores(comp, storage_path):
            # return comp.get_kcores(storage_path)
            return comp.core_number()

        r = Parallel(n_jobs=n_jobs, mmap_mode='r+')(
            delayed(para_cores)(comp, storage_path) for comp in self.c_nodes if comp.size() > 1)
        for l in r:
            for k, v in l.items():
                L[k] += v
        return L

    def k_core(self, k, storage_path=None, n_jobs=-1):
        L = []

        def para_cores(comp, storage_path):
            return comp.k_core(k)

        r = Parallel(n_jobs=n_jobs, mmap_mode='r+')(
            delayed(para_cores)(comp, storage_path) for comp in self.c_nodes if comp.size() > 1)
        for l in r:
            L += l
        return L

    def all_cliques(self, storage_path=None, n_jobs=-1):
        L = defaultdict(list)

        def para_cliques(comp, storage_path):
            # return comp.get_kcliques(storage_path)
            return comp.all_cliques()

        r = Parallel(n_jobs=n_jobs, mmap_mode='r+')(
            delayed(para_cliques)(comp, storage_path) for comp in self.c_nodes if comp.size() > 1)
        for l in r:
            for k, v in l.items():
                L[k] += v
        return L

    def k_clique(self, k, storage_path=None, n_jobs=-1):
        L = []

        def para_cliques(comp, storage_path):
            # return comp.get_kcliques(storage_path)
            return comp.k_clique(k)

        r = Parallel(n_jobs=n_jobs, mmap_mode='r+')(
            delayed(para_cliques)(comp, storage_path) for comp in self.c_nodes if comp.size() > 1)
        for l in r:
            L += l
        return L

    def get_stable_parts_inplace(self):
        L = []
        for comp in self.c_nodes:
            L += comp.get_stable_parts()
        return L

    def get_stable_dag(self):
        # Rajouter les stables parts comme des components
        stable_DAG = sdag.stable_dag()
        stable_DAG.set_id(self.id)
        stable_DAG.times = self.times
        cnt_c_nodes = 0
        for comp in self.c_nodes:
            stable_comps = comp.get_stable_components(format="object")
            for c in stable_comps:
                c.id = cnt_c_nodes
                cnt_c_nodes += 1
            stable_DAG.add_nodes(stable_comps)
            # new_cnt_c_nodes = cnt_c_nodes + len(stable_comps)
            # if new_cnt_c_nodes - cnt_c_nodes > 1:
            #     for j in range(cnt_c_nodes, new_cnt_c_nodes - 1):
            #         stable_DAG.add_link((j, j + 1))
            # cnt_c_nodes = new_cnt_c_nodes
        return stable_DAG

    def compute_links_inplace(self):

        chrono = time.time()
        # IF there's is some links, discard them.
        if self.c_links:
            self.c_links = []
        dict_begin_time_to_cnodes = defaultdict(list)
        dict_end_time_to_cnodes = defaultdict(list)

        if type(self.c_nodes[0]) == list:
            for i, c in self.id_scc_to_scc.items():
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
                                cn_child_nodes = set([c[2] for c in self.id_scc_to_scc[cn_child_id]])
                                cn_parent_nodes = set([c[2] for c in self.id_scc_to_scc[cn_parent_id]])
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
                                cn_child = self.id_scc_to_scc[cn_child_id]
                                cn_parent = self.id_scc_to_scc[cn_parent_id]
                                if not cn_child.nodes.isdisjoint(cn_parent.nodes):
                                    set_links.add((cn_parent_id, cn_child_id))
                                    a_l[cn_parent_id].add(cn_child_id)
        #
        print("[DAG] compute link inplace :", time.time() - chrono)
        chrono = time.time()
        # Remove useless links:
        to_remove = set()
        for l in set_links:
            i,j = l
            for k in a_l[i]:
                if k in a_l:
                    if j in a_l[k]:
                        to_remove.add((i,j))
                        break
        set_links -= to_remove
        print("[DAG] remove useless links :", time.time() - chrono)
        self.c_links = list(set_links)


    def compute_links(self, path_scc, id_wcc, storage_path_scc_dag=None):
        # IF there's is some links, discard them.

        dict_begin_time_to_cnodes = defaultdict(list)
        dict_end_time_to_cnodes = defaultdict(list)
        dict_id_set_nodes = defaultdict(set)

        # Get times and set of nodes corresponding to each c_nodes
        for filename in pathlib.Path(path_scc).iterdir():
            id_scc = int(re.search(r'\d+', os.path.basename(filename)).group())
            with open(filename, 'rb') as input:
                unpacker = msgpack.Unpacker(input, use_list=False)
                times = unpacker.unpack()
                dict_begin_time_to_cnodes[times[0]].append(id_scc)
                dict_end_time_to_cnodes[times[1]].append(id_scc)
                for i in unpacker:
                    if len(i) == 2:
                        u, v = i
                    else:
                        u, v = i[0]
                    dict_id_set_nodes[id_scc].add(u)
                    dict_id_set_nodes[id_scc].add(v)
        # For each time_end_find a time_begin (if it exists) then compute the intersection of ndoes
        for te in dict_end_time_to_cnodes:
            if te in dict_begin_time_to_cnodes:
                # Double for loop in order to match pairwise
                for cn_parent in dict_end_time_to_cnodes[te]:
                    for cn_child in dict_begin_time_to_cnodes[te]:
                        if not dict_id_set_nodes[cn_child].isdisjoint(dict_id_set_nodes[cn_parent]):
                            self.c_links.append((cn_parent, cn_child))
        # print("C Links :",self.c_links)
        if storage_path_scc_dag:
            self.store_links(storage_path_scc_dag, id_wcc)

    def adjacency_list(self):
        if not self.a_l:
            a_l = defaultdict(list)
            for l in self.c_links:
                # ONLY IF the destination is accesible from begin time
                # c = self.index_id_scc_to_nodes[l[1]]
                # if c.times[0] >= threshold:
                a_l[l[0]].append(l[1])
            self.a_l = a_l
        return self.a_l

    ################################
    #       FORMAT                 #
    ################################

    def cluster_to_object(self):
        new_cnodes = []
        for id_cc, cc in self.id_scc_to_scc.items():
            assert type(cc) == list
            new_cnodes = comp.strongly_connected_component(id=id_cc, times=(cc[0][0], cc[0][1]),
                                                           nodes=set([c[2] for c in cc]))
        self.c_nodes = new_cnodes
        self.id_scc_to_scc = {cc.id: cc for cc in new_cnodes}

    ###############################
    #       Paths Methods         #
    ###############################

    def path_induced_substream(self, path, node_to_label=None,
                               path_bounds=False):
        '''
        Transform a path in the condensation dag into a substream
        :param path: Sequence of SCC (object) in the condensation DAG
        :return:
        '''

        if type(path[0]) is int:
            path = [self.id_scc_to_scc[id_scc] for id_scc in path]

        new_nodes = {}
        nodes_to_new_nodes = defaultdict(lambda: len(nodes_to_new_nodes))
        new_node_to_label = {}
        new_node_to_id = {}
        new_links = {}
        t0_min, t1_max = math.inf, -math.inf
        for c in path:
            # print("Comp :",c.id)
            # print("c times :",c.times)
            # print("comp nodes :",c.nodes)
            # print("comp links :")
            # if c.links:
            #     for l in c.links:
            #         print(l)
            #     print()
            t0, t1 = c.times # Initial Bounds
            if path_bounds:
                t0, t1 = max(path_bounds[0],t0), min(path_bounds[1],t1) # Bounds
            t0_min,t1_max = min(t0, t0_min), max(t1, t1_max)

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
                    lt0,lt1 = l[0],l[1]
                    if lt1 < t0 or lt0 > t1:
                        # We do not consider links that end or begin outside the bounds
                        continue
                    lt0,lt1 = max(t0,lt0),min(t1,lt1)
                    u,v = l[2],l[3]
                    new_u, new_v = nodes_to_new_nodes[u], nodes_to_new_nodes[v]
                    if (new_u,new_v) in new_links:
                        if lt0 <= new_links[(new_u,new_v)][-1]:
                            new_links[(new_u,new_v)][-1] = lt1
                        else:
                            new_links[(new_u,new_v)] += [lt0,lt1]
                    else:
                        new_links[(new_u,new_v)] = [lt0,lt1]
                    # print("l :",(node_to_label[new_u],node_to_label[new_v]))

        F = sg.stream_graph(times=[t0_min, t1_max],
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
        '''
        return the Strongly Connected component containing the temporal source *node*.
        :param node:
        :return:
        '''
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
        '''
        Return the time to reach the *destination* from the temporal source node *source* in the SG.
        :param source:
        :param desination:
        :return:
        '''
        ttr = math.inf
        a_l = self.adjacency_list()
        start_comp = self.temporal_node_to_scc(source)
        if start_comp is None:
            return ttr
        st = source[0]
        print("Start comp:", start_comp)
        queue = deque([start_comp])  # comp
        visited = set([start_comp.id])
        while queue:
            comp = queue.popleft()
            if destination[2] in comp.nodes and \
                    destination[0] <= comp.times[0] <= comp.times[1] <= destination[1]:
                ttr = max(min(ttr, comp.times[0] - st), 0)
            if comp.id in a_l and a_l[comp.id]:
                for c_id in a_l[comp.id]:
                    if c_id not in visited:
                        c = self.id_scc_to_scc[c_id]
                        if c.times[0] <= st + ttr:
                            queue.append(c)
                        visited.add(c_id)
        return ttr

    def latency(self, source, destination):

        if type(self.c_nodes[0]) == list:
            self.cluster_to_object()

        if type(source) is int:
            return self._latency(source, destination)
        else:
            return self._latency_temporal_nodes(source, destination)

    def _latency(self, source, destination):
        # TODO : To finish
        return

    def _latency_temporal_nodes(self, source, destination):
        '''
        Return the latency between the temporal node *source* and the temporal node *destination* in the SG.
        :param source:
        :param destination:
        :return:
        '''
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
            start_comp = self.id_scc_to_scc[unvisited.pop()]
            st = start_comp.times[1]  # starting time
            queue = deque([(start_comp, st)])  # comp, starting time
            visited = set([start_comp.id])
            while queue:
                comp, st = queue.popleft()
                # We can reset the starting time
                if source in comp.nodes and source[0] <= comp.times[0] <= comp.times[1] <= source[1]:
                    st = comp.times[1]
                    unvisited.discard(comp.id)
                # Update latencies
                if destination[2] in comp.nodes and destination[0] <= comp.times[0] <= comp.times[1] <= destination[1]:
                    latency = max(min(latency, comp.times[0] - st), 0)

                if comp.id in a_l and a_l[comp.id]:
                    for c_id in a_l[comp.id]:
                        if c_id not in visited:
                            c = self.id_scc_to_scc[c_id]
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
        for n in self.node_to_id_scc:
            potential_ttr = [ttr_comp[c] for c in self.node_to_id_scc[n] if c in ttr_comp]
            if potential_ttr:
                ttr[n] = min(potential_ttr)
        ttr[source] = 0
        return ttr

    def _times_to_reach_ss(self, source):
        id_start_comp = self.node_to_id_scc[source][0]  # The first SCC where the source appears
        start_comp = self.id_scc_to_scc[id_start_comp]
        ttr_comp = self._times_to_reach_comp_ss(id_start_comp, start_time=start_comp.times[0])
        ttr = self.postprocess_ttr(source, ttr_comp)
        return ttr

    def _times_to_reach_temporal_nodes_ss(self, source):
        potential_start_comp = self.segmented_node_to_id_scc[source]
        # #
        # if source == (1386320340.0, 1386329460.0, 1):
        #     for id_comp in potential_start_comp:
        #         if self.id_scc_to_scc[id_comp].times[0] == 1386320340.0:
        #             print(self.id_scc_to_scc[id_comp].nodes)
        #         if 236 in self.id_scc_to_scc[id_comp].nodes:
        #             print(self.id_scc_to_scc[id_comp].times)
        #
        # #

        for id_comp in potential_start_comp:
            comp_t0, comp_t1 = self.id_scc_to_scc[id_comp].times
            if comp_t0 <= source[0] <= comp_t1:
                id_start_comp = id_comp
                break
        ttr_comp = self._times_to_reach_comp_ss(id_start_comp, start_time=source[0])
        ttr = self.postprocess_ttr(source[2], ttr_comp)
        return ttr


    def _times_to_reach_comp_ss(self, id_start_comp, start_time):
        '''
        Return the times to reach from the temporal source node :*source* to every other node in the SG.
        :param source:
        :return:
        '''
        a_l = self.adjacency_list()
        times_to_reach_comp = {id_start_comp: 0}
        queue = deque([id_start_comp])
        visited = {id_start_comp}
        while queue:
            id_comp = queue.popleft()
            comp = self.id_scc_to_scc[id_comp]
            times_to_reach_comp[id_comp] = comp.times[0] - start_time
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
        for n in self.node_to_id_scc:
            potential_latencies = [latencies_comp[c] for c in self.node_to_id_scc[n] if c in latencies_comp]
            if potential_latencies:
                latencies[n] = min(potential_latencies)
        latencies[source] = 0
        return latencies


    def _latencies_ss(self, source):
        ids_start_comp = self.node_to_id_scc[source]  # SCCs where the source appears
        latencies_comp = self._latencies_comp_ss(ids_start_comp)
        latencies = self.postprocess_latencies(latencies_comp, source)
        return latencies


    def _latencies_temporal_nodes_ss(self, source):
        ids_start_comp = self.segmented_node_to_id_scc[source]  # SCC where the source appears
        latencies_comp = self._latencies_comp_ss(ids_start_comp)
        latencies = self.postprocess_latencies(latencies_comp, source[2])
        return latencies

    def _latencies_comp_ss(self, ids_start_comp):
        '''
        Return the latencies from the temporal source node: *source* to every other node in the SG.
        :param source:
        :return:
        '''
        a_l = self.adjacency_list()
        latencies = {i: 0 for i in ids_start_comp}  # clefs : destination nodes ; valeurs : latency
        set_start_comps = set(ids_start_comp)
        unvisited = set(ids_start_comp)
        visited_to_st = {}
        # TODO : unvisited and set_start _comps necessaire ?
        while len(unvisited) != 0:
            id_start_comp = unvisited.pop()
            st = self.id_scc_to_scc[id_start_comp].times[1]  # starting time
            visited_to_st[id_start_comp] = st
            queue = deque([(id_start_comp, st)])  # comp, starting time
            while queue:
                id_comp, st = queue.popleft()
                comp = self.id_scc_to_scc[id_comp]
                if id_comp in set_start_comps:
                    # We can reset the starting time
                    st = comp.times[1]
                    unvisited.discard(id_comp)
                else:
                    if id_comp in latencies:
                        latencies[id_comp] = min(comp.times[0] - st, latencies[id_comp])
                    else:
                        latencies[id_comp] = comp.times[0] - st

                if id_comp in a_l and a_l[id_comp]:
                    for c_id in a_l[id_comp]:

                        if c_id not in visited_to_st or st > visited_to_st[c_id]:
                            # TODO : Verifier que la condition s'applique bien

                            # On ne doit pas avoir c_id plusieurs fois dans la queue !!
                            # Ou tester en sortie de pile !
                            # We leave later :)
                            queue.append((c_id, st))
                            visited_to_st[c_id] = st
        return latencies


    ##############################################
    #   3. Arbitrary Foremost and fastest Paths  #
    ##############################################

    def foremost_path(self, source, destination):
        '''
        Compute a foremost path between 'source' and 'destination" starting at 'begin time'.
        '''
        st = source[0]
        start_comp = self.temporal_node_to_scc(source)
        a_l = self.adjacency_list()

        # Custom BFS on DAG
        def bfs_scc(a_l, start_comp, destination):
            path_queue = deque([(start_comp, [start_comp.id])])
            ttr = math.inf  # This variable, once assigned, is used as a threshold (yeah)
            visited = set([start_comp.id])
            while path_queue:
                (comp, path) = path_queue.popleft()
                # print(" v : ",v)
                # print(" len path queue :",len(path_queue))
                if comp.id in a_l and a_l[comp.id]:
                    for c_id in a_l[comp.id]:
                        if c_id not in visited:
                            c = self.id_scc_to_scc[c_id]
                            visited.add(c_id)
                            # print(" comp nodes :",comp.nodes)
                            if destination[2] in c.nodes and destination[0] <= c.times[0] <= c.times[1] <= destination[
                                1]:
                                # print("PATH Found :",path)
                                ttr = max(min(ttr, c.times[0] - st), 0)
                                yield path + [c_id]
                            elif c.times[0] <= st + ttr:
                                path_queue.append((c, path + [c_id]))

        # print(" Start BFS")
        foremost_paths = list(bfs_scc(a_l, start_comp, destination))
        if not foremost_paths:
            return None, None

        path_times = []
        for p in foremost_paths:
            t = []
            for c_id in p:
                c = self.id_scc_to_scc[c_id]
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
        '''
        Compute the fastest path between 'source' and 'destination'
        :param source:
        :param destination:
        :return:
        '''
        a_l = self.adjacency_list()
        ids_start_comp = set(self.segmented_node_to_id_scc[source])
        unvisited = copy.copy(ids_start_comp)  # SCC where the source appears
        visited_to_st = {}
        latency = math.inf
        fastest_paths = None
        while len(unvisited) != 0:
            id_start_comp = unvisited.pop()
            start_comp = self.id_scc_to_scc[id_start_comp]
            st = start_comp.times[1]  # starting time
            visited_to_st[id_start_comp] = st
            path_queue = deque([((id_start_comp, st), [id_start_comp])]) # ((id comp, start time),path)
            while path_queue:
                e, path = path_queue.popleft()
                id_comp, st = e
                comp = self.id_scc_to_scc[id_comp]

                if id_comp in ids_start_comp:
                    # We can reset the starting time
                    st = comp.times[1]
                    path = [id_comp]
                    unvisited.discard(id_comp)
                if destination[2] in comp.nodes and destination[0] <= comp.times[0] <= comp.times[1] <= destination[1]:
                    new_latency = max(comp.times[0] - st, 0)
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

    def compute_all_foremost_paths(self, dict_id_wcc_to_dag, index_node_to_scc, source, destination, dict_offset=None,
                                   duration_threshold=None, start_comp=None):
        '''
        Compute a foremost path between 'source' and 'destination" starting at 'begin time'.
        '''
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
                c = dict_id_wcc_to_dag[id_wcc].id_scc_to_scc[id_scc]
                if c.times[0] <= start <= c.times[1]:
                    if destination in c.nodes:
                        return [[id_scc]], 0, id_wcc
                    start_comp = id_scc
                    break

        if start_comp is None:
            print("Node : " + str(source) + " does not exist in the Stream Graph at time " + str(start) + " !")
        a_l = defaultdict(list)
        end_time_comp = []
        G = dict_id_wcc_to_dag[id_wcc]
        if G.c_links:
            for l in G.c_links:
                # ONLY IF the destination is accesible from begin time
                c = G.id_scc_to_scc[l[1]]
                if c.times[0] >= start:
                    a_l[l[0]].append(l[1])
                    end_time_comp.append(c.times[1])

        if duration_threshold is None:
            if not end_time_comp:
                return None, None, None
            duration_threshold = max(end_time_comp) + 1 - start

        # print("duration threshold :",duration_threshold)
        # print(" a_l : ", a_l)

        # Custom BFS on DAG
        def bfs_scc(a_l, start_comp, destination):
            path_queue = [(start_comp, [start_comp])]
            foremost_duration = duration_threshold  # This variable, once assigned, is used as a threshold (yeah)
            while path_queue:
                (v, path) = path_queue.pop(0)
                # print(" v : ",v)
                if len(path_queue) > 0 and len(path_queue) % 10000 == 0:
                    print(" len path queue :", len(path_queue))
                if v in a_l and a_l[v]:
                    for c in a_l[v]:
                        comp = dict_id_wcc_to_dag[id_wcc].id_scc_to_scc[c]
                        print(" comp id :", c)
                        print(" comp times :", comp.times)
                        if destination in comp.nodes:
                            # print("PATH Found :",path)
                            foremost_duration = min(foremost_duration, comp.times[0] - start)
                            yield path + [c]
                            # return path + [c]
                        elif comp.times[0] <= start + foremost_duration:
                            path_queue.append((c, path + [c]))

        # print(" Start BFS")
        foremost_paths = list(bfs_scc(a_l, start_comp, destination))
        #
        # foremost_path = bfs_scc(a_l, start_comp, destination)
        # if foremost_path:
        #     return [foremost_path],dict_id_wcc_to_dag[id_wcc].index_id_scc_to_nodes[foremost_path[-1]].times[0]-start,id_wcc
        # else:
        #     return None, None, None

        if not foremost_paths:
            return None, None, None

        path_times = []
        for p in foremost_paths:
            t = []
            for c in p:
                t.append(dict_id_wcc_to_dag[id_wcc].id_scc_to_scc[c].times[0])
            path_times.append(t)
        # print("Path times :", path_times)
        min_t = min([p[-1] for p in path_times])
        # print("Foremost time : ", min_t)
        fm_paths = []
        for p, t in zip(foremost_paths, path_times):
            if t[-1] == min_t:
                fm_paths.append(p)
        return fm_paths, min_t - start, id_wcc


    def all_fastest_paths(self,source,destination):
        '''
        Compute the fastest path between 'source' and 'destination'
        :param source:
        :param destination:
        :return:
        '''
        a_l = self.adjacency_list()
        ids_start_comp = set(self.segmented_node_to_id_scc[source])
        unvisited = copy.copy(ids_start_comp)  # SCC where the source appears
        visited_to_st = {} # We store visited component along with the starting time corresponding with the path that reached them;
        latency = math.inf
        fastest_paths = None
        while len(unvisited) != 0:
            id_start_comp = unvisited.pop()
            start_comp = self.id_scc_to_scc[id_start_comp]
            st = start_comp.times[1]  # starting time
            visited_to_st[id_start_comp] = st
            path_queue = deque([((id_start_comp, st), [id_start_comp])]) # ((id comp, start time),path)
            while path_queue:
                e, path = path_queue.popleft()
                id_comp, st = e
                comp = self.id_scc_to_scc[id_comp]

                if id_comp in ids_start_comp:
                    # We can reset the starting time
                    st = comp.times[1]
                    path = [id_comp]
                    unvisited.discard(id_comp)
                if destination[2] in comp.nodes and destination[0] <= comp.times[0] <= comp.times[1] <= destination[1]:
                    new_latency = max(comp.times[0] - st, 0)
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



    def all_fastest_paths_ss(self,source):
        '''
        Compute the fastest path between 'source' and all other nodes (in condensation DAG)

        :param source:
        :return:
        '''
        a_l = self.adjacency_list()
        ids_start_comp = set(self.segmented_node_to_id_scc[source])
        unvisited = copy.copy(ids_start_comp)  # SCC where the source appears
        latencies_comp = {i: 0 for i in ids_start_comp}  # keys : destination nodes ;values : latency
        fastest_paths = {i:{i} for i in ids_start_comp}
        while len(unvisited) != 0:
            id_start_comp = unvisited.pop()
            start_comp = self.id_scc_to_scc[id_start_comp]
            st = start_comp.times[1]  # starting time
            path_queue = deque([((id_start_comp,st), [id_start_comp])]) # ((id comp,start_time),path)
            while path_queue:
                e, path = path_queue.popleft()
                id_comp,st = e
                comp = self.id_scc_to_scc[id_comp]
                new_latency = comp.times[0] - st
                if id_comp in ids_start_comp:
                    # We can reset the starting time
                    st = comp.times[1]
                    path = [id_comp]
                    unvisited.discard(id_comp)
                else :
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
                        comp = self.id_scc_to_scc[c_id]
                        if c_id not in latencies_comp or latencies_comp[c_id] >= comp.times[0]-st:
                            path_queue.append(((c_id,st), path + [c_id]))

        fastest_paths_nodes, latencies_nodes = self.postprocess_fastest_paths(fastest_paths,latencies_comp)
        del fastest_paths_nodes[source[2]]
        return fastest_paths_nodes, latencies_nodes

    def postprocess_fastest_paths(self,fastest_paths_comp,latencies_comp):
        # FP and Latencies in DAG to FP and Lat in Stream :
        latencies_nodes = defaultdict(lambda:math.inf)
        fastest_paths_nodes = {}
        for n in self.node_to_id_scc:
            lat = math.inf
            for c in self.node_to_id_scc[n]:
                if c in latencies_comp:
                    if latencies_comp[c] < lat:
                        fastest_paths_nodes[n] = fastest_paths_comp[c].copy()
                        latencies_nodes[n] = latencies_comp[c]
                        lat = latencies_comp[c]
                    elif latencies_comp[c] == lat:
                        fastest_paths_nodes[n] |= fastest_paths_comp[c]

        return fastest_paths_nodes,latencies_nodes
    ###############################################
    #   5. Shortest Fastest Path                  #
    #    (Hybrids Methods)                        #
    ###############################################

    def shortest_fastest_path(self,source,destination):
        '''
        TODO : FINISH
        :param source:
        :param destination:
        :return:
        '''
        latencies = {}
        lengths = {}
        list_fp,l = self.all_fastest_path(source, destination)
        latencies[destination] = l
        for fp in list_fp:
            substream = self.path_induced_substream(fp)
            D = substream.distances(source)
            lengths[destination] = min(D[destination],lengths[destination])
        return latencies,lengths


    def shortest_fastest_path_ss(self,source,node_to_label=None):
        '''
        Compute SFP in DAG !!
        :param source:
        :return:
        '''
        lengths = defaultdict(lambda:math.inf)
        fastest_paths_nodes, latencies_nodes = self.all_fastest_paths_ss(source)
        lengths[source[2]] = 0

        for destination in fastest_paths_nodes:
            list_fp = fastest_paths_nodes[destination]
            for fp in list_fp:
                if type(fp) == int:
                    fp = [fp]

                # print("\nfp :",fp)
                # print("destination :",destination)
                # assert source[2] in self.id_scc_to_scc[fp[0]].nodes
                # assert destination in self.id_scc_to_scc[fp[-1]].nodes
                # If length of the cdag path > 1 : last time of the first comp; first time of the last comp
                # Because they are SCC and all stream nodes are reachables at any instant.
                path_bounds = False
                if len(fp) > 1:
                    path_bounds = (self.id_scc_to_scc[fp[0]].times[1],self.id_scc_to_scc[fp[-1]].times[0])
                    # print("Path_bounds :",path_bounds)
                    # print("latencies destination :",latencies_nodes[destination])
                    # assert path_bounds[1] - path_bounds[0] == latencies_nodes[destination]

                # else:
                #     assert latencies_nodes[destination] == 0
                #     print("Path_bounds :",path_bounds)
                substream = self.path_induced_substream(fp, path_bounds=path_bounds,
                                                        node_to_label=node_to_label)  #  Get Substream for DAG paths
                stream_node_to_substream_node = {v : k for k,v in substream.node_to_id.items()}
                substream_source = stream_node_to_substream_node[source[2]]
                substream_destination = stream_node_to_substream_node[destination]
                # assert destination in stream_node_to_substream_node
                # assert source[2] in stream_node_to_substream_node
                # Get source in Substream
                b,e = None,None
                npres = substream.node_presence[substream_source]
                for nt0,nt1 in zip(npres[::2],npres[1::2]):
                    if source[0]<=nt0<=nt1<=source[1]:
                        b,e = nt0,nt1
                        break
                # print("Source in stream :",source)
                # print("Source in substream :",(b,e,substream_source))
                # print("Destination in substream :",substream_destination)
                # substream.plot()
                # plt.show()
                # Compute distances in Substream
                # print("Latency in Stream :",latencies_nodes[destination])
                # L,D = substream.latencies_and_lengths(((b, e, substream_source)))
                D = substream.distances(((b, e, stream_node_to_substream_node[source[2]])))
                # if substream_destination not in D:
                #     print("nb wcc substream :",substream.number_of_weakly_connected_component())
                #     wcc = substream.weakly_connected_components()
                #     a_l = substream.instant_graph(substream.times[0], label=False)
                #     plot_adjacency_list(substream, a_l, label = False)
                #     substream.plot()
                #     plt.show()


                # if fp == [8, 25, 12, 13] and destination == 3 and source[2] == 7:
                #     # a_l = substream.instant_graph(substream.times[0], label=False)
                #     # plot_adjacency_list(substream, a_l, label = False)
                #     self.plot()
                #     substream.plot()
                #     plt.show()

                # print("Latency in substream :",L[substream_destination])
                # print("length in substream :",D[substream_destination])
                # assert L[substream_destination] == latencies_nodes[destination]
                lengths[destination] = min(lengths[destination],D[substream_destination])
            # print("Length SFP HYBRID :",lengths[destination])
        return latencies_nodes,lengths

    ################### END PATHS #################################

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
                           col="#8f246b", label_nodes_by_letters=False):
        '''
        Draw a path on the Stream Graph S
        :param S: A Stream Graph
        :return:
        '''

        if type(path[0]) is int:
            path = [self.id_scc_to_scc[id_scc] for id_scc in path]

        S = self.path_induced_substream(path)
        S.plot(label_nodes_by_letters=label_nodes_by_letters)
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
        plt.tight_layout()

    def plot_fastest_path(self, path, latency, source, destination,
                          col="#8f246b", label_nodes_by_letters=False):
        '''
        Draw a path on the Stream Graph S
        :param S: A Stream Graph
        :return:
        '''
        if type(path[0]) is int:
            path = [self.id_scc_to_scc[id_scc] for id_scc in path]
        S = self.path_induced_substream(path)
        S.plot(label_nodes_by_letters=label_nodes_by_letters)
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
        plt.tight_layout()

    def refactor_path(self, path_comp, path_times, source, destination):
        print("Path comp :", path_comp)
        print("Path Time :", path_times)
        P = pt.path([], [])
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

    def plot(self, colors=None, saving_path=None, format='pdf', title=None,
             node_to_label=None, fontsize=26):
        fig, ax = plt.subplots()
        rectangles = []
        labels = []
        unique_nodes = set()
        min_t = math.inf
        max_t = -math.inf
        icol = 0
        if not colors:
            colors = sg.get_cmap(len(self.c_nodes))
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
            eps = random.choice([1, -1]) * (random.random() / 10)
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
        ax.set_xlim((min_t, max_t))
        ax.set_ylim((min(unique_nodes), max(unique_nodes) + 1))

        # if label_nodes_by_letters == True and len(unique_nodes) <= 26:
        #     plt.yticks(range(len(unique_nodes)), ascii_uppercase[0:len(unique_nodes)], fontname='Ubuntu',
        #                fontsize=12, color='#666699')
        if node_to_label:
            plt.yticks(range(len(unique_nodes)),
                       [node_to_label[i] for i in range(len(unique_nodes))],
                       fontname='Ubuntu',
                       fontsize=fontsize, color='#666699')
        else:
            plt.yticks(range(len(unique_nodes)), range(len(unique_nodes)),
                       fontname='Ubuntu', fontsize=fontsize, color='#666699')
        plt.ylim((-1, len(unique_nodes)))
        plt.ylabel("Nodes", fontname='Ubuntu', fontsize=fontsize, color='#666699')

        plt.xticks(np.linspace(int(min_t), int(max_t), 11),
                   fontname='Ubuntu', fontsize=fontsize, color='#476b6b')
        plt.xlabel("t", fontname='Ubuntu', fontsize=fontsize, color='#476b6b')
        if title:
            plt.title(title, fontname='Ubuntu', fontsize=fontsize)
        for place,spine in plt.gca().spines.items():
            if place != 'bottom':
                spine.set_visible(False)
            else:
                spine.set_bounds(self.times[0],self.times[1])
                spine.set_color('#476b6b')

        plt.tick_params(top=False, bottom=True, right=False, left=True, labelleft=True, labelbottom=True)
        plt.tight_layout()

        if saving_path:
            fig.savefig(saving_path + "." + format, format=format)
        return fig


def compute_links_scf(path_scc_scf, dict_scc_dag, storage_path_scc_dag):
    id_wcc_old = None
    cnt_lines = 0
    t_begin = time.time()
    dict_scc_to_lineno = {}
    dict_begin_times_to_scc = defaultdict(list)
    dict_end_times_to_scc = defaultdict(list)
    with open(path_scc_scf + 'comp.parquet', 'rb') as input, open(storage_path_scc_dag + 'scc_dag_with_links.scf',
                                                                 'wb') as output_file:
        pq_reader = pa.RecordBatchFileReader(input)
        n_scc = pq_reader.num_record_batches
        print("N SCC :", n_scc)
        for i in range(n_scc):  # SCC components are ordered by their end time
            current_batch = pq_reader.get_batch(i).to_pydict()
            id_wcc, id_scc = current_batch['id'][0]
            t0, t1 = current_batch['times'][0]
            links = current_batch['links'][0]
            if id_wcc in dict_scc_dag:
                if not id_wcc_old:  # INIT
                    id_wcc_old = id_wcc
                if (id_wcc % 100000 == 0 and id_wcc != 0) or (id_scc % 100000 == 0 and id_scc != 0):
                    print("id wcc :", id_wcc, "id comp :", id_scc)
                    print("comp :", t0, t1, len(links), " t :", time.time() - t_begin)
                if id_wcc != id_wcc_old:
                    # When we change of wcc we reset data structure
                    S = dict_scc_dag[id_wcc_old]
                    # For each time_end_find a time_begin (if it exists) then compute the intersection of nodes
                    # print("end times :",dict_end_times_to_scc)
                    # print("begin times :",dict_begin_times_to_scc)
                    for te in dict_end_times_to_scc:
                        if te in dict_begin_times_to_scc:
                            # CATCH nodes in comp for dict_end_times_to_cnodes[te] and dict_begin_times_to_cnodes[te]
                            list_current_scc = dict_end_times_to_scc[te] + dict_begin_times_to_scc[te]
                            dict_id_set_nodes = load_sccs_from_lineno(pq_reader, list_current_scc, dict_scc_to_lineno)
                            for cn_parent in dict_end_times_to_scc[te]:
                                for cn_child in dict_begin_times_to_scc[
                                    te]:  # Double for loop in order to match pairwise
                                    if not dict_id_set_nodes[cn_child].isdisjoint(dict_id_set_nodes[cn_parent]):
                                        S.c_links.append((cn_parent, cn_child))
                    S.store(output_file)  # Store the dag.
                    # print("STORE S :", id_wcc_old)
                    del dict_scc_dag[id_wcc_old]
                    id_wcc_old = id_wcc
                    dict_scc_to_lineno = {}
                    dict_begin_times_to_scc = defaultdict(list)
                    dict_end_times_to_scc = defaultdict(list)

                dict_scc_to_lineno[id_scc] = cnt_lines
                dict_begin_times_to_scc[t0].append(id_scc)
                dict_end_times_to_scc[t1].append(id_scc)

            cnt_lines += 1

        print("LAST ITERATION :", id_wcc)
        if id_wcc in dict_scc_dag:
            S = dict_scc_dag[id_wcc]
            # For each time_end_find a time_begin (if it exists) then compute the intersection of ndoes
            # print("end times :", dict_end_times_to_scc)
            # print("begin times :", dict_begin_times_to_scc)
            for te in dict_end_times_to_scc:
                if te in dict_begin_times_to_scc:
                    # CATCH nodes in comp for dict_end_times_to_cnodes[te] and dict_begin_times_to_cnodes[te]
                    list_current_scc = dict_end_times_to_scc[te] + dict_begin_times_to_scc[te]
                    dict_id_set_nodes = load_sccs_from_lineno(pq_reader, list_current_scc, dict_scc_to_lineno)
                    for cn_parent in dict_end_times_to_scc[te]:
                        for cn_child in dict_begin_times_to_scc[te]:  # Double for loop in order to match pairwise
                            if not dict_id_set_nodes[cn_child].isdisjoint(dict_id_set_nodes[cn_parent]):
                                S.c_links.append((cn_parent, cn_child))
            S.store(output_file)
            print("STORE S :", id_wcc)


def load_sccs_from_lineno(pq_reader, list_current_scc, dict_scc_to_lineno):
    dict_id_set_nodes = defaultdict(set)
    for id_scc in list_current_scc:
        current_batch = pq_reader.get_batch(dict_scc_to_lineno[id_scc]).to_pydict()
        links = current_batch['links'][0]
        for j in links:
            dict_id_set_nodes[id_scc].add(j[0])
            dict_id_set_nodes[id_scc].add(j[1])
    return dict_id_set_nodes


def compute_links_msgpack(path_scc_scf, dict_scc_dag, storage_path_scc_dag):
    t_begin = time.time()
    id_wcc_old = None
    dict_begin_times_to_scc = defaultdict(list)
    dict_end_times_to_scc = defaultdict(list)
    scc_2_offset = {}
    ##############
    gc.disable()
    #############
    with open(path_scc_scf + 'comp.scf', 'rb') as input, open(storage_path_scc_dag + 'scc_dag_with_links.scf',
                                                             'wb') as output_file:
        U = msgpack.Unpacker(input, use_list=False)
        offset = U.tell()
        for i in U:  # SCC components are ordered by their end time
            id_wcc, id_scc, t0, t1 = i[0], i[1], i[2], i[3]
            if id_wcc in dict_scc_dag:
                if id_wcc_old is None:  # INIT
                    id_wcc_old = id_wcc
                if id_wcc % 100000 == 0 or (id_scc % 100000 == 0 and id_scc != 0):
                    print("id wcc :", id_wcc, "id comp :", id_scc)
                    print("comp :", t0, t1, len(i[4]), " t :", time.time() - t_begin)
                if id_wcc != id_wcc_old:
                    # When we change of wcc we reset data structure
                    S = dict_scc_dag[id_wcc_old]
                    # For each time_end_find a time_begin (if it exists) then compute the intersection of nodes
                    # print("end times :",dict_end_times_to_scc)
                    # print("begin times :",dict_begin_times_to_scc)
                    for te in dict_end_times_to_scc:
                        if te in dict_begin_times_to_scc:
                            # CATCH nodes in comp for dict_end_times_to_cnodes[te] and dict_begin_times_to_cnodes[te]
                            list_current_scc = dict_end_times_to_scc[te] + dict_begin_times_to_scc[te]
                            dict_id_set_nodes = load_set_nodes_from_msgpack(path_scc_scf, list_current_scc,
                                                                            scc_2_offset)
                            for cn_parent in dict_end_times_to_scc[te]:
                                for cn_child in dict_begin_times_to_scc[
                                    te]:  # Double for loop in order to match pairwise
                                    if not dict_id_set_nodes[cn_child].isdisjoint(dict_id_set_nodes[cn_parent]):
                                        S.c_links.append((cn_parent, cn_child))
                    S.store(output_file)  # Store the dag.
                    # print("STORE S :", id_wcc_old)
                    del dict_scc_dag[id_wcc_old]
                    id_wcc_old = id_wcc
                    dict_begin_times_to_scc = defaultdict(list)
                    dict_end_times_to_scc = defaultdict(list)
                    scc_2_offset = {}

                scc_2_offset[id_scc] = offset
                dict_begin_times_to_scc[t0].append(id_scc)
                dict_end_times_to_scc[t1].append(id_scc)

            offset = U.tell()

        # print("LAST ITERATION :", id_wcc)
        if id_wcc in dict_scc_dag:
            S = dict_scc_dag[id_wcc]
            # For each time_end_find a time_begin (if it exists) then compute the intersection of ndoes
            # print("end times :", dict_end_times_to_scc)
            # print("begin times :", dict_begin_times_to_scc)
            for te in dict_end_times_to_scc:
                if te in dict_begin_times_to_scc:
                    # print("current time :",te)
                    # CATCH nodes in comp for dict_end_times_to_cnodes[te] and dict_begin_times_to_cnodes[te]
                    list_current_scc = dict_end_times_to_scc[te] + dict_begin_times_to_scc[te]
                    # print("Current SCC :",list_current_scc)
                    dict_id_set_nodes = load_set_nodes_from_msgpack(path_scc_scf, list_current_scc, scc_2_offset)
                    for cn_parent in dict_end_times_to_scc[te]:
                        for cn_child in dict_begin_times_to_scc[te]:  # Double for loop in order to match pairwise
                            if not dict_id_set_nodes[cn_child].isdisjoint(dict_id_set_nodes[cn_parent]):
                                S.c_links.append((cn_parent, cn_child))
            S.store(output_file)
            # print("STORE S :", id_wcc)
        ##############
        gc.enable()
        ##############


def load_scc_from_msgpack(id_wcc, id_scc, scc_storage, dict_offset):
    gc.disable()
    with open(scc_storage + 'comp.scf', 'rb') as file_input:
        file_input.seek(dict_offset[(id_wcc, id_scc)])
        i = msgpack.Unpacker(file_input, use_list=False).unpack()
        # print("i :", i)
    gc.enable()
    return comp.strongly_connected_component(id=id_scc,
                                             times=(i[2], i[3]),
                                             nodes=set([u for l in i[4] for u in l[-2:]]),
                                             links=i[4])


def load_set_nodes_from_msgpack(storage_scc, list_current_scc, scc_2_offset):
    '''
    Get nodes inside each considered comp using an offset to fast random access (via peek())
    :param storage_scc: Storage of SCC (msgpack format (scf ?))
    :param list_current_scc: List of comp id's
    :param scc_2_offset: dictionary id_scc <-> offset (in the msgpack file)
    :return: A dictionary id_scc <-> nodes (inside the comp)
    '''
    dict_id_2_set_nodes = defaultdict(set)
    list_current_scc.sort()
    with open(storage_scc + "comp.scf", 'rb') as ipt:
        for id_scc in list_current_scc:
            ipt.seek(scc_2_offset[id_scc])
            i = msgpack.Unpacker(ipt).unpack()
            for j in i[4]:
                if len(j) == 2:
                    dict_id_2_set_nodes[id_scc].add(j[0])
                    dict_id_2_set_nodes[id_scc].add(j[1])
                else:
                    dict_id_2_set_nodes[id_scc].add(j[2])
                    dict_id_2_set_nodes[id_scc].add(j[3])
    return dict_id_2_set_nodes


if __name__ == '__main__':

    #

    __directory__ = "/home/leo/Dev/Data_Stream/" + "Socio_Patterns/Workplace/"
    __filedate__ = "Workplace"

    # __directory__ = "/home/leo/Dev/Data_Stream/"
    # __file__ = "generated"
    # __file__ ="example"
    # __file__ = "High_School_2013"

    __directory__ = "/home/leo/Dev/Data_Stream/" + "Socio_Patterns/High_School_2013/"
    __filedate__ = "High_School_2013"

    __directory__ = "/home/leo/Dev/Data_Stream/Crawdad/Rollernet/"
    __filedate__ = "rollernet"
    # __directory__ = "/home/leo/Dev/Data_Stream/Crawdad/Infocom/"
    # __file__ = "infocom"

    # __directory__ = "/home/leo/Dev/Data_Stream/Bitcoin/bitcoin_otc/"
    # __file__ = "bitcoin_otc"

    # __directory__ = "/home/leo/Dev/Data_Stream/Bitcoin/bitcoin_alpha/"
    # __file__ = "bitcoin_alpha"

    # __directory__ = "/home/leo/Dev/Data_Stream/CollegeMsg/"
    # __file__ = "CollegeMsg"

    # __directory__ = "/home/leo/Dev/Data_Stream/Enron/"
    # __file__ = "enron"

    # __directory__ = "/home/leo/Dev/Data_Stream/parse_net/11September/"
    # __file__ = "Days"

    # __directory__ = "/home/leo/Dev/Data_Stream/askubuntu/"  # LIP6
    # __file__ = "askubuntu" # Open With Crawdad

    # __directory__ = "/home/leo/Dev/Data_Stream/DBLP/"  # LIP6
    # __file__ = "dblp" # Open With Crawdad

    # __directory__ = "/home/leo/Dev/Data_Stream/example_paper/"
    # __file__ ="example_paper"
    #
    __directory__ = "/home/leo/Dev/Data_Stream/"
    __filedate__ = "sg_generated"



    # __directory__ = "/home/leo/Dev/Straph/examples/"
    # __file__ = "example_path"
    #
    # __directory__ = "/home/leo/Dev/Data_Stream/slashdot/"
    # __file__ = "slashdot-threads"

    # __directory__ = "/home/leo/Dev/Data_Stream/ITS/"
    # __file__ = "Suricata"
    # __directory__ = "/home/leo/Dev/Data_Stream/path_examples/"
    # __file__ = "path_example"

    storage_path_scc = __directory__ + __filedate__ + "_scc_storage/"
    storage_path_scc_dag = __directory__ + __filedate__ + "_scc_dag_storage/"

    pathlib.Path(storage_path_scc_dag).mkdir(parents=True, exist_ok=True)

    dict_scc_dag = load_scc_dag_scf(storage_path_scc_dag + 'scc_dag.scf')
    print("number of dag with at least 2 nodes :", len(dict_scc_dag))
    print("dict comp dag :", dict_scc_dag)
    print("total number of nodes :", sum([len(i.c_nodes) for i in dict_scc_dag.values()]))

    ##########################
    #  COMPUTE SCC DAG LINKS #
    ##########################
    profile_shit(None,
                 "compute_links_msgpack(storage_path_scc,dict_scc_dag,storage_path_scc_dag)")
    compute_dict_offset(storage_path_scc + 'comp.scf', storage_path_scc_dag + "dict_offset_scc.mpck")
    #
    # TODO : - Mettre le compute dict offset, une seule fois, garder les offsets une fois pour toute.
    # TODO : - Checker les foremost paths ~ pour workplace and high_school == tps minimal pour qu'une information transite.
    # TODO : - Rajouter + Checker Fastest Path.
    exit()
    ########################
    #       PATHS          #
    ########################
    dict_scc_dag = load_scc_dag_scf(storage_path_scc_dag + 'scc_dag_with_links.scf')
    # Get interesting WCC :
    for k, v in dict_scc_dag.items():
        if len(v.c_nodes) > 300:
            print("\nid_wcc : ", k)
            print("n links :", len(v.c_links))
            print("n nodes :", len(v.c_nodes))

    id_wcc = 7987
    # id_wcc = 1864

    print(dict_scc_dag[id_wcc].c_nodes)
    print(dict_scc_dag[id_wcc].c_links)
    dict_offset = msgpack.unpack(open(storage_path_scc_dag + 'dict_offset_scc.mpck', 'rb'), use_list=False)
    for id_scc in dict_scc_dag[id_wcc].c_nodes:
        c = load_scc_from_msgpack(id_wcc, id_scc, storage_path_scc, dict_offset)
        print(" id comp :", id_scc)
        for n in c.nodes:
            print("n :", n)
    # paths,duration =dict_scc_dag[id_wcc].foremost_path(source=1412,destination=1218,begin_time=1385982620,
    #                                    scc_storage=storage_path_scc,dict_offset_scc=dict_offset_scc)
    # print("Foremost Paths :",paths)
    # print("Duration :",duration)

    paths, latency = dict_scc_dag[id_wcc].fastest_path(source=205, destination=275,
                                                       scc_storage=storage_path_scc, dict_offset=dict_offset)
    print("Fastest Paths :", paths)
    print("Nb Fastest Paths :", len(paths))
    print("Latency :", latency)

    # load_dag_from_ids(1,0,storage_path_scc+'comp.scf')
    # load_dag_from_ids(3,11,storage_path_scc+'comp.scf')
    # list_scc_dag = compute_links_scf(storage_path_scc,dict_scc_dag,storage_path_scc_dag)

    # p = pathlib.Path(storage_path_scc_dag)
    # for filename in list(p.glob('*_nodes.mspk')):
    #     n_wcc = re.search(r'\d+', os.path.basename(filename)).group()
    #     S = load_scc_dag(filename)
    #     S.compute_links(storage_path_scc+n_wcc+"/",n_wcc,storage_path_scc_dag)
    #     # print("  nodes :",S.c_nodes)
    #     # print("  links :",S.c_links)
    #     if 10000 > len(S.c_nodes) > 100:
    #         print("node path :", filename)
    #         # S.plot()
    #         # plt.show()
