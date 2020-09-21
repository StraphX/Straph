import copy
from joblib import Parallel, delayed
from collections import defaultdict
import matplotlib.pyplot as plt
import random
import math
import msgpack
from straph.utils import get_cmap

def algo_kcores_batagelj(a_l, degrees, core_ordering=False):
    '''
    Compute k_cores of a static graph from its adjacency list and nodes degrees

    :param a_l: Adjacency list of the graph
    :param degrees: degree of each node
    :param core_ordering: Optionnal parameter, if True return a coreness ordering
    :return: the core number of each node in the graph
    References
    ----------
    [1] An O(m) Algorithm for Cores Decomposition of Networks
    Vladimir Batagelj and Matjaz Zaversnik, 2003.
    http://arxiv.org/abs/cs.DS/0310049
    '''
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
        u = (t0, t1,l[0])
        v = (t0, t1,l[1])
        nodes.add(u)
        nodes.add(v)
        if ordering[v] < ordering[u]:
            a_l[v].add(u)
        elif ordering[u] < ordering[v]:
            a_l[u].add(v)
    return a_l

def algo_kcliques_KCList(k, a_l, node_label, C=[], R=[]):
        '''
        Compute k_cliques of a static graph from its core odering dag.
        :param k: The parameter k, number of nodes in the considered cliques
        :param a_l: Adjacency list of the core ordering dag
        :param node_label: label of each node
        :param C: List of current cliques
        :param R: List of completed cliques
        :param DEBUG:
        :return: The k_cliques of the graph.
        References
        ---------
        [2] Listing k-cliques in Sparse Real-World Graphs
        Maximilien Danisch, et al., 2018
        https://dl.acm.org/citation.cfm?id=3186125
        '''
        if k == 2:
            for u in node_label:
                if a_l[u]:
                    print(" u : ", u)
                    print(" current C :", C)
                    for v in a_l[u]:
                        if v in node_label:
                            if node_label[v] == k:
                                print("\t v : ",v)
                                final_C = copy.copy(C) + [u, v]
                                print("\t final C:",final_C,"\n")
                                R.append(final_C)
        else:
            for u in node_label:
                # INDUCED NODES
                new_node_label = defaultdict(int)
                if a_l[u]:
                    print(" u :", u)
                    print(" node label :", node_label)
                    for v in a_l[u]:
                        print("\t v :",v)

                        if v in node_label:
                            if node_label[v] == k:
                                new_node_label[v] = k - 1
                    # RECURSION
                    print("\t new node label :",new_node_label)

                    new_C = copy.copy(C) + [u]
                    print("\t new C :",new_C)

                    if new_node_label:
                        algo_kcliques_KCList(k - 1, a_l, new_node_label, C=new_C, R=R)
        return R


def neighborhood_and_degrees_from_links(t0, t1, links):
    a_l = defaultdict(set)
    for l in links:
        u = (t0, t1,l[0])
        v = (t0, t1,l[1])
        a_l[u].add(v)
        a_l[v].add(u)
    degrees = {n: len(a_l[n]) for n in a_l}
    return a_l, degrees


class connected_component:
    def __init__(self,
                 id=None,
                 times=None,
                 nodes=None,
                 links=None,
                 link_presence=None
                 ):
        '''
        A basic constructor for a connected component object

        :param id : id in SCC STREAM (a string)
        :param times = [beginning time, ending time]
        :param nodes : A set of nodes present in the component
        :param links : A set of links present in the component (Only useful during construction)
        :param link_presence : a default dictionary (with list as default) link to time of presence
        '''
        self.id = id
        self.times = times
        self.nodes = nodes
        self.links = links
        self.link_presence = link_presence

    def __copy__(self):
        t = copy.copy(self.times)
        n = copy.copy(self.nodes)
        lp = defaultdict(list, {k: copy.copy(v) for k, v in self.link_presence.items()})
        return connected_component(times=t,
                                   nodes=n,
                                   link_presence=lp)

    # def __del__(self):
    #     del self.id
    #     del self.times
    #     del self.nodes
    #     del self.links
    #     del self.link_presence


    def set_id(self,id):
        self.id=id

    def dump_link_presence(self,output):
        '''
        Dump the connected component link_presence to a msgpack_file
        :param output: Writable bytes file open(...,'wb')
        :return:
        '''
        with open(output+self.id+'.pck','wb') as output:
            msgpack.dump(self.link_presence,output)
        del self.link_presence
        # del self.links
        # self.link_presence.clear()
        # gc.collect()

    def load_link_presence(self,input):
        '''
        Retrieve the connected component link_presence from a msgpack file
        :param input: Readable bytes file open(...,'rb')
        :return:
        '''
        if self.id:
            with open(input+self.id+'.pck','rb') as input:
                self.link_presence = msgpack.load(input,use_list=False)


    def clear(self):
        self.nodes.clear()
        self.links.clear()
        self.times.clear()
        self.link_presence.clear()

    def surface(self):
        return sum([lt1 - lt0 for v in self.link_presence.values()
                    for lt1, lt0 in zip(v[::2], v[1::2])])

    def size(self):
        return len(self.nodes)

    def set_begin_time(self, t):
        self.times[0] = t
        to_del = []
        for k, v in self.link_presence.items():
            if v[-1] <= t:
                to_del.append(k)
            else:
                v[0] = t
        for l in to_del:
            del self.link_presence[l]

    def set_end_time(self, t):
        self.times[1] = t
        for v in self.link_presence.values():
            if v[-1] > t:
                v[-1] = t

    def add_link(self, link):
        u, v = link[0]
        if u not in self.nodes:
            self.nodes.add(u)
        if v not in self.nodes:
            self.nodes.add(v)
        self.links.add(link[0])
        self.link_presence[link[0]] += [link[1], link[2]]

    def merge(self, C):
        self.nodes |= C.nodes
        self.links |= C.links
        self.link_presence.update(C.link_presence)

    def remove_link(self, link):
        self.links.discard(link)

    def remove_node(self, node):
        self.nodes.discard(node)

    def split_on_link(self, link):
        self.links.discard(link)
        R = self.BFS_split(link)
        if R:
            return R
        else:
            return False

    def BFS_split(self, link):
        # CUSTOM BFS
        # In a component which is being splitted, the output can
        # only be 1 or 2 components.
        # If at any moment u and v find themselves in the same component <- stop
        u, v = link
        comp_bfs = []
        link_comp_bfs = []
        dict_node_2_component_bfs = {}
        for l in self.links:
            n1, n2 = l
            # logger.debug("\t current link : "+str(l))
            if n1 not in dict_node_2_component_bfs and n2 not in dict_node_2_component_bfs:
                # logger.debug("\t Create comp")
                dict_node_2_component_bfs[n1] = len(comp_bfs)
                dict_node_2_component_bfs[n2] = len(comp_bfs)
                comp_bfs.append(set([n1, n2]))
                link_comp_bfs.append(set([l]))

            elif n1 in dict_node_2_component_bfs and n2 not in dict_node_2_component_bfs:
                # logger.debug("\t ADD "+str(n2)+" in "+str(n1)+" component")
                nc1 = dict_node_2_component_bfs[n1]
                comp_bfs[nc1].add(n2)
                link_comp_bfs[nc1].add(l)
                dict_node_2_component_bfs[n2] = nc1

            elif n1 not in dict_node_2_component_bfs and n2 in dict_node_2_component_bfs:
                # logger.debug("\t ADD "+str(n1)+" in "+str(n2)+" component")
                nc2 = dict_node_2_component_bfs[n2]
                comp_bfs[nc2].add(n1)
                link_comp_bfs[nc2].add(l)
                dict_node_2_component_bfs[n1] = nc2

            elif dict_node_2_component_bfs[n1] != dict_node_2_component_bfs[n2]:
                # logger.debug("\t Merge "+str(n2)+" component into "+str(n1)+" component")
                nc1, nc2 = dict_node_2_component_bfs[n1], dict_node_2_component_bfs[n2]
                comp_bfs[nc1] |= comp_bfs[nc2]
                link_comp_bfs[nc1] |= link_comp_bfs[nc2]
                link_comp_bfs[nc1].add(l)
                link_comp_bfs[nc2].clear()
                for node in comp_bfs[nc2]:
                    dict_node_2_component_bfs[node] = nc1
                comp_bfs[nc2].clear()
            else:
                link_comp_bfs[dict_node_2_component_bfs[n1]].add(l)
            # Test if we can stop the bfs_with_path
            if u in dict_node_2_component_bfs and v in dict_node_2_component_bfs and \
                    dict_node_2_component_bfs[u] == dict_node_2_component_bfs[v]:
                # logger.debug("\t Test succesfull : " + str(u) + "is with " + str(v))
                return False
        # logger.debug("\t comp bfs_with_path : " + str(comp_bfs))
        # logger.debug("\t link bfs_with_path : " + str(link_comp_bfs))
        C = []
        for c_nodes, c_links in zip(comp_bfs, link_comp_bfs):
            C.append(
                connected_component(times=copy.copy(self.times),
                                    nodes=c_nodes,
                                    links=c_links,
                                    link_presence=defaultdict(list, {l: self.link_presence[l] for l in c_links if
                                                                     l in self.link_presence}
                                                              )))
        return C

    def get_kcores(self,storage_path=None):
        # Get interaction times
        # Divides into slices
        # Compute cores with Batagelj (can be parallelized)
        if storage_path:
            self.load_link_presence(storage_path)
        L = defaultdict(list)
        if self.size() == 2:
            L[1] = [( self.times[0], self.times[1],list(self.nodes)[0])]
            L[1].append(( self.times[0], self.times[1],list(self.nodes)[1]))
            return L
        interact_times = sorted(set([t for v in self.link_presence.values()
                                     for t in v]))
        time_2_pos = {t: i for t, i in zip(interact_times, range(len(interact_times)))}
        inter_links = [[] for k in range(len(interact_times) - 1)]
        for l in self.link_presence:
            for lp0, lp1 in zip(self.link_presence[l][::2], self.link_presence[l][1::2]):
                for i in range(time_2_pos[lp0], time_2_pos[lp1]):
                    inter_links[i].append((l[0], l[1]))

        def para_core(i):
            d = defaultdict(list)
            t0, t1 = interact_times[i], interact_times[i + 1]
            current_links = inter_links[time_2_pos[t0]]
            a_l, degrees = neighborhood_and_degrees_from_links(t0, t1, current_links)
            cores = algo_kcores_batagelj(a_l, degrees)  # Algo Batelgej
            for v in cores:
                d[cores[v]].append(v)
            return d

        r = Parallel(n_jobs=1)(delayed(para_core)(i) for i in range(len(interact_times) - 1))
        for l in r:
            for k, v in l.items():
                L[k] += v
        if storage_path:
            self.link_presence.clear()
        return L

    def get_stable_parts(self):
        interact_times = sorted(set([t for v in self.link_presence.values()
                                     for t in v]))
        time_2_pos = {t: i for t, i in zip(interact_times, range(len(interact_times)))}
        inter_links = [set() for k in range(len(interact_times) - 1)]
        for l in self.link_presence:
            for lp0, lp1 in zip(self.link_presence[l][::2], self.link_presence[l][1::2]):
                for i in range(time_2_pos[lp0], time_2_pos[lp1]):
                    inter_links[i].add(l[0])
                    inter_links[i].add(l[1])
        L = []
        for j in range(len(interact_times)-1):
            L.append([(interact_times[j],interact_times[j+1],u) for u in inter_links[time_2_pos[interact_times[j]]]])
        return L

    def get_stable_components(self):
        '''
        :return: stable components stemmed from the current component
        '''
        interact_times = sorted(set([t for v in self.link_presence.values()
                                     for t in v]))
        time_2_pos = {t: i for t, i in zip(interact_times, range(len(interact_times)))}
        inter_nodes = [set() for k in range(len(interact_times) - 1)]
        inter_links = [[] for k in range(len(interact_times) - 1)]
        for l in self.link_presence:
            for lp0, lp1 in zip(self.link_presence[l][::2], self.link_presence[l][1::2]):
                for i in range(time_2_pos[lp0], time_2_pos[lp1]):
                    inter_nodes[i].add(l[0])
                    inter_nodes[i].add(l[1])
                    inter_links[i].append((l[0],l[1]))
        stable_components = []
        for j in range(len(interact_times)-1):
            c = connected_component(id=j,
                                    times=(interact_times[j],interact_times[j+1]),
                                    nodes = set([u for u in inter_nodes[time_2_pos[interact_times[j]]]]),
                                    links = [l for l in inter_links[time_2_pos[interact_times[j]]]]
                                    )
            stable_components.append(c)
        return stable_components


    def get_kcliques(self,storage_path=None):
        # Get interaction times
        # Divides into slices
        # Compute cores with KCList (can be parallelized)
        if storage_path:
            self.load_link_presence(storage_path)
        L = defaultdict(list)
        if self.size() == 2:
            # L[1] = [(list(self.nodes)[0],self.times[0], self.times[1] )]
            # L[1].append((list(self.nodes)[1],self.times[0], self.times[1]))
            return L
        interact_times = sorted(set([t for v in self.link_presence.values()
                                     for t in v]))
        time_2_pos = {t: i for t, i in zip(interact_times, range(len(interact_times)))}
        inter_links = [[] for k in range(len(interact_times) - 1)]
        for l in self.link_presence:
            for lp0, lp1 in zip(self.link_presence[l][::2], self.link_presence[l][1::2]):
                for i in range(time_2_pos[lp0], time_2_pos[lp1]):
                    inter_links[i].append((l[0], l[1]))

        def para_clique(i):
            cliques = {}
            t0, t1 = interact_times[i], interact_times[i + 1]
            current_links = inter_links[time_2_pos[t0]]
            a_l, degrees = neighborhood_and_degrees_from_links(t0, t1, current_links)
            cores, core_ordering = algo_kcores_batagelj(a_l, degrees, core_ordering=True)
            max_core_number = max(cores.values())
            #print("\n t0,t1 :",t0,t1)
            #print("Max core number :",max_core_number)
            #print("Core ordering :",core_ordering)
            K = 3
            a_l = get_graph_from_ordering(t0, t1, current_links, core_ordering)
            #print(" a_l :",a_l)
            #print("degrees :",degrees)
            while K <= max_core_number + 1 :
                node_label = defaultdict(int, {n: K for n in degrees})
                cliques[K] = algo_kcliques_KCList(K, a_l, node_label, R=[])
                K += 1
            #print("kcliques :",cliques)
            return cliques

        r = Parallel(n_jobs=1)(delayed(para_clique)(i) for i in range(len(interact_times) - 1))
        for l in r:
            for k, v in l.items():
                L[k] += v
        if storage_path:
            self.link_presence.clear()
        return L

    def adjacency_list_at_t(self, t):
        a_l = defaultdict(set)
        for l, lp in self.link_presence.items():
            for t0, t1 in zip(lp[::2], lp[1::2]):
                if t0 <= t <= t1:
                    a_l[l[0]].add(l[1])
                    a_l[l[1]].add(l[0])
        return a_l

    def adjacency_list(self):
        a_l = defaultdict(list)
        for l, lp in self.link_presence.items():
            for t0, t1 in zip(lp[::2], lp[1::2]):
                a_l[l[0]].append((t0, t1, l[1]))
                a_l[l[1]].append((t0, t1, l[0]))
        for k in a_l.keys():
            a_l[k].sort()
        return a_l

    def random_path(self, source, destination, rand_time=False):
        '''
        Return a random path between the source and the destination inside the component
        :param source:
        :param destination:
        :return:
        '''
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

    def shortest_path(self, source, destination):
        '''
        Return the shortest path between the source and the destination inside the component
        :param source:
        :param destination:
        :return:
        '''
        return

    def random_path_ss(self, source):
        '''
        ss : single source
        Return random paths between the source and all the other nodes in the component
        :param source:
        :return:
        '''
        t = random.random() * (self.times[1] - self.times[0]) + self.times[0]  # Random time
        print("Choosen time :", t)
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


    def shortest_path_ss(self, source,storage_path = None):
        '''
        weight\forall link == 1
        :param source:
        :return:
        '''
        if storage_path:
            self.load_link_presence(storage_path)
        S = set() # Visited Nodes
        L = set() # Visited Links
        D = {n : math.inf for n in self.nodes}
        P = {n: None for n in self.nodes}           # Previous node
        T = {n: self.times[0] for n in self.nodes}  # Time of arrival from previous node
        a_l = self.adjacency_list()                 #(t0,t1,node)
        D[source] = 0
        while S != self.nodes:
            v = min(set(D.keys()) - S, key=D.get)
            for (t0,t1,w) in set(a_l[v])-L:
                if t1 < T[w]:   # Si le lien est déjà fini, on ne peut plus le prendre (sans déconner)!
                    continue
                new_path = D[v] +1
                if new_path< D[w]:
                    D[w] = new_path
                    P[w] = v
                    T[w] = t0
                elif new_path == D[w] and t0 < T[w]:
                    P[w] = v
                    T[w] = t0
            S.add(v)
            L |= set([el for el in a_l[v]])
        print("D :", D)
        print("P :", P)
        print("T :", T)
        paths = {}
        targets = copy.copy(self.nodes)
        targets.discard(source)
        while targets:
            v = targets.pop()
            destination = v
            path = []
            times = []
            print("v : ",v)
            while v is not None:
                path.append(v)
                times.append(T[v])
                v = P[v]
            times = list(reversed(times))
            path = list(reversed(path))
            for i in range(len(times)-1):
                if times[i+1] < times[i]:
                    times[i+1] = times[i]
            print("path :", path)
            paths[destination] = (times,path)
        print("paths :", paths)
        if storage_path:
            self.link_presence.clear()
        return paths


    def random_path_pw(self):
        '''
        pw : pairwise
        Return random paths between every node (pairwise)
        :return:
        '''
        paths = []
        for n in self.nodes:
            paths.append(self.random_path_ss(n))
        return paths

    def shortest_path_pw(self):
        '''
        pw : pairwise
        Return all the shortest path between every node (pairwise)
        :return:
        '''

        return

    def plot(self, links=True, link_pimping=False, title=None):
        # todo : Adapt plot function from Stream.
        '''
                Display in an elegant way a small stream graph
                We can also specify a path to save the animated of ploted stream graph
                :param animated: If we want an animation with cursor which move according to the time
                :return: A matplotlib plot
                '''
        lnodes = len(self.nodes)
        c_map = get_cmap(lnodes)
        dict_colors = {n: c_map(i) for n, i in zip(self.nodes, range(lnodes))}
        # random.shuffle(l_colors)
        fig = plt.figure()
        # Plot Clusters
        nodes_list = list(self.nodes)

        for p in self.nodes:
            coln = dict_colors[p]
            # plt.axhline(p, linestyle='--', linewidth=0.7,
            #             color=coln, alpha=0.3)
            plt.hlines([nodes_list.index(p)], xmin=self.times[0], linewidth=2,
                       xmax=self.times[1], colors=coln, alpha=1)
        # Plot Links
        if links:
            for k, lp in self.link_presence.items():
                id1 = nodes_list.index(k[0])
                id2 = nodes_list.index(k[1])
                coln1 = dict_colors[k[0]]
                coln2 = dict_colors[k[1]]
                idmax = max(id1, id2)
                idmin = min(id1, id2)
                eps = random.choice([1, -1]) * (random.random() / 5)
                if link_pimping:
                    plt.hlines([(idmax + idmin) / 2 + eps] * (len(lp) // 2), xmin=lp[::2], xmax=lp[1::2],
                               linestyles=(3, (3, 3)),
                               color=coln1,
                               linewidth=2,
                               alpha=0.6)
                    plt.hlines([(idmax + idmin) / 2 + eps] * (len(lp) // 2), xmin=lp[::2], xmax=lp[1::2],
                               linestyles=(0, (3, 3)),
                               color=coln2,
                               linewidth=2,
                               alpha=0.6)
                else:
                    plt.hlines([(idmax + idmin) / 2 + eps] * (len(lp) // 2), xmin=lp[::2], xmax=lp[1::2],
                               colors='k',
                               linewidth=1.7,
                               alpha=0.5)
                plt.vlines(lp[::2],
                           ymin=idmin, ymax=idmax,
                           linewidth=1.4, alpha=0.15)
                # ARRIVALS
                plt.plot([lp[::2]], [idmin], color='#004d00', marker='^', alpha=1, markersize=7)
                plt.plot([lp[::2]], [idmax], color='#004d00', marker='v', alpha=1, markersize=7)
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
        plt.tight_layout()
        return