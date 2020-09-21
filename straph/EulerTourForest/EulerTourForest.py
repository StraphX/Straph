import random
import numpy as np
import matplotlib.pyplot as plt
import msgpack,gc
from collections import defaultdict
import copy


import math,time
from straph import stream as sg
from straph import generators as gen

from straph.EulerTourForest.EulerTourTree import construct_euler_tour_tree, union_treap


# Source  : https://dl.acm.org/citation.cfm?id=320215


# TODO :? Utiliser les temps de départ des noeuds en lieu et place des priority des TreapNode ?
#  MIEUX : - Lors des stream graphs : prendre le liens qui partent le plus tard comme tree edges
#        - Lorsqu'on remplace un tree edge, prendre un non tree edge qui part le plus tard possible,
#          du coup les choisir dans l'ordre de départ inverse (héhé)


def spanning_forest_from_edge_list(edge_list):
    '''
    Construct a spanning forest from an edge list (graph may be disconnected)
    :param edge_list: [(u1,u2),(u0,u2),....]
    :return: Spanning Forest, which is a list of tuple : [(tree edge, non tree edge),...]
    '''
    a_l = defaultdict(set)
    for l in edge_list:
        u, v = l
        a_l[u].add(v)
        a_l[v].add(u)
    spanning_forest = []
    nodes = set(a_l.keys())
    while nodes:
        tree_edges = set()
        node = nodes.pop()
        visited, queue = set(), [node]
        prev = {}
        # DFS to get spanning tree (tree edges)
        while queue:
            node = queue.pop()
            if node not in visited:
                visited.add(node)
                nodes_to_visit = a_l[node] - visited
                queue.extend(nodes_to_visit)
                for n in nodes_to_visit:
                    prev[n] = node
                if node in prev:
                    tree_edges.add((prev[node], node))
        # To get non tree edges
        non_tree_edges = set()
        for u in visited:
            for v in a_l[u]:
                if (v, u) in tree_edges or (u, v) in tree_edges:
                    continue
                elif (v, u) in non_tree_edges:
                    continue
                else:
                    non_tree_edges.add((u, v))
        nodes -= visited
        spanning_forest.append((tree_edges, non_tree_edges))
    return spanning_forest


def construct_euler_tour_forest(edge_list):
    '''
    Construct an Euler Tour Forest from an edge list (the initial graph may be disconnected)
    :param edge_list:
    :return:
    '''
    SF = spanning_forest_from_edge_list(edge_list)
    tree_edge_2_node = {}
    Trees = []
    non_tree_edges_al = defaultdict(set)
    cnt_trees = 0
    for tree_edges, non_tree_edges in SF:
        print("Tree edges :", tree_edges)
        print("Non tree edges :", non_tree_edges)
        ETT, d = construct_euler_tour_tree(list(tree_edges))
        tree_edge_2_node.update(d)
        for l in non_tree_edges:
            u, v = l
            non_tree_edges_al[u].add(v)
            non_tree_edges_al[v].add(u)
        ETT.root.tree_number = cnt_trees
        Trees.append(ETT)
        cnt_trees +=1
    return EulerTourForest(trees=Trees, tree_edge_2_node=tree_edge_2_node,
                           non_tree_edges_al=non_tree_edges_al)


class ETF_collections(object):
    '''
    Collection of Spanning Euler Forest (1 for each level of the algorithm)
    '''

    def __init__(self):
        '''
        Forests : A list of Spanning Forests of different levels
        Edge_2_level : A dictionary associating an edge to his current level
        '''
        self.forests = []
        self.edge_2_level = {}

    def search(self, e):
        if e in self.edge_2_level:
            return True  #
        else:
            return False  # Rajouter un cas où les noeuds sont présents tous les deux MAIS dans des arbres différents

    def insert(self, e):
        '''
        Insert and edge in the Forest at the level 0
        :param e:
        :return:
        '''
        if not self.search(e):
            self.forests[0].insert_edge(e)  # We insert the edge in the level 0
            self.edge_2_level[e] = 0

    def remove(self, e):
        '''
        REMOVE AN EDGE in the current forest
        :param e:
        :return:
        '''
        if e in self.edge_2_level:
            self.replace(e, self.edge_2_level[e])

    def replace(self, e, level):
        '''
        REPLACE A TREE EDGE
        :param e: edge
        :param level: level of the edge to replace in the forest
        :return:
        '''
        while level > 0:
            status = self.forests[level].replace(e)
            if status:
                return
            level -= 1
            # return No replacement found, du coup ya du split dans l'air, done with .remove(e) ???
            # Move all edges of v_tree to the level i+1
            # Récuperer les non tree edges de v_tree et tester si il y en a qui reconnecte u_tree et v_tree
            #  Soit f un nontree dedges:
            #  - Si f ne connecte pas u_tree et v_tree, le move to the level i+1
            #  - Si f reconnecte u_tree and v_tree insert(f) and in u_tree


class EulerTourForest(object):
    '''
    Custom Data Structure Consisting in a Forest of Euler Tour Trees
    In the same manner as a Spanning Forest if composed of Spanning Trees
    '''
    def __init__(self, trees=[], tree_edge_2_node={}, non_tree_edges_al=defaultdict(set)):
        '''
        :param trees:  List of trees constituting the spanning forest (each Tree an EulerTourTree).
        :param tree_edge_2_tree: Dictionary associating a tree edge -> Euler Tour Tree it belongs
        :param tree_edge_2_node : Dictionary associating a tree edge -> the corresponding node (CTreapNode)
        :param non_tree_edges_al: Adjacency List (set) of non tree edges.
        '''
        self.trees = trees
        self.tree_edge_2_node = tree_edge_2_node
        self.non_tree_edges_al = non_tree_edges_al

    def __repr__(self):
        rep = " Tree edges : " + str([k for k in self.tree_edge_2_node.keys()]) + "\n"
        rep += " Non Tree edges : " + str(self.non_tree_edges_al) + "\n"
        for i in range(len(self.trees)):
            if self.trees[i]:
                rep += str(self.trees[i]) + "\n"
        return rep

    def is_tree_edge(self, e):
        '''
        Return true if *e* is tree edge, false otherwise
        :param e: an edge
        :return:
        '''
        if e in self.tree_edge_2_node:
            return e
        if (e[1], e[0]) in self.tree_edge_2_node:
            return (e[1], e[0])
        return False

    def plot(self, title=None):
        '''
        Plot the Euler Tour Trees in the Euler Tour Forest
        :param title: An optional title
        :return:
        '''
        for i, T in enumerate(self.trees):
            if T:
                if title:
                    T.plot(title)
                else:
                    T.plot(str(i) + " Tree of the Euler Tour Forest ")

    def check_edge_presence(self, e):
        '''
        Check if the edge is in the current spanning forest
        :param e:
        :return:
        '''
        if e[0] not in self.tree_edge_2_node:
            return False
        if e[1] not in self.tree_edge_2_node:
            return False
        return True

    def check_edge_endpoints(self, e):
        '''
        Check if one of the endpoints is present in th forest
        :param e:
        :return:
        '''
        if e[0] in self.tree_edge_2_node:
            return True
        if e[1] in self.tree_edge_2_node:
            return True
        return False

    def add_edge_to_tree(self, E, present_node, node_to_add):
        '''
        We assume that u is in E and we just add v
        :param E:
        :param v:
        :return:
        '''
        u_node = self.tree_edge_2_node[(present_node,present_node)][0]
        # Releaf
        # E.plot(" Inital E")
        if E.last != u_node:
            E = E.releaf(where=u_node)
        uv_node = E.insert(data=(present_node,node_to_add),inlast =True)
        # E.plot(" After insertion of "+str((present_node,node_to_add)))
        # print(" After insertion of "+str((present_node,node_to_add)))
        # print(E)
        v_node = E.insert(data = (node_to_add,node_to_add),inlast=True)
        # E.plot(" After insertion of "+str((node_to_add,node_to_add)))
        # print(" After insertion of "+str((node_to_add,node_to_add)))
        # print(E)
        vu_node = E.insert(data = (node_to_add,present_node),inlast = True)
        # E.plot(" After insertion of "+str((node_to_add, present_node)))
        # print(" After insertion of "+str((node_to_add, present_node)))
        # print(E)
        self.tree_edge_2_node[(present_node,node_to_add)] = [uv_node,vu_node]
        self.tree_edge_2_node[(node_to_add,node_to_add)] = [v_node]
        return E

    def insert_edge(self, l, SCC=None, format="cluster", exception=False, streaming_output=None):
        '''
        Insert an edge in the Euler Tour Forest
        :param e:
        :return:
        '''
        _, t0, t1, u, v = l
        e = (u,v)
        if (u,u) not in self.tree_edge_2_node and (v,v) not in self.tree_edge_2_node:
            # New Nodes, New edge, New Spanning Tree, New EulerTourTree
            E,d = construct_euler_tour_tree([(u,v)])
            self.tree_edge_2_node.update(d)
            self.trees.append(E)
            E.begin_time = t0
            E.root.tree_number = len(self.trees)-1
            # print(" Construct New Tree :",E.root.tree_number," with :",e)

        elif (u,u) not in self.tree_edge_2_node:
            # New Node, New edge
            v_tree_number = self.tree_edge_2_node[(v,v)][0].find_root().tree_number
            v_tree = self.trees[v_tree_number]
            v_tree.end_time = t0
            # self.write_to_msgpack(v_tree)
            self.add_to_scc(v_tree, SCC, format=format, exception=exception, streaming_output=streaming_output)


            E = self.add_edge_to_tree(v_tree,present_node=v,node_to_add=u)
            E.root.tree_number = v_tree_number
            E.begin_time = t0
            self.trees[v_tree_number] = E
            # print(" Add new node :",u," in tree : ",v_tree_number)

        elif (v,v) not in self.tree_edge_2_node:
            # New Node, New edge
            u_tree_number = self.tree_edge_2_node[(u, u)][0].find_root().tree_number
            u_tree = self.trees[u_tree_number]
            u_tree.end_time = t0
            # self.write_to_msgpack(u_tree)
            self.add_to_scc(u_tree, SCC, format=format, exception=exception, streaming_output=streaming_output)

            E = self.add_edge_to_tree(u_tree,present_node=u,node_to_add=v)
            E.root.tree_number = u_tree_number
            E.begin_time = t0
            self.trees[u_tree_number] = E
            # print(" Add new node :",v," in tree : ",u_tree_number)

        else:
            # Test if Merge else Update
            u_tree_number = self.tree_edge_2_node[(u,u)][0].find_root().tree_number
            v_tree_number = self.tree_edge_2_node[(v,v)][0].find_root().tree_number
            if u_tree_number == v_tree_number:  # They are already in the same tree
                # print(" Insert in tree number :", u_tree_number)
                if self.is_tree_edge(e):
                    return  # Nothing to do be do be do
                # Else insert in non tree edges, same cost as checking if its already a non tree                edge
                self.non_tree_edges_al[e[0]].add(e[1])
                self.non_tree_edges_al[e[1]].add(e[0])
            else:  # They are in different trees
                # print(" Merge Trees ", v_tree_number, " and ", u_tree_number)
                u_tree = self.trees[u_tree_number]
                v_tree = self.trees[v_tree_number]
                u_tree.end_time = t0
                v_tree.end_time = t0
                # self.write_to_msgpack(u_tree)
                # self.write_to_msgpack(v_tree)
                self.add_to_scc(u_tree, SCC, format=format, exception=exception, streaming_output=streaming_output)
                self.add_to_scc(v_tree, SCC, format=format, exception=exception, streaming_output=streaming_output)

                uv_tree = self.link_euler_tour_trees(u_tree, v_tree, e)
                uv_tree.begin_time = t0
                uv_tree.root.tree_number = u_tree_number
                # assert uv_tree.root.tree_number == self.tree_edge_2_node[(u,u)][0].find_root().tree_number
                # assert uv_tree.root.tree_number == self.tree_edge_2_node[(v, v)][0].find_root().tree_number
                self.trees[u_tree_number] = uv_tree
                self.trees[v_tree_number] = None
        return

    def replacement_edge(self, E1, E2):
        # We assume that E1 is smaller than E2 (TODO : implement a size of the tree (aka len(E1.nt_a_l))?
        r2 = E2.root
        for n in E1:
            u, v = n.data
            if u == v:
                for v in self.non_tree_edges_al[u]:
                    if self.tree_edge_2_node[(v, v)][0].find_root() == r2:
                        return (u, v)

        # for u in self.non_tree_edges_al:
        #     if self.tree_edge_2_node[(u, u)][0].find_root() == r1:
        #         for v in self.non_tree_edges_al[u]:
        #             if self.tree_edge_2_node[(v, v)][0].find_root() == r2:
        #                 return (u, v)
        return False

    def replace_edge(self, E1, E2):
        '''
        Find is there is a non tree edge linking E1 and E2
        :param E1: An euler tour tree
        :param E2: An euler tour tree
        :return: a replacement edge if found, false otherwise
        '''

        e = self.replacement_edge(E1, E2)
        if e:
            # print("  Replacement edge :", e)
            # print("  Found Replacement Edge :) hamdoulilah")
            E = self.link_euler_tour_trees(E1, E2, e)
            self.non_tree_edges_al[e[0]].remove(e[1])
            self.non_tree_edges_al[e[1]].remove(e[0])
            return E
        # else:
        #     print("  Did not Find Replacement Edge :( starfullah")
        return False

    def remove_edge(self, l, SCC=None, format="cluster", exception=False, streaming_output=None):
        '''
        Remove an edge from the Euler Tour Forest
        :param e:
        :return:
        '''

        _, t1, u, v = l
        e = (u,v)
        if (u,u) not in self.tree_edge_2_node:
            raise KeyError(" Trying to remove the link "+str(e)+ " whereas the node "+str(u)+" isn't even present !")
        if (v,v) not in self.tree_edge_2_node:
            raise KeyError(" Trying to remove the link "+str(e)+ " whereas the node "+str(v)+" isn't even present !")
        tree_number = self.tree_edge_2_node[(u, u)][0].find_root().tree_number
        current_tree = self.trees[tree_number]
        # print(" Tree Number :",tree_number)
        # print(" Remove in tree rooted at : ", current_tree.root.data)
        e = self.is_tree_edge(e)
        if e:
            copy_tree = copy.copy(current_tree)
            # print(" Remove tree edge : ", e)
            nodes = self.tree_edge_2_node[e]
            # Cut the euler tour into two distincts euler tour corresponding
            # to the removal of *e*
            E1, E2 = current_tree.cut(nodes)
            # Try to find a replacement edge among the non tree edges
            if E1 and E2:
                E = self.replace_edge(E1, E2)
            if E:
                self.trees[tree_number] = E
                E.begin_time = copy_tree.begin_time
                E.root.tree_number = tree_number
            else:
                copy_tree.end_time = t1
                # self.write_to_msgpack(copy_tree)
                self.add_to_scc(copy_tree, SCC, format=format, exception=exception, streaming_output=streaming_output)

                self.trees[tree_number] = E1
                E1.root.tree_number = tree_number
                E1.begin_time = t1
                l = len(self.trees)
                self.trees.append(E2)
                E2.root.tree_number = l
                E2.begin_time = t1
                # print(" Separate current tree into tree ",tree_number," and ",l)
                # print("E1 :\n",E1)
                # print("E2 :\n",E2)
            del self.tree_edge_2_node[e]
        else:
            # print(" Remove non tree edge : ", (u,v))
            self.non_tree_edges_al[u].remove(v)
            self.non_tree_edges_al[v].remove(u)

    def link_euler_tour_trees(self, E1, E2, e):
        '''
        Merge two euler tour tree
        :param E1: An euler tour tree
        :param E2: An euler tour tree
        :param e: An edge
        :return:
        '''
        u, v = e  # We know that u is in E1 and v in E2
        u_node = self.tree_edge_2_node[(u, u)][0]
        v_node = self.tree_edge_2_node[(v, v)][0]
        # print("   u pos :", u_node.data, "  v pos :", v_node.data)
        # Releaf
        if E1.last != u_node:
            E1 = E1.releaf(where=u_node)
            # print("  After releafing :")
            # print(E1)
            # E1.plot("  E1 after releafing in :" + repr(u_node.data))
        # Reroot (so releaf in v_node.pred)

        if E2.first != v_node:
            E2 = E2.releaf(where=v_node.pred)
            # print("  After rerooting :")
            # print(E2)
            # E2.plot("  E2 after rerooting in :" + repr(v_node.data))

        # print("###########################################################")
        # print("E1 first :", E1.first.data)
        # print("E1 last :", E1.last.data)
        uv_node = E1.insert(data=e, inlast=True)  # (u,u) is in E1

        # print("After insertion of :", e, " with data :", uv_node.data)
        # E1.plot("E1 after insertion of :" + repr(e))
        # print(E1)
        E = union_treap(E1, E2)
        # print(" After Union :")
        # E.plot("Union of E1 and E2 ")
        # print(E)
        vu_node = E.insert(data=(v, u), inlast=True)  # (v,v) is in E2

        self.tree_edge_2_node[e] = [uv_node]
        self.tree_edge_2_node[e].append(vu_node)
        # print(" After final insertion of :", (v, u), " with data :", vu_node.data)
        # E.plot(" After Final insertion of : " + repr((v, u)))
        # print(E)
        return E

    def add_to_scc(self, T, SCC, format="streaming", exception=False, streaming_output=None):
        '''
        Add the current tree *T* (aka a connected component) to *SCC*
        :param T: An Euler Tour Tree
        :param SCC: list of Strongly Connected Components
        :return: Updated SCC
        '''
        t0, t1 = T.begin_time, T.end_time
        if t0 != t1 or exception:
            L = T.get_data_in_priority_order()
            # nodes_set = set()
            # for l in L:
            #     u, v = l
            #     nodes_set.add(u[2])
            #     nodes_set.add(v[2])
            nodes_set = set([u for l in L for u in l])
            if len(nodes_set) > 1:
                if format == "streaming":
                    if streaming_output:
                        streaming_output.write(str(t0) + ";" + str(t1) + ";" + str(len(nodes_set)))
                        streaming_output.write("\n")
                    else:
                        c = (t0, t1, len(nodes_set))
                        SCC.append(c)
                if format == "cluster":
                    current_comp = [(t0, t1, n) for n in nodes_set]
                    SCC.append(current_comp)

    def write_to_msgpack(self,T):
        '''
        A handler to store a connected component
        :param storage_file:
        :return:
        '''
        t0, t1 = T.begin_time, T.end_time
        if t0 != t1:
            L = T.get_data_in_priority_order()
            links = set()
            # Add tree edges
            for l in L:
                u,v = l
                if l not in links and (v, u) not in links and u != v:
                    links.add(l)
                    # Add Non Tree Edges
                    for n in self.non_tree_edges_al[u]:
                        if (n,u) not in links and (u,n) not in links:
                            links.add((n,u))
                    for n in self.non_tree_edges_al[v]:
                        if (n,v) not in links and (v,n) not in links:
                            links.add((n,v))
            if links:
                # print("Links to store :", tuple(links))
                storage_file.write(packer.pack((t0,t1,tuple(links))))
            # TODO : compteur id comp
        return


# def dynamic_connectivity(E, M):
#     ETF = construct_euler_tour_forest(E)
#     print("Initial Euler Tour Forest :\n", ETF)
#     ETF.plot("Initial ETT")
#     while M:
#         c, u, v = M.pop()
#         if c == -1:  # Deletion
#             print("\nDeletion : ", (u, v))
#             ETF.remove_edge((u, v))
#             # plt.show()
#         if c == 1:  # Insertion
#             print("\nInsertion : ", (u, v))
#             ETF.insert_edge((u, v))
#             # plt.show()
#         print("ETF :\n", ETF)
#     print("ETF after sequence :\n", ETF)
#     # for T in ETF.trees:
#     #     if T:
#     #         ETF.write_to_msgpack(T)
#     return


# def scc_etf(input_file):
#     ##############
#     gc.disable()
#     ##############
#     ETF = construct_euler_tour_forest([])
#
#     with open(input_file,'rb') as input_file:
#         unpacker = msgpack.Unpacker(input_file,use_list=False)
#         for l in unpacker:
#             # print(l)
#             if l[0] == -1:
#                 # print("Deletion")
#                 ETF.remove_edge(l[1:])
#             if l[0] == 1:
#                 # print("Insertion")
#                 ETF.insert_edge(l[1:])
#             # print("ETF :\n",ETF)
#     #############
#     gc.enable()
#     #############
#     return


def strongly_connected_components_ETF(S, isolated_nodes=True, format="cluster", streaming_output=None,
                                      free_memory=False):
    if streaming_output:
        opt = open(streaming_output, 'w')
    else:
        opt = None

    SCC = []
    ETF = construct_euler_tour_forest([])
    # E = S.augmented_ordered_links()
    E = S.ordered_batch_links(free_memory=free_memory)
    for batch in E:
        c = batch[0][0]
        if c == 1:
            ETF.insert_edge(batch[0], SCC, format=format, exception=True, streaming_output=opt)
            for l in batch[1:]:
                ETF.insert_edge(l, SCC, format=format, streaming_output=opt)
        else:
            # for l in b:
            #     print("l :",l)
            ETF.remove_edge(batch[0], SCC, format=format, exception=True, streaming_output=opt)
            for l in batch[1:]:
                ETF.remove_edge(l, SCC, format=format, streaming_output=opt)

    if isolated_nodes:
        for c in S.get_isolated_nodes():
            if format == "streaming":
                if opt:
                    opt.write(str(c[0]) + ";" + str(c[1]) + ";" + str(1))
                    opt.write("\n")
                else:
                    c = (c[0], c[1], 1)
                    SCC.append(c)

            if format == "cluster":
                SCC.append([c])

    return SCC



import time
if __name__ == '__main__':
    random.seed(13)
    np.random.seed(13)
    T = [0, 1000]
    nb_node = 250

    occurrence_law = 'poisson'
    presence_law = 'poisson'

    occurrence_param_node = 2
    presence_param_node = 500
    occurrence_param_link = 20
    presence_param_link = 50
    #
    p_link = 1.3 * np.sqrt(nb_node) / nb_node
    # print("p_link :", p_link)  # " with eps :",eps)
    S = gen.erdos_renyi(T,
                        nb_node,
                        occurrence_law,
                        occurrence_param_node,
                        presence_law,
                        presence_param_node,
                        occurrence_law,
                        occurrence_param_link,
                        presence_law,
                        presence_param_link,
                        p_link)
    S.describe()

    # __directory__ = "/home/leo/Dev/CODE/Straph/examples/"
    # __file__ = "example"
    # S = sg.read_stream_graph(path_nodes=__directory__ + __file__ + "_nodes.sg",
    #                       path_links=__directory__ + __file__ + "_links.sg")


    chrono = time.time()
    scc_direct = S.strongly_connected_components(method="Direct", format="cluster", isolated_nodes=True)
    print("SCC Direct DONE :", time.time() - chrono)
    #
    #
    chrono = time.time()
    scc_etf = S.strongly_connected_components(method="FD", format="cluster", isolated_nodes=True)
    print("SCC ETF DONE in :",time.time()-chrono)
    #
    #
    chrono = time.time()
    scc_uf = S.strongly_connected_components(method="UF", format="cluster")
    print("SCC UF DONE in :", time.time() - chrono)

    # S.plot(clusters=scc_direct,title="SCC Classic")
    # S.plot(clusters=scc_etf, title="SCC UF")
    # S.plot()
    # plt.show()

    set_etf = set()
    set_classic = set()
    set_uf = set()

    for i in scc_uf:
        set_nodes = set()
        t0, t1 = None, None
        for (t0, t1, n) in i:
            set_nodes.add(n)
        set_uf.add((t0, t1, tuple(sorted(set_nodes))))


    for i in scc_etf:
        set_nodes = set()
        t0,t1 = None,None
        for (t0,t1,n) in i:
            set_nodes.add(n)
        set_etf.add((t0, t1, tuple(sorted(set_nodes))))

    for i in scc_direct:
        set_nodes = set()
        t0,t1 = None,None
        for (t0,t1,n) in i:
            set_nodes.add(n)
        set_classic.add((t0,t1,tuple(sorted(set_nodes))))

    print("N comp classic:", len(set_classic))
    print("N comp etf :", len(set_etf))
    print("N comp UF :", len(set_uf))

    for i in set_etf:
        if i not in set_classic:
            print("etf not classic :", i)

    for i in set_classic:
        if i not in set_etf:
            print("classic not etf :",i)

    print("\n" + "######################" + "\n")

    for i in set_uf:
        if i not in set_classic:
            print("UF not classic :", i)

    for i in set_classic:
        if i not in set_uf:
            print("classic not UF :", i)

    print("\n" + "######################" + "\n")

    for i in set_uf:
        if i not in set_etf:
            print("uf not etf :", i)

    for i in set_etf:
        if i not in set_uf:
            print("etf not uf :", i)

    # for j in set_uf:
    #     if j not in set_classic:
    #         print("UF not classic :", j)
    #     if j not in set_etf:
    #         print("UF not ETF :", j)
