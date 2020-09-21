import math
import time
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from heapq import *
from sortedcollections import SortedList,SortedSet
import sortedcontainers
from collections import defaultdict, deque


# TODO: - Add detailed paths + plots : !!!!!!! (lot of work)


##########################################
#       Visualisation for L-Algorithm    #
##########################################

def plot_adjacency_list(S, a_l):
    '''
    Plot the current adjacency list *a_l*.

    :param S: A stream graph (we get its labels)
    :param a_l: an adjacency list
    :return: Plot of adjacency list
    '''
    fig = plt.figure()
    ax = plt.axes()

    G = nx.Graph()
    set_edge = set()
    for k, v in a_l.items():
        for i in v:
            if (k, i[1]) not in set_edge and (i[1], k) not in set_edge:
                G.add_edge(S.node_to_label[k[2]], S.node_to_label[i[1][2]], t1=i[0])
                set_edge.add((k, i[1]))
    pos = nx.circular_layout(G)
    nx.draw_networkx_nodes(G, pos, node_size=700,
                           node_color='#5a5ff5', alpha=0.5, ax=ax)
    nx.draw_networkx_edges(G, pos, edge_color='#2d5986',
                           alpha=0.3, width=5, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=20, ax=ax)
    nx.draw_networkx_edge_labels(G, pos, font_size=20, ax=ax)
    plt.xlabel("t", fontname='Ubuntu', fontsize=20, color='#476b6b')
    for spine in plt.gca().spines.values():
        spine.set_visible(False)
    plt.tick_params(top=False, bottom=False, right=False, left=False, labelbottom=False, labelleft=False)
    plt.tight_layout()
    plt.show()


def plot_stream(S, t):
    '''
    Plot the stream graph and a marker for the instant *t*.

    :param S: A stream Graph
    :param t:  An instant
    :return:  Plot of the stream graph
    '''
    S.plot()
    ax = plt.gca()
    y = np.linspace(0, len(S.nodes), 1000)
    x = [t] * 1000
    line, = ax.plot(x, y, color='#476b6b', lw=8, alpha=0.3)


def get_temporal_sources(S, source):
    '''
    Return the maximal segmented nodes (t0,t1,u) s.t u = source.

    :param S: Stream Graph
    :param source: source node
    :return:
    '''
    return [(t0, t1, source) for t0, t1 in zip(S.node_presence[source][::2], S.node_presence[source][1::2])]


#######################################
#           F-Algorithm               #
#######################################


def F_Algorithm(S, L_functions):
    """
    An implementation of F-Algorithm (add ref).
    Pairwise Algorithm to compute temporal paths in Stream Graph.

    :param S: A Stream Graph
    :param L_functions: Functions according to the path's type (supported :
        - foremost path
        - shortest foremost path
        - shortest path
        - fastest shortest path
        - fastest path
        - shortest fastest path)
    :return:
    """
    # Initialisation
    F = {}
    for n, np in zip(S.nodes, S.node_presence):
        for t0, t1 in zip(np[::2], np[1::2]):
            s = (t0, t1, n)
            F[s] = L_functions['initialisation'](s)
    temporal_adjacency_list = defaultdict(set)  # currently present links in the form of adjacency list
    nodes_to_udpate = set()  # Reached nodes so far
    E = S.augmented_ordered_links()
    for i in E:
        c = i[0]
        if c == 1:  # Arrival
            _, t0, t1, u, v = i
            l = (t0, t1, u, v)
            temporal_adjacency_list[u].add((t1, v))
            temporal_adjacency_list[v].add((t1, u))
            nodes_to_udpate.add(u)
            nodes_to_udpate.add(v)
            F_update(F, [l], temporal_adjacency_list, nodes_to_udpate, L_functions)
        else:
            _, t1, u, v = i
            temporal_adjacency_list[u].remove((t1, v))
            temporal_adjacency_list[v].remove((t1, u))
    return F


def F_update(F, batch_arrival, temporal_adjacency_list, nodes_to_update, L_functions):
    """
    Proceed to BFS_Update W_n for every encounter node so far (i.e. in nodes_to_update).
    :param F:
    :param batch_arrival:
    :param temporal_adjacency_list:
    :param nodes_to_update:
    :param L_functions:
    :return:
    """
    for n in nodes_to_update:
        bfs_update(F[n], n, batch_arrival, temporal_adjacency_list, L_functions)


#######################################
#           L-Algorithm               #
#######################################

def L_Algorithm(S, source,
                initialisation_function,
                update_function,
                source_init_function,
                start_time=None,
                is_temporal_source=False,
                sfp_special_case=False,
                verbose=False,
                E=None):
    """
    An implementation of L-Algorithm (described in ~add ref)
    Single Source Algorithm to compute temporal paths in Stream Graph

    :param S: A Stream Graph
    :param source: A temporal source node
    :param L_functions: Functions according to the path's type (supported :
            - foremost path
            - shortest foremost path
            - shortest path
            - fastest shortest path
            - fastest path
            - shortest fastest path)
    :return: R containing the best results
    # IMPORTANT : We suppose that we are in the WCC of 'source' otherwise it's fucking expensive !
    """

    L, R = initialisation_function(S,source)
    temporal_adjacency_list = defaultdict(set)  # key : temporal node, value : (temporal node, end time of link)

    if E is None:
        E = S.ordered_events()

    batch_arrival = []
    t_last_arrival = None

    if sfp_special_case:
        batch_departure = []
        t_last_departure = None

    # ADD Counter With %
    cnt_event = 0

    for i in E:
        c = i[0]
        #
        # if verbose:
        #     cnt_event += 1
        #     if cnt_event % 100000 == 0:
        #         print((cnt_event / len(E)) * 100, "% done")
        #
        if c == 1:  # Link ARRIVAL
            _, t0, t1, u, v = i

            if t1 >= start_time:
                l = (t0, t1, u, v)

                if sfp_special_case and batch_departure:
                    if temporal_adjacency_list[source[1]]:
                        bfs_update(L, R, source,
                                   batch_departure,
                                   temporal_adjacency_list,
                                   update_function,
                                   sfp_special_case=sfp_special_case)
                    for t1d,ud,vd in batch_departure:
                        temporal_adjacency_list[ud].remove((t1d, vd))
                        temporal_adjacency_list[vd].remove((t1d, ud))
                    batch_departure = []

                if t_last_arrival is None:
                    t_last_arrival = t0

                if t0 == t_last_arrival:
                    batch_arrival.append(l)

                else:
                    t_last_arrival = t0
                    if batch_arrival:
                        bfs_update(L, R, source, batch_arrival,
                                   temporal_adjacency_list,
                                   update_function,
                                   sfp_special_case=sfp_special_case)
                    batch_arrival = [l]

                temporal_adjacency_list[u].add((t1, v))
                temporal_adjacency_list[v].add((t1, u))


        elif c == -1:  # LINK DEPARTURE
            _, t1, u, v = i

            if t1 >= start_time:
                if batch_arrival:
                    bfs_update(L, R, source,
                               batch_arrival,
                               temporal_adjacency_list,
                               update_function,
                               sfp_special_case=sfp_special_case)
                    batch_arrival = []

                if sfp_special_case:
                    if t_last_departure is None:
                        t_last_departure = t1
                    if t1 == t_last_departure:
                        batch_departure.append((t1, u, v))
                    else:
                        t_last_departure = t1
                        if batch_departure:
                            if temporal_adjacency_list[source[1]]:
                                bfs_update(L, R, source, batch_departure,
                                           temporal_adjacency_list,
                                           update_function,
                                           sfp_special_case=sfp_special_case)
                            for t1d,ud,vd in batch_departure:
                                temporal_adjacency_list[ud].remove((t1d, vd))
                                temporal_adjacency_list[vd].remove((t1d, ud))
                        batch_departure = [(t1, u, v)]

                if not sfp_special_case:
                    temporal_adjacency_list[u].remove((t1, v))
                    temporal_adjacency_list[v].remove((t1, u))

        elif c == -2: # Node departure

            if sfp_special_case and batch_departure:
                if temporal_adjacency_list[source[1]]:
                    bfs_update(L, R, source,
                               batch_departure,
                               temporal_adjacency_list,
                               update_function,
                               sfp_special_case=sfp_special_case)
                for t1d, ud, vd in batch_departure:
                    temporal_adjacency_list[ud].remove((t1d, vd))
                    temporal_adjacency_list[vd].remove((t1d, ud))
                batch_departure = []


            t1,n = i[1:]
            if n in L:
                L.pop(n)
                if (len(L) == 0 and is_temporal_source == True) or (
                        len(L) == 0 and S.node_presence[source[1]][-1] == t1):
                    # If there is no possibility to extend any paths : we return the fuck out
                    return R

        # Node arrival (c == 2)
        else:
            if is_temporal_source == False:
                t0, t1, n = i[1:]
                if n == source[1]:
                    source_init_function(L,source,t0,t1)
    return R


def bfs_update(L, R, source, batch_links,
               temporal_adjacency_list,
               update_function,
               sfp_special_case = False):
    """
    Proceeds to browse every links present at instant :math:'t_0' in order to propagate the update on current possible paths

    :param L:
    :param source:
    :param batch_links:
    :param temporal_adjacency_list:
    :param L_functions:
    :return:
    """

    if sfp_special_case:
        # print()
        # for no in L:
        #     print("node :",no,"\t L[",no,"] :",L[no])
        # print("Current Batch :",batch_links)
        Q = []
        current_time = batch_links[0][0]
        for el in batch_links:
            e, u1, u2 = el[-3:]
            if u1 in L:
                su1, du1, au1 = L[u1]
                su1 = -su1
                priority = (-min(current_time, su1), du1, au1)
                heappush(Q, (priority, u1, u2, e))

            if u2 in L:
                su2, du2, au2 = L[u2]
                su2 = -su2
                priority = (-min(current_time, su2), du2, au2)
                heappush(Q, (priority, u2, u1, e))

        # Check in reachable nodes if there is one that leave later than current time
        # Only in the link departure case.
        if len(batch_links[0]) == 3:
            # We restart from the source
            if temporal_adjacency_list[source[1]]:
                for t, next in temporal_adjacency_list[source[1]]:
                    scur, dcur, acur = L[source[1]]
                    scur = -scur
                    priority = (-min(current_time, scur), dcur, acur)
                    heappush(Q, (priority, source[1], next, t))

        if Q:
            visited = set()
            while Q:
                _, pred, cur, end_link = heappop(Q)
                visited.add(pred)
                updated = update_function(L, R, current_time, end_link, pred, cur)
                if updated:
                    for t, next in temporal_adjacency_list[cur]:
                        if next not in visited:
                            scur, dcur, acur = L[cur]
                            scur = -scur
                            priority = (-min(current_time, scur), dcur, acur)
                            heappush(Q, (priority, cur, next, t))
    else:
        Q = []
        for _, e, u1, u2 in batch_links:
            if u1 in L:
                heappush(Q, (L[u1], u1, u2, e))
            if u2 in L:
                heappush(Q, (L[u2], u2, u1, e))

        if Q:
            visited = set()
            current_time = batch_links[0][0]
            while Q:
                _, pred, cur, end_link = heappop(Q)

                visited.add(pred)
                updated = update_function(L, R, current_time, end_link, pred, cur)
                if updated:
                    for t, next in temporal_adjacency_list[cur]:
                        if next not in visited:
                            heappush(Q, (L[cur], cur, next, t))


######################################
#           Foremost Paths (FoP)     #
######################################

def FoP_initialisation(S,source):
    """
    'L' initialisation for FoP

    :param source:
    :return:
    """
    L, R = {},{}
    L[source[1]] = source[0]  # (a_u)
    R[source[1]] = source[0]
    return L, R


def FoP_source_initialisation(L,source,t0,t1):
    L[source[1]] = t0

def FoP_update(L, R, t0, t1, u, v):
    """
    We update 'L' according to the link (t0,t1,u,v) and FoP properties.

    :param L:
    :param t0:
    :param t1:
    :param u:
    :param v:
    :return:
    """
    new_arrival = t0
    # If v is already in L, we've reached v at t <= t0.
    if v not in L:
        L[v] = new_arrival
        if v not in R:
            R[v] = new_arrival
    # elif t0 > L[v]:
    #     return False
    else:
        return False
    return True


def FoP_postprocess(R, source, destination=None, start_time=None):
    """
    Post process elements in 'L' to obtain time to reach 'ttr'

    :param L:
    :param source:
    :param destination: optional
    :return:
    """
    if destination is None:  #  Single source
        ttr = {k:max(v-start_time,0) for k,v in R.items()}
    else:                    # Source - Target
        ttr = {destination: max(R[destination]-start_time,0) if destination in R else math.inf}
    return ttr


def FoP(S, source, destination=None, start_time=None, E=None):
    """
    If *destination* is specified: Return the time to reach *destination* from *source* in *S*.
    Otherwise: Return the time to reach every reachable node in *S* from *source*.

    :param S: Stream graph
    :param source: temporal source node
    :param destination: optional
    :param start_time : optional, if not specified assume that the start time is source[0]
    :return:

    """

    # Checking the type of the input 'source":

    if type(source) == int:
        source = (S.node_presence[source][0], source)

    if start_time is None:
        start_time = source[0]

    is_temporal_source = True  # Times to reach always takes a temporal source

    R = L_Algorithm(S, source, FoP_initialisation, FoP_update, FoP_source_initialisation, start_time=start_time,
                    is_temporal_source=is_temporal_source, E=E)
    ttr = FoP_postprocess(R,source,destination,start_time)
    return ttr


def FoP_postprocess_pw(W):
    ttr_pw = {}
    for n in W:
        ttr_pw[n] = FoP_postprocess(W[n], n)
    return ttr_pw


def FoP_pw(S):
    """
    Return the time to reach every reachable node in *S* from every other node.
    (infinity if not reachable)
    :param S: Stream graph
    :return:
    """
    L_functions = {'initialisation': FoP_initialisation,
                   'update': FoP_update,
                   }
    W = F_Algorithm(S, L_functions)
    if W is not None:
        times_to_reach_pw = FoP_postprocess_pw(W)
        return times_to_reach_pw


######################################
#           Shortest Foremost Paths  #
######################################

def SFoP_initialisation(S,source):
    L, R = {}, {}
    L[source[1]] = (0, source[0])  # (d_v,a_v)
    R[source[1]] = (0, source[0])
    return L, R

def SFoP_source_initialisation(L,source,t0,t1):
    L[source[1]] = (0, t0)

def SFoP_update(L, R, t0, t1, u, v):
    du, au = L[u]
    new_arrival = t0
    new_distance = du + 1

    if v in R:
        current_best_distance, current_best_arrival = R[v]

        if new_arrival == current_best_arrival and new_distance < current_best_distance:
            R[v] = (new_distance, new_arrival)

        if v in L:
            dv, av = L[v]
            if new_distance < dv:
                L[v] = (new_distance, new_arrival)
            # elif new_distance > dv:
            #    return False
            else:
                return False
        else:
            L[v] = (new_distance, new_arrival)
    else:
        L[v] = (new_distance, new_arrival)
        R[v] = (new_distance, new_arrival)
    return True


def SFoP_postprocess(R, source, destination=None, start_time=None):
    '''
    Post Process elements in 'L' to obtain times to reach and lengths of shortest foremost path

    :param L:
    :param source:
    :param destination:
    :param start_time:
    :param ttr:
    :param lengths:
    :return:
    '''

    if destination is None:
        lengths = {k:v[0] for k,v in R.items()}
        ttr = {k:max(v[1]-start_time,0) for k,v in R.items()}

    else:
        lengths = {destination: R[destination][0] if destination in R else math.inf}
        ttr = {destination: max(R[destination][1] - start_time, 0) if destination in R else math.inf}

    return ttr, lengths


def SFoP(S, source, destination=None, start_time=None, E=None):


    # Checking the type of the input 'source":
    if type(source) == int:
        source = (S.node_presence[source][0], source)

    is_temporal_source = True  #  Times to reach always takes a temporal source

    if start_time is None:
        start_time = source[0]

    R = L_Algorithm(S, source, SFoP_initialisation, SFoP_update,
                    SFoP_source_initialisation, start_time=start_time,
                    is_temporal_source=is_temporal_source, E=E)
    ttr,lengths = SFoP_postprocess(R,source,destination,start_time)
    return ttr,lengths


def SFoP_W_postprocess(W):
    ttr_pw = {}
    length_pw = {}
    for n in W:
        ttr_pw[n], length_pw[n] = SFoP_postprocess(W[n], n)
    return ttr_pw, length_pw


def SFoP_pw(S):
    L_functions = {'initialisation': SFoP_initialisation,
                   'update': SFoP_update,
                   }
    W = F_Algorithm(S, L_functions)
    if W is not None:
        ttr_pw, lengths_pw = SFoP_W_postprocess(W)
        return ttr_pw, lengths_pw


######################################
#           Shortest Paths           #
######################################

def SP_initialisation(S,source):
    L,R = {},{}
    L[source[1]] = (0, source[0])  # (d,a)
    R[source[1]] = (0,source[0])
    return L, R

def SP_source_initialisation(L,source,t0,t1):
    L[source[1]] = (0, t0)  # (d,a)


def SP_update(L, R, t0, t1, u, v):
    '''
    Update Lv with the best element of Lu !
    :param L:
    :param t0:
    :param t1:
    :param u:
    :param v:
    :return:
    '''
    du, au = L[u]

    new_distance = du + 1
    new_arrival = t0

    if v in R:
        current_best_distance = R[v][0]
        if new_distance < current_best_distance:
            R[v] = (new_distance,new_arrival)

        if v in L:
            dv, av = L[v]
            if new_distance < dv:
                L[v] = (new_distance, new_arrival)
            # elif new_distance > dv: # Maybe >= (have to test with many arrival at the same time)
            #     return False
            else:
                return False
        else:
            L[v] = (new_distance, new_arrival)
    else:
        L[v] = (new_distance, new_arrival)
        R[v] = (new_distance, new_arrival)
    return True


def SP_postprocess(R,source, destination=None,
                   start_time=None):

    if destination is None:
        distances = {k:v[0] for k,v in R.items()}
    else:
        distances = {destination: R[destination][0] if destination in R else math.inf}
    return distances


def SP(S, source, destination=None, start_time=None, E=None):


    # Checking the type of the input 'source":
    is_temporal_source = True
    if type(source) == int:
        source = (S.node_presence[source][0], source)
        is_temporal_source = False

    if start_time is None:
        start_time = source[0]

    R = L_Algorithm(S, source, SP_initialisation, SP_update,
                    SP_source_initialisation, start_time=start_time,
                    is_temporal_source=is_temporal_source, E=E)
    distances = SP_postprocess(R,source,destination,start_time)
    return distances


def SP_W_postprocess(W):
    distances_pw = {}
    for n in W:
        distances_pw[n] = SP_postprocess(W[n])
    return distances_pw


def SP_pw(S):
    L_functions = {'initialisation': SP_initialisation,
                   'update': SP_update,
                   }
    W = F_Algorithm(S, L_functions)
    if W is not None:
        distances_pw = SP_W_postprocess(W)
        return distances_pw


######################################
#           Fastest Shortest Paths   #
######################################

def FSP_initialisation(S,source):
    L, R = {}, {}
    source_segment = [(t0, t1, source[1])
                      for t0, t1 in zip(S.node_presence[source[1]][::2],
                                        S.node_presence[source[1]][1::2])
                      if t0 <= source[0] <= t1][0]

    L[source[1]] = (0, -source_segment[1], source[0])  # (d_u,s_u,a_u)
    R[source[1]] = (0, source_segment[1], source[0])
    return L, R

def FSP_source_initialisation(L,source,t0,t1):
    L[source[1]] = (0,-t1,t0)

def FSP_update(L, R, t0, t1, u, v):
    du, su, au = L[u]
    su = -su
    new_arrival = t0
    new_distance = du + 1
    new_start = min(t1, su)
    new_duration = max(new_arrival - new_start, 0)

    if v in R:
        current_best_distance = R[v][0]
        current_best_duration = max(R[v][2] - R[v][1], 0)
        if new_distance < current_best_distance or (new_distance == current_best_distance and
                                                    new_duration < current_best_duration):
            R[v] = (new_distance, new_start, new_arrival)

        if v in L:
            dv, sv, av = L[v]
            sv = -sv
            if new_distance < dv or (new_distance == dv and new_start > sv):
                L[v] = (new_distance, -new_start, new_arrival)
            # elif new_distance > dv or (new_distance == dv and new_start < sv):
            #    return False
            else:
                return False
        else:
            L[v] = (new_distance, -new_start, new_arrival)
    else:
        L[v] = (new_distance, -new_start, new_arrival)
        R[v] = (new_distance, new_start, new_arrival)
    return True


def FSP_postprocess(R,source, destination=None, start_time=None):

    if destination is None:
        distances = {k:v[0] for k,v in R.items()}
        durations = {k:max(v[2]-v[1],0) for k,v in R.items()}
    else:
        distances = {destination: R[destination][0] if destination in R else math.inf}
        durations = {destination: max(R[destination][2] - R[destination][1],0) if destination in R else math.inf}
    return distances, durations


def FSP(S, source, destination=None, start_time=None, E=None):


    # Checking the type of the input 'source":
    is_temporal_source = True
    if type(source) == int:
        source = (S.node_presence[source][0], source)
        is_temporal_source = False

    if start_time is None:
        start_time = source[0]

    R = L_Algorithm(S, source, FSP_initialisation, FSP_update,
                    FSP_source_initialisation, start_time=start_time, is_temporal_source=is_temporal_source,
                    E=E)
    distances, durations = FSP_postprocess(R,source,
                                           destination,start_time)
    return distances, durations


def FSP_W_postprocess(W):
    distances_pw = {}
    duration_pw = {}
    for n in W:
        distances_pw[n], duration_pw[n] = FSP_postprocess(W[n])
    return distances_pw, duration_pw


def FSP_pw(S):
    L_functions = {'initialisation': FSP_initialisation,
                   'update': FSP_update,
                   }
    W = F_Algorithm(S, L_functions)
    if W is not None:
        distances_pw, duration_pw = FSP_W_postprocess(W)
        return distances_pw, duration_pw


######################################
#           Fastest Paths            #
######################################

#  IMPORTANT : On indexe L par -s_u au lieu de s_u pour gérer les batch arrival !!

def FP_initialisation(S,source):
    L, R = {}, {}
    source_segment = [(t0, t1, source[1])
                      for t0, t1 in zip(S.node_presence[source[1]][::2],
                                        S.node_presence[source[1]][1::2])
                      if t0 <= source[0] <= t1][0]

    L[source[1]] = (-source_segment[1], source[0])  # (-s_u,a_u)
    R[source[1]] = (source_segment[1], source[0])
    return L, R


def FP_source_initialisation(L,source,t0,t1):
    L[source[1]] = (-t1,t0)

def FP_update(L, R, t0, t1, u, v):
    '''
    Lu is sorted by the value of su and we analyze temporal links in temporal order.
    :param L:
    :param t0:
    :param t1:
    :param u:
    :param v:
    :return:
    '''
    su, au = L[u]
    su = -su  #  Trick for the sorted queue in bfs update
    new_arrival = t0
    new_start = min(t1, su)  # Remember that the aim is to leave the source node the latest possible (thus s_u
    #  can be > to t_1).
    if v in R:
        current_best_duration = max(R[v][1] - R[v][0], 0)
        if max(new_arrival - new_start, 0) < current_best_duration:
            R[v] = (new_start, new_arrival)

        if v in L:
            sv, av = L[v]
            sv = -sv
            if new_start > sv:
                L[v] = (-new_start, new_arrival)
            # elif new_start < sv:
            #     # This path is longer : don't propagate the path !
            #     return False
            else:
                return False
        else:
            L[v] = (-new_start, new_arrival)
    else:
        L[v] = (-new_start, new_arrival)
        R[v] = (new_start, new_arrival)
    return True


def FP_postprocess(R,source, destination=None, start_time=None):

    if destination is None:
        latencies = {k:max(v[1]-v[0],0) for k,v in R.items()}
    else:
        latencies = {destination:max(R[destination][1]-R[destination][0],0) if destination in R else math.inf}
    return latencies


def FP(S, source, destination=None, start_time=None, E=None):


    # Checking the type of the input 'source":
    is_temporal_source = True
    if type(source) == int:
        source = (S.node_presence[source][0], source)
        is_temporal_source = False

    if start_time is None:
        start_time = source[0]

    R = L_Algorithm(S, source, FP_initialisation, FP_update,
                    FP_source_initialisation, start_time=start_time,
                    is_temporal_source=is_temporal_source, E=E)
    latencies = FP_postprocess(R,source,destination,start_time)
    return latencies


def FP_W_postprocess(W):
    latencies_pw = {}
    for n in W:
        latencies_pw[n] = FP_postprocess(W[n])
    return latencies_pw


def FP_pw(S):
    L_functions = {'initialisation': FP_initialisation,
                   'update': FP_update,
                   }
    W = F_Algorithm(S, L_functions)
    if W is not None:
        latencies_pw = FP_W_postprocess(W)
        return latencies_pw


######################################
#           Shortest Fastest Paths   #
######################################

def SFP_initialisation(S,source):

    L, R = {},{}
    source_segment = [(t0, t1, source[1])
                      for t0, t1 in zip(S.node_presence[source[1]][::2],
                                        S.node_presence[source[1]][1::2])
                      if t0 <= source[0] <= t1][0]
    #
    # L[source[1]] = SortedSet([(source_segment[1], 0, source[0])])  # (s_u, d_u, a_u)

    # TEST SFP #
    L[source[1]] = (-source_segment[1], 0, source[0])  # (-s_u, d_u, a_u)
    # END TEST SFP#

    R[source[1]] = (source[0], source_segment[1], 0)  # (a_u, s_u, d_u)
    return L, R


def SFP_source_initialisation(L, source, t0, t1):
    # L[source[1]] = SortedSet([(t1, 0, t0)])  # (s_u, d_u, a_u)

    # TEST SFP #
    L[source[1]] = (-t1, 0, t0) #(-s_u, d_u, a_u)
    # END TEST SFP#


def SFP_update(L, R, t0, t1, u, v):
    '''
    Update L and R according to the link (t0,t1,u,v)
    :param L:
    :param R:
    :param t0:
    :param t1:
    :param u:
    :param v:
    :return:
    '''
    su, du, au = L[u]
    su = -su
    new_start = min(t1, su)
    new_arrival = t0
    new_distance = du + 1
    new_duration = max(new_arrival - new_start, 0)
    if v in R:
        current_best_duration = max(R[v][0] - R[v][1], 0)
        current_best_length = R[v][2]
        if new_duration < current_best_duration or \
                (new_duration == current_best_duration
                 and new_distance < current_best_length):
            R[v] = (new_arrival, new_start, new_distance)
        if v in L:
            sv, dv, av = L[v]
            sv = -sv
            if sv >= t0:  # Starting time in the future
                # In this mode we consider best distances for nodes in the future
                # and we ignore their potential future starting time
                if new_distance < dv or (new_distance == dv and new_start > sv):
                    L[v] = (-new_start, new_distance, new_arrival)
                elif new_distance > dv or (new_distance == dv and new_start < sv):
                    return False
                # else:
                #    return False
            else:
                if new_start > sv or (new_start == sv and new_distance < dv):
                    L[v] = (-new_start, new_distance, new_arrival)
                # elif new_start < sv:# or (new_start == sv and new_distance > dv):
                #    return False
                else:
                    return False
        else:
            L[v] = (-new_start, new_distance, new_arrival)
    else:
        L[v] = (-new_start, new_distance, new_arrival)
        R[v] = (new_arrival, new_start, new_distance)
    return True


# BEFORE TEST
def SFP_update_old(L, R, t0, t1, u, v):
    '''
    Update Lv from (t0,t1,u,v) according to SFP constraints.
    Dans L_u pour les SFP on a : soit un seul triplet et a_u > s_u either many triplets and a_u <= s_u,     they are sorted by their departure time
    :param L:
    :param t0:
    :param t1:
    :param u:
    :param v:
    :return:
    '''
    # Cleaning L[u] :
    if len(L[u]) > 1:
        bu_pos = min(L[u].bisect_left((t0, 0, 0)), len(L[u]) - 1)
        if bu_pos > 0:
            # print(" u :", u)
            # print(" current time :", t0)
            # print("\t b pos cleaning u:", bu_pos)
            # print("\t before L u:", L[u])
            L[u].difference_update(L[u][:bu_pos])
            # print("\t after L u:", L[u])
        # TEST !! #
        if len(L[u]) > 1:
            # print("current time :", t0)
            # print("u , L[u] :", u, L[u])
            for i in range(len(L[u]) - 1):
                s_pred, d_pred, a_pred = L[u][i]
                s_next, d_next, a_next = L[u][i + 1]
                assert max(t0 - s_pred, 0) == max(t0 - s_next, 0)
                assert s_pred < s_next
                assert d_pred < d_next


    if v in L and len(L[v]) > 1:
        # Cleaning L[v]
        bv_pos = min(L[v].bisect_left((t0, 0, 0)), len(L[v]) - 1)
        if bv_pos > 0:
            # print("v :", v)
            # print(" current time :", t0)
            # print("\t b pos cleaning v:", bv_pos)
            # print("\t before L v:", L[v])
            L[v].difference_update(L[v][0:bv_pos])
            # print("\t after L v:", L[v])
        # TEST !! #
        if len(L[v]) > 1:
            # print("current time :", t0)
            # print("u , L[u] :", v, L[v])
            for i in range(len(L[v]) - 1):
                s_pred, d_pred, a_pred = L[v][i]
                s_next, d_next, a_next = L[v][i + 1]
                assert max(t0 - s_pred, 0) == max(t0 - s_next, 0)
                assert s_pred < s_next
                assert d_pred < d_next


    updated = False
    # print("\t u,v :",u,v)
    for su, du, au in L[u]: # As we have clean L_u we can consider every elements

        # print("\t au,su,du :",au,su,du)
        new_arrival = t0
        new_start = min(t1, su)
        new_distance = du + 1
        new_duration = max(new_arrival - new_start, 0)
        # print("new a, new s, new d :",new_arrival,new_start,new_distance)

        if v in R:

            current_best_duration = max(R[v][0] - R[v][1], 0)
            current_best_length = R[v][2]
            if new_duration < current_best_duration or \
                    (new_duration == current_best_duration
                     and new_distance < current_best_length):
                R[v] = (new_arrival, new_start, new_distance)
            if v in L:

                to_add = False
                to_remove = set()
                for sv,dv,av in L[v]: # as we have clean L_v we can consider every elements

                    # if new_start > sv or (new_start == sv and new_distance < dv) or \
                    #         (new_duration == 0 and new_distance < dv):
                    #             L[v].add((new_start,new_distance,new_arrival))
                    if new_start > sv:
                        if new_distance <= dv:
                            to_add = True
                            to_remove.add((sv,dv,av))

                        elif (sv,dv,av) == L[v][-1]: # IMPORTANT take care of the else
                            to_add = True

                    elif new_start < sv:
                        if new_distance < dv:
                            to_add = True
                            break
                        elif new_distance == dv:
                            updated = True
                            break
                        else:
                            break

                    else: # new_start = sv
                        if new_distance < dv:
                            to_remove.add((sv,dv,av))
                            to_add = True
                        elif new_distance == dv:
                            updated = True
                            break
                        else:
                            break

                if to_remove:
                    L[v].difference_update(to_remove)
                if to_add:
                    L[v].add((new_start, new_distance, new_arrival))
                    updated = True
            else:
                L[v] = SortedSet([(new_start, new_distance, new_arrival)])
                updated = True
        else:
            L[v] = SortedSet([(new_start, new_distance, new_arrival)])
            R[v] = (new_arrival, new_start, new_distance)
            updated = True



    return updated


def SFP_postprocess(R,source, destination=None,start_time = None):

    if destination is None:
        latencies = {k: max(v[0]-v[1],0) for k,v in R.items()}
        lengths = {k:v[2] for k,v in R.items()}
    else:
        latencies = {destination:max(R[destination][0]-R[destination][1],0) if destination in R else math.inf}
        lengths = {destination:R[destination][2] if destination in R else math.inf}
    return latencies, lengths


def SFP(S, source, destination=None, start_time=None, E=None):


    # Checking the type of the input 'source":
    is_temporal_source = True
    if type(source) == int:
        source = (S.node_presence[source][0], source)
        is_temporal_source = False

    if start_time is None:
        start_time = source[0]

    R = L_Algorithm(S, source, SFP_initialisation, SFP_update,
                    SFP_source_initialisation, start_time=start_time,
                    is_temporal_source=is_temporal_source, sfp_special_case=True,
                    E=E)
    latencies, lengths = SFP_postprocess(R,source,destination,start_time)

    assert latencies[source[1]] == 0
    assert lengths[source[1]] == 0
    return latencies, lengths


def SFP_W_postprocess(W):
    latencies_pw = {}
    lengths_pw = {}
    for n in W:
        latencies_pw[n], lengths_pw[n] = SFP_postprocess(W[n])
    return latencies_pw, lengths_pw


def SFP_pw(S):
    L_functions = {'initialisation': SFP_initialisation,
                   'update': SFP_update,
                   }
    W = F_Algorithm(S, L_functions)
    if W is not None:
        latencies_pw, lengths_pw = SFP_W_postprocess(W)
        return latencies_pw, lengths_pw
