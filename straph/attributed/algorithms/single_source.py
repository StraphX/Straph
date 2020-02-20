from collections import defaultdict
from heapq import *
import math

from sortedcollections import SortedSet


##################################################
#           Temporary Functions                  #
##################################################


def get_temporal_sources(S, source):
    '''
    Return the maximal segmented nodes (t0,t1,u) s.t u = source.

    :param S: Stream Graph
    :param source: source node
    :return:
    '''
    return [(t0, t1, source) for t0, t1 in zip(S.node_presence[source][::2], S.node_presence[source][1::2])]


#####################################
#   L-Algorithms source node        #
#####################################

def L_Algorithm(S, source, initialisation_function, update_function, start_time=None):
    """
        An implementation of L-Algorithm (described in ~add ref) !! DIRECTED VERSION !!!
        Single Source Algorithm to compute temporal paths in Stream Graph

        :param S: A Stream Graph
        :param source: A temporal node (t,x) \in W or a node x \in V
        :param initialisation_function: Initialisation according to the path's type (supported :
                - foremost path
                - shortest foremost path
                - shortest path
                - fastest shortest path
                - fastest path
                - shortest fastest path)
        :param update_function: Update according to the path's type (supported :
                - foremost path
                - shortest foremost path
                - shortest path
                - fastest shortest path
                - fastest path
                - shortest fastest path)
        :return: R containing the best results
        # IMPORTANT : We suppose that we are in the WCC of 'source' otherwise it's fucking expensive !
        """

    # Checking the type of the input 'source":
    is_temporal_source = True
    if type(source) == int:
        source = (S.node_presence[source][0], source)
        is_temporal_source = False

    if start_time is None:
        start_time = source[0]
    #
    # Initialisation:
    E = S.ordered_events()
    temporal_adjacency_list = defaultdict(set)  # key : temporal node, value : (temporal node, end time of link)
    batch_arrival = []
    L, R = initialisation_function(S, source)
    t_last_arrival = None

    for i in E:
        # print()
        # print("L :",L)
        # print("R :",R)
        # print("i :",i)

        c = i[0]

        # Link ARRIVAL
        if c == 1:
            t1 = i[2]
            if t1 >= start_time:
                t0, t1, u, v, w, d = i[1:]
                l = i[1:]
                if t_last_arrival is None:
                    t_last_arrival = t0
                if t0 == t_last_arrival:
                    temporal_adjacency_list[u].add((t1, v, w, d))
                    # temporal_adjacency_list[v].add((t1, u, w, d))
                    batch_arrival.append(l)

                else:
                    if batch_arrival:
                        bfs_update(L, R, source, batch_arrival,
                                   temporal_adjacency_list,
                                   update_function)
                    temporal_adjacency_list[u].add((t1, v, w, d))
                    # temporal_adjacency_list[v].add((t1, u, w, d))
                    batch_arrival = [l]
                    t_last_arrival = t0

        # Link DEPARTURE
        elif c == -1:
            t1 = i[1]
            if t1 >= start_time:
                t1, u, v, w, d = i[1:]
                #
                if batch_arrival:
                    bfs_update(L, R, source, batch_arrival,
                               temporal_adjacency_list,
                               update_function)
                    batch_arrival = []
                #
                temporal_adjacency_list[u].remove((t1, v, w, d))
                # temporal_adjacency_list[v].remove((t1, u, w, d))

        # Node Departure
        elif c == -2:
            t1, n = i[1:]
            if n in L:
                L.pop(n)  # The node is not present anymore, we remove it
                if len(L) == 0 and (is_temporal_source == True or S.node_presence[source[1]][-1] == t1):
                    # If there is no possibility to extend any paths : we return the fuck out
                    return R

        # Node arrival (c == 2)
        else:
            if is_temporal_source == False:
                t0, t1, n = i[1:]
                if n == source[1]:
                    L[n] = SortedSet([(t0, t1, 0, 0)], key=lambda x: x[0])  # a,s,d,l
    return R


def bfs_update(L, R, source, batch_arrival, temporal_adjacency_list, update_function):
    """
    Proceeds to browse every links present at instant :math:'t_0' in order to propagate the update on current possible paths

    :param L:
    :param source:
    :param batch_arrival:
    :param temporal_adjacency_list:
    :param L_functions:
    :return:
    """
    Q = []
    visited = {source[1]}  #  Its pointless to come back to the source !
    begin_link = batch_arrival[0][0]
    for _, e, u1, u2, w, d in batch_arrival:
        #  We check if there is an arrival time compatible the the ending of the current link !!
        #  there is 'a' such as bu1 <= a <= e !!
        if u1 in L:
            # print("u1 :",u1," Lu1 :",L[u1])
            #  (maximal compatible time of arrival, minimal distance according to 'e')
            e_pos = L[u1].bisect_right((e, 0, 0))
            if e_pos > 0:
                for au1, su1, du1, lu1 in L[u1][:e_pos]:
                    #assert au1 <= e  #  There exists an arrival before the end of the link
                    priority = (du1, au1, su1, lu1)  # TODO : L'ordre DEPEND DU TYPE DE CHEMIN !!!
                    heappush(Q, (priority, u1, u2, e, w, d))

        # if u2 in L:
        #     for du2,au2 in L[u2]: # remplacer par des bissect_right
        #         if au2 <= e:
        #             heappush(Q, ((du2,au2), u2, u1, e, w, d))
        #             break
    if Q:
        while Q:
            priority, pred, cur, end_link, w, d = heappop(Q)
            visited.add(pred)

            cur_link = (begin_link, end_link, pred, cur, w, d)
            updated = update_function(L, R, cur_link, priority)
            if updated:
                for e, next, w, d in temporal_adjacency_list[cur]:
                    if next not in visited:
                        #  (maximal compatible time of arrival, minimal distance according to 'e')
                        e_pos = L[cur].bisect_right((e, 0, 0))
                        if e_pos > 0:
                            # acur, scur, dcur, lcur = L[cur][e_pos - 1]
                            for acur, scur, dcur, lcur in L[cur][:e_pos]:
                                #assert acur <= e  #  There exists an arrival before the end of the link
                                priority = (dcur, acur, scur, lcur) # TODO : L'ordre DEPEND DU TYPE DE CHEMIN !!!
                                heappush(Q, (priority, cur, next, e, w, d))


#####################
#  Paths functions  #
#####################

######################################
#       Shortest Paths
#####################################

def SP_initialisation(S,source):
    '''
    :param source: temporal source node (t,x) \in W
    :return:
    '''
    L, R = {}, {}

    source_segment = [(t0, t1, source[1]) for t0, t1 in zip(S.node_presence[source[1]][::2],
                                                            S.node_presence[source[1]][1::2])
                      if t0 <= source[0] <= t1][0]

    L[source[1]] = SortedSet([(source[0], source_segment[1], 0, 0)], key=lambda x: x[0])  # (a,s,d,l)
    R[source[1]] = (0, source[0], source_segment[1], 0)  # (d,a,s,l)
    return L, R


def SP_update(L, R, link, priority):
    '''
    Update Lv with the best element of Lu !
    :param L: Bigga Structure
    :param t0: current time
    :param t1: end of link
    :param u: node already updated (b0,e0,n0)
    :param v: node to update (b0,e0,n1)
    :param w: weight
    :param tr: traversal time
    :return:
    '''
    b, e, u, v, w, d = link
    # print("current link in update : ",link)
    du, au, su, lu = priority  #  Priority is like R
    # assert au <= e
    new_arrival = max(au, b) + d
    new_distance = du + w
    new_start = min(su, e - lu) # End of link minus the travel time to reach it (maximal starting time).
    new_trip_time = lu + d
    # print("au,su,du :",au,su,du)
    # print("av,sv,dv ? :",new_arrival,new_start,new_distance)
    if v in R:
        best_distance = R[v][0]  # in R[v] -> (d,a,s)
        if new_distance < best_distance:
            R[v] = (new_distance, new_arrival, new_start, new_trip_time)

        if v in L:
            #  (maximal compatible time of arrival, minimal distance according to 'e')
            pos = L[v].bisect_right((new_arrival, 0, 0))
            if pos > 0:
                best_dv = L[v][pos - 1][2]  #  in L[v] -> (a,s,d,l)
                if new_distance < best_dv:
                    # We remove useless elements before the current time: b != new_arrival as d can be > 0
                    pos_clean = L[v].bisect_right((b, 0, 0))
                    L[v].difference_update(L[v][0:pos_clean])
                    #

                    L[v].add((new_arrival, new_start, new_distance, new_trip_time))
                else:
                    return False

            else:  # This new_arrival comes before all the other: we must add it
                L[v].add((new_arrival, new_start, new_distance, new_trip_time))
        else:
            L[v] = SortedSet([(new_arrival, new_start, new_distance, new_trip_time)], key=lambda x: x[0])
    else:
        L[v] = SortedSet([(new_arrival, new_start, new_distance, new_trip_time)], key=lambda x: x[0])
        R[v] = (new_distance, new_arrival, new_start, new_trip_time)
    return True


# def SP_postprocess(R, destination=None,
#                    distances=None):
#     if destination is None:
#         if distances is None:
#             distances = {}
#         for k, v in R.items():
#             n = k[2]
#             if n in distances:
#                 distances[n] = min(distances[n], v[0])
#             else:
#                 distances[n] = v[0]
#     else:
#         if distances is None:
#             distances = math.inf
#         if type(destination) == int:  # destination node
#             for k, v in R.items():
#                 if k[2] == destination:
#                     distances = min(distances, v[0])
#         else:  # temporal destination node
#             for k, v in R.items():
#                 if k[2] == destination[2] and k[0] <= destination[0] <= destination[1] <= k[1]:
#                     distances = min(distances, v[0])
#
#     return distances


# Shortest paths
def SP(S, source, destination=None,start_time = None):
    distances = None
    R = L_Algorithm(S, source, SP_initialisation, SP_update, start_time)
    if R is not None:
        print("R :", R)
        print("Number of destination :",len(R)," nb nodes :",len(S.nodes))
        print("Max of distances from ", source, " :",max([v[0] for v in R.values()]))
        if destination and destination in R:
            return R[destination][0]
        distances = {k: v[0] for k, v in R.items()}
    return distances


import matplotlib.pyplot as plt
import numpy as np
import random
import straph.generators as gsg
from straph.attributed.attributed_test import plot
from straph.stream import read_stream_graph

from straph.utils import profile_shit

if __name__ == '__main__':
    ## TEST single source node
    T = [0, 1000]
    nb_node = 1000
    occurrence_law = 'poisson'
    presence_law = 'poisson'
    occurrence_param_node = 5
    presence_param_node = 200
    occurrence_param_link = 15
    presence_param_link = 20

    weights_law = 'poisson'
    weights_param = 2
    trips_law = 'uniform'
    trips_law_param = 2

    p_link = 1 * np.sqrt(nb_node) / nb_node

    random.seed(1)
    np.random.seed(1)

    S = gsg.erdos_renyi(T,
                        nb_node,
                        occurrence_law,
                        occurrence_param_node,
                        presence_law,
                        presence_param_node,
                        occurrence_law,
                        occurrence_param_link,
                        presence_law,
                        presence_param_link,
                        p_link,
                        # directed=True,
                        weights_law=weights_law,
                        weights_law_param=weights_param,
                        trips_law=trips_law,
                        trips_law_param=trips_law_param,
                        )

    __directory__ = "/home/leo/Dev/CODE/Straph/examples/path_examples/"
    __file__ = "attributed_example_2"
    #
    S = read_stream_graph(path_nodes=__directory__ + __file__ + "_nodes.sg",
                          path_links=__directory__ + __file__ + "_links.sg",
                          path_weights=__directory__ + __file__ + "_weights.sg",
                          path_trips=__directory__ + __file__ + "_trips.sg")

    # weights = []
    # trips = []
    # for l,lp in zip(S.links,S.link_presence):
    #     w_list = []
    #     t_list = []
    #     for _ in range(len(lp)//2):
    #         w_list.append(np.random.randint(1,10))
    #         t_list.append(np.random.randint(1,5))
    #     weights.append(w_list)
    #     trips.append(t_list)
    # S.weights = weights
    # S.trips = trips
    # print("link :", S.links)
    # print("link presence :", S.link_presence)
    # print("weigths :", S.weights)
    # print("link presence :", S.trips)
    S.describe()
    print()
    print("SOURCE :", 0)

    # profile_shit("SP(S,0)",snakeviz=True)

    distances = SP(S, 0)  # source x = 0
    print()
    # print("SOURCE :", (0, 0))
    # distances = SP(S, (0, 0))  #  source (t,x) = (1,0)

    plot(S)
    plt.show()

    ## Test single source temporal node
