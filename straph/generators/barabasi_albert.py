import time, random, warnings
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from collections import defaultdict

from straph import stream as sg
from straph.utils import profile_shit




def random_node_presence(t_windows, rep, plac, dur):
    """
    Generate the occurrence and the presence of a node given occurrence_law(occurrence_param)
    and presence_law(presence_param).


    :param t_windows: Time window of the Stream Graph
    :param rep: Number of segmented nodes
    :param plac: Emplacement of each segmented node (sorted array)
    :param dur: Length of each interval corresponding to a segmented node
    :return: node presence
    """
    acc = plac[0] + dur[0]
    if acc > t_windows[1]:
        # Si on dépasse la borne supérieure de temps de notre Stream Graph
        acc = t_windows[1]
        n_presence = [plac[0], acc]
        return n_presence
    n_presence = [plac[0], acc]  #  Initialisation, [t0,t0+d0]
    for i in range(1, rep):
        acc = plac[i] + dur[i]
        if acc > t_windows[1]:
            # Si on dépasse la borne supérieure de temps de notre Stream Graph
            acc = t_windows[1]
        if acc <= n_presence[-1]:
            # Cas ou l'intervalle est inclus dans le précédent
            continue
        if plac[i] <= n_presence[-1]:
            # Cas de l'intersection : on fusionne
            n_presence[-1] = acc
        else:
            # Cas disjoint
            n_presence += [plac[i], acc]
    return n_presence


def random_link_presence(intersec, rep, dur):
    """
    Generate the occurrence and the presence of a link given occurrence_law(occurrence_param)
    and presence_law(presence_param).

    :param intersec: Intersection between the prense of extremities (realisable interval for link))
    :param rep: Number of segmented links
    :param dur: Length of each interval corresponding to a segmented link
    :return: link presence
    """
    plac = []
    id_intersec = {}  #  Pour retenir l'intervalle d'intersection ou se situe le lien (#BlackMagic)
    id_plac = np.random.randint(len(intersec) // 2, size=rep)  # On choisit aléatoirement un temps commun aux 2 noeuds
    for i in id_plac:
        i = 2*i
        T0 = intersec[i]
        T1 = intersec[i + 1]
        plac.append(np.random.uniform(T0, T1))
        id_intersec[plac[-1]] = i + 1
    plac.sort()
    # Initialisation avant check puis fusion si necessaire
    acc = plac[0] + dur[0]
    if acc > intersec[id_intersec[plac[0]]]:
        acc = intersec[id_intersec[plac[0]]]
    l_presence = [plac[0], acc]

    for i in range(1, rep):
        acc = plac[i] + dur[i]
        if acc > intersec[id_intersec[plac[i]]]:
            # Si on dépasse la borne supérieure de temps de notre fenêtre
            acc = intersec[id_intersec[plac[i]]]
        if acc <= l_presence[-1]:
            # Cas ou l'intervalle est inclus dans le précédent (cas d'un accroissement faible)
            continue
        if plac[i] <= l_presence[-1]:
            # Cas de l'intersection : on fusionne
            l_presence[-1 ] = acc
        else:
            # Cas disjoint
            l_presence += [plac[i], acc]
    return l_presence


def get_intersection(u, v, node_presence):
    """
    Get the intersection between the presence of u and v.

    :param u: First Node
    :param v: Second Node
    :param node_presence: Node presence
    :return: Interection
    """
    intersec = []
    for ut0, ut1 in zip(node_presence[u][::2], node_presence[u][1::2]):
        for vt0, vt1 in zip(node_presence[v][::2], node_presence[v][1::2]):
            if ut0 <= vt1 and vt0 <= ut1:
                intersec += [max(ut0, vt0), min(vt1, ut1)]
    return intersec


def barabasi_albert(t_window,
                    nb_node,
                    occurrence_law_node = 'poisson',
                    occurrence_param_node = None,
                    presence_law_node = 'poisson',
                    presence_param_node = None,
                    occurrence_law_link = 'poisson',
                    occurrence_param_link= None,
                    presence_law_link = 'poisson',
                    presence_param_link = None,
                    initial_nb_node = 1,
                    m_link = 1,
                    ):
    """
    Stream Graph generator following a Barabasi-Albert like behavior. The stream graph begins with *initial_nb_node*
    all connected, then *nb_node*-*initial_nb_node* are added they are connected with *m_link* existing node
    with probability p = degree_of_node/sum_of_degrees.
    Each node occurs following *occurrence_law_node(occurrence_param_node)*, each segmented node has a presence
    which length follows *presence_law_node(presence_param_node)*, each link occurs following *occurrence_law_link(
    occurrence_param_link)* and each segmented link has a presence which length follows
     *presence_law_link(presence_param_link)*.

    :param t_window: Time windows of the Stream Graph
    :param nb_node: Desired Number of nodes in the stream graphs
    :param occurrence_law_node: Random variable for node occurrence (numpy function or 'poisson'
    :param occurrence_param_node: Parameter of the node occurrence law
    :param presence_law_node: Random variable for node presence (numpy function or 'poisson' or 'uniform')
    :param presence_param_node: Parameter of the node presence law
    :param occurrence_law_link: Random variable for link occurrence (numpy function or 'poisson')
    :param occurrence_param_link: Parameter of the link occurrence law
    :param presence_law_link: Random variable for link presence (numpy function or 'poisson' or 'uniform')
    :param presence_param_link: Parameter of the link presence law
    :param initial_nb_nodes: Parameter of the Barabasi-Albert model : Number of connected nodes initial present (m0)
    :param m_link: Parameter of the Barabasi-Albert model : the number of edge when we (preferentially) attach it.
    :return:
    """
    if nb_node > 200000:
        warnings.warn("The number of nodes is probably to big, it may lead to a memory error !")

    # Default parameters
    if occurrence_param_node is None:
        occurrence_param_node = 1
    if presence_param_node is None:
        presence_param_node = (t_window[1] - t_window[0]) / 2
    if occurrence_param_link is None:
        occurrence_param_link = 1
    if presence_param_link is None:
        presence_param_link = (t_window[1] - t_window[0]) / 3

    nodes,node_presence = generate_node_presence(t_window, nb_node,
                                                 occurrence_law_node,
                                                 occurrence_param_node,
                                                 presence_law_node,
                                                 presence_param_node)

    links,link_presence = generate_link_presence(m_link, initial_nb_node, nb_node, node_presence,
                                                 occurrence_law_link,
                                                 occurrence_param_link, presence_law_link,
                                                 presence_param_link)

    S = sg.stream_graph(times=t_window,
                        nodes=nodes,
                        node_presence=node_presence,
                        links=links,
                        link_presence=link_presence,
                        node_to_label={n: str(n) for n in nodes}

                        )
    return S


def generate_link_presence(m_link,
                           initial_nb_nodes,
                           nb_node, node_presence,
                           occurrence_law_link,
                           occurrence_param_link,
                           presence_law_link,
                           presence_param_link):
    '''
    Generate links presence and occurrence.

    :param m_link: Parameter of the Barabasi-Albert model : the number of edge when we (preferentially) attach it.
    :param initial_nb_nodes: Parameter of the Barabasi-Albert model : Number of connected nodes initial present (m0)
    :param nb_node: Number of Nodes
    :param node_presence: Node presence
    :param occurrence_law_link: Random variable for link occurrence (numpy function or 'poisson')
    :param occurrence_param_link: Parameter of the link occurrence law
    :param presence_law_link: Random variable for link presence (numpy function or 'poisson' or 'uniform')
    :param presence_param_link: Parameter of the link presence law
    :return:


    References
    ----------
    .. [1] A. L. Barabási and R. Albert "Emergence of scaling in
       random networks", Science 286, pp 509-512, 1999.
    '''

    if m_link > initial_nb_nodes:
        raise ValueError("The number of link to generate :" + str(m_link) +
                         " is higher than the initial number of nodes :" + str(initial_nb_nodes))

    links = []
    link_presence = []


    # Add the number of initial nodes (m0 in barabasi model)
    # Initially the graph is connected
    edges = [(i,j) for i,j in zip(range(initial_nb_nodes-1),range(1,initial_nb_nodes))]
    #  List of existing nodes, with nodes repeated a number of time equal to their degree
    repeated_nodes = list(range(initial_nb_nodes))
    for i in range(initial_nb_nodes, nb_node):
        # i: the current number of node in the graph
        targets = set()
        while len(targets) != m_link:
            id = random.sample(repeated_nodes,1)
            targets.add(repeated_nodes[id[0]])
        edges += [(i,t) for t in targets]
        repeated_nodes+=list(targets)+[i]*m_link
    nb_edge = len(edges)

    # Random sampling for the link's occurrence.
    if type(occurrence_law_link) == str:
        if occurrence_law_link == 'poisson':
            occurrences = np.maximum(np.random.poisson(occurrence_param_link, nb_edge), np.ones(nb_edge))
        else:
            raise ValueError("The random distribution " + str(occurrence_law_link) + " is not supported."
                                                                                   " Try a numpy function directly.")
    else:
        occurrences = occurrence_law_link(occurrence_param_link, nb_edge)
    #  For each interval of presence we sample a random duration.
    if type(presence_law_link) == str:
        if presence_law_link == 'uniform':
            durations = list(np.random.uniform(0, 2 * presence_param_link, int(sum(occurrences))))
        elif presence_law_link == 'poisson':
            durations = list(np.random.poisson(presence_param_link, int(sum(occurrences))))
        else:
            raise ValueError("The random distribution "+str(presence_law_link)+" is not supported."
                                                                               " Try a numpy function directly.")
    else:
        durations = list(presence_law_link(presence_param_link, sum(int(occurrences))))

    cnt = 0
    for i, (u, v) in enumerate(edges):
        rep = int(occurrences[i])
        dur = durations[cnt:cnt + rep]
        cnt += rep
        intersec = get_intersection(u, v, node_presence)
        if intersec:
            links.append((u, v))
            link_presence.append(random_link_presence(intersec, rep, dur))
    return links,link_presence


def generate_node_presence(t_window, nb_node, occurrence_law_node, occurrence_param_node,
                           presence_law_node, presence_param_node):
    '''
    Generate nodes presence and occurrence.

    :param t_window: Time windows of the Stream Graph
    :param nb_node: Number of Nodes
    :param occurrence_law_node: Random variable for node occurrence (numpy function or 'poisson'
    :param occurrence_param_node: Parameter of the node occurrence law
    :param presence_law_node: Random variable for node presence (numpy function or 'poisson' or 'uniform')
    :param presence_param_node: Parameter of the node presence law
    :return:
    '''

    # Nodes initialisation
    nodes = list(range(nb_node))
    # Node presence initialisation
    node_presence = [[] for n in nodes]
    # Random sampling for the node's occurrence.
    if type(occurrence_law_node) == str:
        if occurrence_law_node == 'poisson':
            occurrences = np.maximum(np.random.poisson(occurrence_param_node, nb_node), np.ones(nb_node))
        else:
            raise ValueError("The random distribution " + str(occurrence_law_node) + " is not supported."
                                                                                   " Try a numpy function directly.")
    else:
        occurrences = np.maximum(occurrence_law_node(occurrence_param_node, nb_node), np.ones(nb_node))
    #  For each interval of presence we sample a random duration.
    if type(presence_law_node) == str:
        if presence_law_node == 'uniform':
            durations = np.random.uniform(0, 2 * presence_param_node, int(sum(occurrences)))
        elif presence_law_node == 'poisson':
            durations = np.random.poisson(presence_param_node, int(sum(occurrences)))
        else:
            raise ValueError("The random distribution "+str(presence_law_node)+" is not supported."
                                                                               " Try a numpy function directly.")
    else:
        durations = presence_law_node(presence_param_node, int(sum(occurrences)))
    # Placement aléatoire des intervalles sur la fenêtre de temps, fusion des intervalles si nécessaire
    # X ~ (Uniform(T0,T1))
    emplacements = np.random.uniform(t_window[0], t_window[1], int(sum(occurrences)))
    cnt = 0
    for n in nodes:
        rep = int(occurrences[n])
        dur = durations[cnt:cnt + rep]
        plac = np.sort(emplacements[cnt:cnt + rep])
        node_presence[n] = random_node_presence(t_window, rep, plac, dur)
        cnt += rep

    return nodes,node_presence