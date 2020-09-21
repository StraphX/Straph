import time, random, warnings
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from straph import stream as sg
from straph.utils import profile_shit


def random_link(p_link):
    """
    Generate a random link with probability p_link.

    :param p_link: probability (0<= p_link <1).
    :return: True if the link exists False otherwise
    """
    if random.random() < p_link:
        return True
    return False


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
        i = 2 * i
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
            l_presence[-1] = acc
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


def erdos_renyi(t_window,
                nb_node,
                occurrence_law_node='poisson',
                occurrence_param_node=None,
                presence_law_node='poisson',
                presence_param_node=None,
                occurrence_law_link='poisson',
                occurrence_param_link=None,
                presence_law_link='poisson',
                presence_param_link=None,
                p_link=None,
                directed=False,
                weights_law=False,
                weights_law_param=False,
                trips_law=False,
                trips_law_param=False,
                ):
    """
    Stream Graph generator following an Erdos-Renyi like behavior. Each link is sampled with probability *p_link*,
    each node occurs following *occurrence_law_node(occurrence_param_node)*, each segmented node has a presence
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
    :param p_link: Probability of the existence of a link between two random nodes (0<=p_link<1), Erdos Renyi parameter.
    :param directed: True for a random directed Stream Graph, False for an undirected Stream Graphs
    :return:
    """
    if nb_node > 200000:
        warnings.warn("The number of nodes is probably to big, it may lead to a memory error !")

    # Default parameters
    if p_link is None:
        p_link = np.sqrt(nb_node) / nb_node
    if occurrence_param_node is None:
        occurrence_param_node = 1
    if presence_param_node is None:
        presence_param_node = (t_window[1] - t_window[0]) / 2
    if occurrence_param_link is None:
        occurrence_param_link = 1
    if presence_param_link is None:
        presence_param_link = (t_window[1] - t_window[0]) / 3

    nodes, node_presence = generate_node_presence(t_window, nb_node,
                                                  occurrence_law_node,
                                                  occurrence_param_node,
                                                  presence_law_node,
                                                  presence_param_node)

    links, link_presence = generate_link_presence(nb_node, node_presence,
                                                  occurrence_law_link,
                                                  occurrence_param_link, p_link,
                                                  presence_law_link, presence_param_link,
                                                  directed, weights_law,
                                                  weights_law_param,
                                                  trips_law,
                                                  trips_law_param, )
    weights = None
    if weights_law:
        weights = generate_weights(link_presence, weights_law, weights_law_param)

    trips = None
    if type(trips_law) == str:
        trips = generate_trips(links, link_presence, node_presence, trips_law, trips_law_param)

    S = sg.stream_graph(times=t_window,
                        nodes=nodes,
                        node_presence=node_presence,
                        links=links,
                        link_presence=link_presence,
                        node_to_label={n: str(n) for n in nodes},
                        weights=weights,
                        trips=trips,
                        )
    return S


def get_node_presence_from_link(node_presence, n, b, e):
    """
    Return the maximal temporal node corresponding to (b,e,n)
    :param n: node
    :param b: beginning of the interval (time)
    :param e: ending of the interval (time)
    :return: Maximal temporal node presence corresponding to (b,e,n) : (t0,t1)
    """
    for t0, t1 in zip(node_presence[n][::2], node_presence[n][1::2]):
        if t0 <= b and e <= t1:
            return (t0, t1)
    return None


def generate_trips(links, link_presence, node_presence, trips_law, trips_law_param):
    occurences = sum([len(lp) // 2 for lp in link_presence])
    if trips_law == 'uniform':
        trips = list(np.random.uniform(0, 2 * trips_law_param, occurences))
    elif trips_law == 'poisson':
        trips = list(np.random.poisson(trips_law_param, occurences))
    else:
        raise ValueError("The random distribution " + str(trips_law) + " is not supported for traversal times."
                                                                       " Try a numpy function directly.")
    #  We need to assert that the node are present during the trip+duration of the link !!!
    trips_list = []
    i = 0
    for l, lp in zip(links, link_presence):
        u, v = l
        # Check node presence of u
        current_trips = []

        for lt0, lt1, tr in zip(lp[::2], lp[1::2], trips[i:i + len(lp) // 2]):
            ut0, ut1 = get_node_presence_from_link(node_presence, u, lt0, lt1)
            vt0, vt1 = get_node_presence_from_link(node_presence, v, lt0, lt1)
            current_trips.append(min(tr, min(ut1 - lt1, vt1 - lt1)))
        i += len(lp) // 2
        trips_list.append(current_trips)

    return trips_list


def generate_weights(link_presence, weights_law, weights_law_param):
    occurences = sum([len(lp) // 2 for lp in link_presence])
    if type(weights_law) == str:
        if weights_law == 'uniform':
            weights = list(np.maximum(np.random.uniform(0, 2 * weights_law_param, occurences),np.ones(occurences)))
        elif weights_law == 'poisson':
            weights = list(np.maximum(np.random.poisson(weights_law_param, occurences),np.ones(occurences)))
        else:
            raise ValueError("The random distribution " + str(weights_law) + " is not supported for weights."
                                                                             " Try a numpy function directly.")
        # We distribute these weights into a list of list corresponding to the links
        weights_list = []
        i = 0
        for lp in link_presence:
            weights_list.append(weights[i:i + len(lp) // 2])
            i += len(lp) // 2

        return weights_list


def generate_link_presence(nb_node, node_presence, occurrence_law_link,
                           occurrence_param_link, p_link, presence_law_link, presence_param_link,
                           directed=False,
                           weights_law=False,
                           weights_law_param=False,
                           trips_law=False,
                           trips_law_param=False, ):
    '''
    Generate links presence and occurrence.

    :param nb_node: Number of Nodes
    :param node_presence: Node presence
    :param occurrence_law_link: Random variable for link occurrence (numpy function or 'poisson')
    :param occurrence_param_link: Parameter of the link occurrence law
    :param p_link: Probability of the existence of a link between two random nodes (0<=p_link<1), Erdos Renyi parameter
    :param presence_law_link: Random variable for link presence (numpy function or 'poisson' or 'uniform')
    :param presence_param_link: Parameter of the link presence law
    :param directed: True for a random directed Stream Graph, False for an undirected Stream Graphs
    :return:

    References
    ----------
    .. [1] P. Erdős and A. Rényi, On Random Graphs, Publ. Math. 6, 290 (1959).
    .. [2] E. N. Gilbert, Random Graphs, Ann. Math. Stat., 30, 1141 (1959).
    '''

    links = []
    link_presence = []

    # Sample links according to p_link
    if directed:
        max_nb_edge = (nb_node * (nb_node - 1))
    else:
        max_nb_edge = (nb_node * (nb_node - 1)) // 2
    nb_edge = np.random.binomial(max_nb_edge, p_link)  # Nb of distincts edges
    edges = []
    seen = set()
    cnt_edge = 0
    list_u = np.random.randint(nb_node, size=nb_edge)
    list_v = np.random.randint(nb_node, size=nb_edge)
    i = 0
    if directed:
        while cnt_edge < nb_edge:
            if i == len(list_u):
                i = 0
                list_u = np.random.randint(nb_node, size=nb_edge - len(edges))
                list_v = np.random.randint(nb_node, size=nb_edge - len(edges))
            u = list_u[i]
            v = list_v[i]
            i += 1
            if u != v and (u, v) not in seen:
                seen.add((u, v))
                edges.append((u, v))
                cnt_edge += 1
    else:
        while cnt_edge < nb_edge:
            if i == len(list_u):
                i = 0
                list_u = np.random.randint(nb_node, size=nb_edge - len(edges))
                list_v = np.random.randint(nb_node, size=nb_edge - len(edges))
            u = list_u[i]
            v = list_v[i]
            i += 1
            if u != v and (u, v) not in seen and (v, u) not in seen:
                seen.add((u, v))
                edges.append((u, v))
                cnt_edge += 1
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
            raise ValueError("The random distribution " + str(presence_law_link) + " is not supported."
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

    return links, link_presence



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
            raise ValueError("The random distribution " + str(presence_law_node) + " is not supported."
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

    return nodes, node_presence


def generate_kcore_example(nb_nodes, k):
    nodes = [i for i in range(nb_nodes)]
    node_presence = [[0, nb_nodes] for i in range(1, nb_nodes + 1)]
    links = [(i, i + 1) for i in range(nb_nodes - 1)]
    link_presence = [[0, i] for i in range(1, len(links) + 1)]
    d = 1 / nb_nodes
    if k > 1:
        for c in range(2, k + 1):
            links += [(i, i + c)
                      if i + c < nb_nodes
                      else (i, (i + c) % (nb_nodes))
                      for i in range(nb_nodes - 1)]
            for i in range(nb_nodes - 1):
                link_presence.append([d / nb_nodes, d])
                d += 1 / nb_nodes

    T = [0, nb_nodes]
    S = sg.stream_graph(times=T,
                        nodes=nodes,
                        node_presence=node_presence,
                        links=links,
                        link_presence=link_presence)
    print("Nodes :", len(nodes), " n_prez :", len(node_presence))
    print("Links :", len(links), " l_prez :", len(link_presence))
    return S


def generate_2core_example(nb_nodes):
    nodes = [i for i in range(nb_nodes)]
    node_presence = [[0, nb_nodes] for i in range(1, nb_nodes + 1)]
    print("Presence :", node_presence)
    links = []
    link_presence = []
    k = 2
    d = 1
    links += [(i, i + k)
              if i + k < nb_nodes
              else (i, (i + k) % (nb_nodes))
              for i in range(nb_nodes)]
    for i in range(nb_nodes):
        link_presence.append([0, d])
        d += 1
    # for l,lp in zip(L,L_presence):
    #     print("link :",l," : ",lp)

    T = [0, nb_nodes]
    S = sg.stream_graph(times=T,
                        nodes=nodes,
                        node_presence=node_presence,
                        links=links,
                        link_presence=link_presence)
    print("Nodes :", len(nodes), " n_prez :", len(node_presence))
    print("Links :", len(links), " l_prez :", len(link_presence))
    return S
