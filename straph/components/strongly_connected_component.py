import copy, pathlib

from joblib import Parallel, delayed

from collections import defaultdict
import msgpack


from straph.condensation import condensation_dag as cdag




def is_connected(comp):
    '''
    Return true if a component is connected, False otherwise:
    :param comp:
    :return:
    '''
    set_comp = {}
    node_to_id_comp = {}
    for l in comp.links:
        # print("set_comp :",set_comp)
        _,_,u,v = l
        if u in node_to_id_comp and v in node_to_id_comp:
            id_u_comp = node_to_id_comp[u]
            id_v_comp = node_to_id_comp[v]
            if id_u_comp != id_v_comp:
                # MERGE u comp into v comp
                u_comp = set_comp[id_u_comp]
                v_comp = set_comp[id_v_comp]
                u_comp |= v_comp
                for n in v_comp:
                    node_to_id_comp[n] = id_u_comp
                set_comp[id_v_comp] = None
        elif u in node_to_id_comp and v not in node_to_id_comp:
            id_u_comp = node_to_id_comp[u]
            u_comp = set_comp[id_u_comp]
            node_to_id_comp[v] = id_u_comp
            u_comp.add(v)

        elif u not in node_to_id_comp and v in node_to_id_comp:
            id_v_comp = node_to_id_comp[v]
            v_comp = set_comp[id_v_comp]
            node_to_id_comp[u] = id_v_comp
            v_comp.add(u)
        else:
            new_comp = set([u,v])
            node_to_id_comp[u],node_to_id_comp[v] = len(set_comp),len(set_comp)
            set_comp[len(set_comp)] = new_comp
    if sum([1 for c in set_comp if set_comp[c] is not None]) == 1:
        return True
    return False





def store_scc_to_sgf(scc, storage_path):
    '''
    Dump the strongly connected component to a msgpack_file (sgf format)
    :param output: Writable bytes file open(...,'wb')
    :return:
    '''
    # print("STORE COMP")
    pathlib.Path(storage_path).mkdir(parents=True, exist_ok=True)
    packer = msgpack.Packer(use_bin_type=True)
    links = []
    links_to_sort = []
    begin_time, end_time = scc.times
    if not scc.links:
        print("Empty SCC :", storage_path, " do not store !")
    else:
        for l, t0, t1 in scc.links:
            if t0 == begin_time and t1 == end_time:
                links.append(l)
            else:
                links_to_sort.append((l, t0, t1))
        links_to_sort = sorted(links_to_sort, key=lambda x: (x[1], -x[0]))  # Sort the link
        links += links_to_sort
        with open(storage_path + "scc_" + str(scc.id) + '.sccf', 'wb') as output:
            output.write(packer.pack((begin_time, end_time)))
            for l in links:
                output.write(packer.pack(l))


def write_scc_to_msgpack(scc, id_wcc, output_file, packer=None):
    '''
    Dump the strongly connected component to a msgpack_file (sgf format)
    :param output: Writable bytes file open(...,'wb')
    :return:
    '''
    # print("STORE COMP")
    links = []
    links_to_sort = []
    if not scc.links:
        print("Empty SCC do not store !")
        return
    begin_time, end_time = scc.times
    for t0, t1, u, v in scc.links:
        if t0 == begin_time and t1 == end_time:
            links.append((u, v))
        else:
            links_to_sort.append((t0, t1, u, v))
    links_to_sort = sorted(links_to_sort, key=lambda x: x[0])  # Sort the link
    links += links_to_sort
    if not packer:
        packer = msgpack.Packer(use_bin_type=True)
    output_file.write(packer.pack((id_wcc, scc.id, begin_time, end_time, links)))


def neighborhood_and_degrees_from_links(t0, t1, links):
    a_l = defaultdict(set)
    for l in links:
        u = (t0, t1, l[0])
        v = (t0, t1, l[1])
        a_l[u].add(v)
        a_l[v].add(u)
    degrees = {n: len(a_l[n]) for n in a_l}
    return a_l, degrees


def algo_kcores_batagelj(a_l, degrees, core_ordering=False):
    '''
    Compute k_cores of a static graph from its adjacency list and nodes degrees
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
        u = (t0, t1, l[0])
        v = (t0, t1, l[1])
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
                for v in a_l[u]:
                    if v in node_label:
                        if node_label[v] == k:
                            final_C = copy.copy(C) + [u, v]
                            R.append(final_C)

    else:
        for u in node_label:
            # INDUCED NODES
            new_node_label = defaultdict(int)
            if a_l[u]:
                for v in a_l[u]:
                    if v in node_label:
                        if node_label[v] == k:
                            new_node_label[v] = k - 1
                new_C = copy.copy(C) + [u]
                if new_node_label:
                    algo_kcliques_KCList(k - 1, a_l, new_node_label, C=new_C, R=R)
    return R


def catch_scc_from_msgpack(id_wcc, id_scc, storage_file):
    with open(storage_file + 'scc.scf', 'rb') as input:
        unpacker = msgpack.Unpacker(input)
        for l in unpacker:
            if l[0] == (id_wcc, id_scc):
                return strongly_connected_component(id=id_scc, times=l[1],
                                                    links=l[2])


def load_scc_from_sgf(storage_path, scc_id):
    '''
    Load a weakly connected component from a msgpack_file (sgf format).
    :param storage_path: Readable bytes file open(...,'rb')
    :return:
    '''
    unpacker = msgpack.Unpacker()
    links = []
    with open(storage_path + "scc_" + str(scc_id) + '.sgf', 'rb') as input:
        t0, t1 = unpacker.unpack(input.readline())
        for i in input:
            l = unpacker.unpack(i)
            if len(l) <= 2:
                links.append((l, t0, t1))
            else:
                links.append(l)
    return strongly_connected_component(id=scc_id,
                                        times=(t0, t1),
                                        links=links)


def compute_stable_connected_components(S, format="cluster", stable_dag=False):
    '''
    Compute Stable Connected Components
    :param S:
    :param format:
    :param stable_dag:
    :return:
    '''
    stable_comps = []
    scc, c_dag = compute_strongly_connected_components(S, format="object", condensation_dag=True)

    for c in scc:
        # print()
        # print("scc id:",c.id," scc times :",c.times)
        # print("scc.nodes :",c.nodes)
        # print("scc.links :",c.links)
        stable_comps += c.get_stable_components(format=format)

    if stable_dag:
        stable_dag = c_dag.get_stable_dag()
        stable_dag.compute_links_inplace()
        return stable_comps, stable_dag

    return stable_comps


def compute_strongly_connected_components(S, format="cluster", condensation_dag=False):
    '''
    Compute Strongly Connectec Components (SCC) of a Stream Graph.
    :param S: A Stream Graph
    :param format: Format of the output can be "cluster" or "scc_object"
    :param condensation_dag: Boolean, true if we want to output the Condensation DAG, false otherwise
    :return:
    '''
    node_2_status = {}  # Dictionary associating a node to his current status : (current degree, number current comp)
    tmp_components = []  # List of current strongly connected components (object)
    final_components = []  # Clusters : [(t0,t1,u)]...
    cnt_scc_id = 0
    # Condensation DAG
    scc_dag = cdag.scc_dag()
    predecessor_in_dag_tmp = defaultdict(list)
    predecessor_in_dag_final = {}
    #
    id_wcc = S.id
    E = S.ordered_links()
    t_last_departure = None
    batch_departure = []
    for i in E:
        # print("i :",i)
        c= i[0]
        l = i[1:]
        if c == 1: # ARRIVAL
            if batch_departure:
                cnt_scc_id = new_departure_procedure(batch_departure, node_2_status, tmp_components,
                                                     final_components, scc_dag,
                                                     format, cnt_scc_id, id_wcc,
                                                     predecessor_in_dag_tmp,
                                                     predecessor_in_dag_final)
                batch_departure = []

            cnt_scc_id = process_arrival(l, node_2_status, tmp_components,
                                         final_components, scc_dag,
                                         format, cnt_scc_id, id_wcc,
                                         predecessor_in_dag_tmp, predecessor_in_dag_final)
        else: # DEPARTURE
            t = l[0]
            if t == t_last_departure:
                batch_departure.append(l)
            else:
                if batch_departure:
                    cnt_scc_id = new_departure_procedure(batch_departure, node_2_status, tmp_components,
                                                         final_components, scc_dag,
                                                         format, cnt_scc_id, id_wcc,
                                                         predecessor_in_dag_tmp,
                                                         predecessor_in_dag_final)
                batch_departure = [l]
                t_last_departure = t

    # Process Last links :
    if batch_departure:
        cnt_scc_id = new_departure_procedure(batch_departure, node_2_status, tmp_components,
                                             final_components, scc_dag,
                                             format, cnt_scc_id, id_wcc,
                                             predecessor_in_dag_tmp,
                                             predecessor_in_dag_final)


    # Add isolated Nodes
    for c in S.get_isolated_nodes():
        isolated_scc = strongly_connected_component(id=cnt_scc_id,
                                                    times=[c[0], c[1]],
                                                    nodes=set([c[2]]))
        if format == "cluster":
            final_components.append([c])
        if format == "object":
            final_components.append(isolated_scc)
        if condensation_dag:
            scc_dag.add_node(isolated_scc)
        cnt_scc_id += 1

    if condensation_dag:
        scc_dag.set_id(id_wcc)
        return final_components, scc_dag
    else:
        return final_components


def process_arrival(l, node_2_status, tmp_components, final_components, scc_dag,
                    format, cnt_scc_id, id_wcc, predecessor_in_dag_tmp, predecessor_in_dag_final):
    # ARRIVAL

    u,v = l[2],l[3]
    if u not in node_2_status and v not in node_2_status:
        create_scc(l, node_2_status, tmp_components, predecessor_in_dag_tmp)

    elif u in node_2_status and v not in node_2_status:
        cnt_scc_id = update_scc(u, v, l, node_2_status, tmp_components,
                                final_components,
                                scc_dag,
                                predecessor_in_dag_tmp,
                                predecessor_in_dag_final,
                                cnt_scc_id=cnt_scc_id,
                                id_wcc=id_wcc,
                                format=format
                                )

    elif u not in node_2_status and v in node_2_status:
        cnt_scc_id = update_scc(v, u, l, node_2_status, tmp_components,
                                final_components,
                                scc_dag,
                                predecessor_in_dag_tmp,
                                predecessor_in_dag_final,
                                cnt_scc_id=cnt_scc_id,
                                id_wcc=id_wcc,
                                format=format
                                )

    elif node_2_status[u][1] != node_2_status[v][1]:
        cnt_scc_id = merge_scc(l, node_2_status, tmp_components,
                               final_components,
                               scc_dag,
                               predecessor_in_dag_tmp,
                               predecessor_in_dag_final,
                               cnt_scc_id=cnt_scc_id,
                               id_wcc=id_wcc,
                               format=format
                               )
    else:
        node_2_status[u][0] += 1
        node_2_status[v][0] += 1
        current_comp = tmp_components[node_2_status[u][1]]
        current_comp.add_link(l)
    return cnt_scc_id




def merge_scc(l, node_2_status, tmp_components,
              final_components,
              scc_dag,
              predecessor_in_dag,
              predecessor_in_dag_final,
              cnt_scc_id=None,
              id_wcc=None,
              format="cluster"):
    t0, t1, u, v = l
    n_comp_1 = node_2_status[u][1]
    n_comp_2 = node_2_status[v][1]
    comp_1 = tmp_components[n_comp_1]
    comp_2 = tmp_components[n_comp_2]
    if comp_1.times[0] != t0:
        cnt_scc_id = close_component(comp_1, n_comp_1, t0, final_components, cnt_scc_id, scc_dag,
                                     predecessor_in_dag,
                                     predecessor_in_dag_final,
                                     id_wcc=id_wcc,
                                     format=format)
        # predecessor_in_dag[n_comp_1] += [cnt_scc_id - 1]  # previous closed component
        comp_1.set_begin_time(t0)

    if comp_2.times[0] != t0:
        cnt_scc_id = close_component(comp_2, n_comp_2, t0, final_components, cnt_scc_id, scc_dag,
                                     predecessor_in_dag,
                                     predecessor_in_dag_final,
                                     id_wcc=id_wcc,
                                     format=format)
        # predecessor_in_dag[n_comp_1] += [cnt_scc_id - 1]  # previous closed component (n_comp_1 because we merge)

    comp_1.merge(comp_2)
    comp_1.add_link(l)  # Add the current link

    for n in comp_2.get_nodes():
        node_2_status[n][1] = n_comp_1  # Actualize referencement before deletion of 2nd comp

    node_2_status[u][0] += 1
    node_2_status[v][0] += 1

    tmp_components[n_comp_2] = None
    # print("after merge :",comp_1.links)
    return cnt_scc_id


def update_scc(node_to_update, node_to_add, l, node_2_status, tmp_components,
               final_components,
               scc_dag,
               predecessor_in_dag,
               predecessor_in_dag_final,
               cnt_scc_id=None,
               id_wcc=None,
               format="cluster"
               ):
    '''

    :param node_to_update:
    :param node_to_add:
    :param l:
    :param node_2_status:
    :param tmp_components:
    :param scc_dag:
    :param cnt_scc_id:
    :param write_to_msgpack:
    :param id_wcc:
    :return:
    '''
    t0, t1 = l[0], l[1]
    n_current_comp = node_2_status[node_to_update][1]
    current_comp = tmp_components[n_current_comp]
    if current_comp.times[0] != t0:
        cnt_scc_id = close_component(current_comp, n_current_comp,
                                     t0, final_components, cnt_scc_id, scc_dag,
                                     predecessor_in_dag,
                                     predecessor_in_dag_final,
                                     id_wcc=id_wcc,
                                     format=format)
        # predecessor_in_dag[n_current_comp] += [cnt_scc_id - 1]  # previous closed component
        current_comp.set_begin_time(t0)  # Input a new begining time
    current_comp.add_link(l)  # Actualize the component (will add the node and the link)
    node_2_status[node_to_add] = [1, n_current_comp]
    node_2_status[node_to_update][0] += 1
    # print("after update :",current_comp.links)
    return cnt_scc_id


def create_scc(l, node_2_status, tmp_components, predecessor_in_dag):
    '''
    Create a Strongly Connected Component from the link *l*
    :param l:
    :param node_2_status:
    :param tmp_components:
    :return:
    '''
    t0, t1, u, v = l
    n_comp = len(tmp_components)
    node_2_status[u] = [1, n_comp]
    node_2_status[v] = [1, n_comp]
    tmp_components.append(strongly_connected_component(times=[t0, t0],
                                                       set_links=set([(u, v)]),
                                                       links=[[t0, t1, u, v]]))
    # print("creation :",{u,v})


# def process_batch_departure(batch_departure, node_2_status, tmp_components, final_components, scc_dag,
#                             format, cnt_scc_id, id_wcc, predecessor_in_dag):
#     # DEPARTURE
#     for l in batch_departure:
#         u, v = l[1], l[2]
#         node_2_status[u][0] -= 1
#         node_2_status[v][0] -= 1
#         # # Regrouper les liens partant concernant la même composantes
#         if node_2_status[u][0] == 0 and node_2_status[v][0] == 0:
#             cnt_scc_id = double_departure(l, node_2_status, tmp_components,
#                                           final_components,
#                                           scc_dag,
#                                           predecessor_in_dag,
#                                           cnt_scc_id=cnt_scc_id,
#                                           id_wcc=id_wcc,
#                                           format=format
#                                           )
#         elif node_2_status[u][0] == 0:
#             cnt_scc_id = single_departure(u, l, node_2_status, tmp_components,
#                                           final_components,
#                                           scc_dag,
#                                           predecessor_in_dag,
#                                           cnt_scc_id=cnt_scc_id,
#                                           id_wcc=id_wcc,
#                                           format=format
#                                           )
#         elif node_2_status[v][0] == 0:
#             cnt_scc_id = single_departure(v, l, node_2_status, tmp_components,
#                                           final_components,
#                                           scc_dag,
#                                           predecessor_in_dag,
#                                           cnt_scc_id=cnt_scc_id,
#                                           id_wcc=id_wcc,
#                                           format=format
#                                           )
#         else:
#             # In this procedure we can use the BFS to compute any kind of distance, which can be of use further on.
#             cnt_scc_id = split_procedure(l, node_2_status, tmp_components,
#                                          final_components,
#                                          scc_dag,
#                                          predecessor_in_dag,
#                                          cnt_scc_id=cnt_scc_id,
#                                          id_wcc=id_wcc,
#                                          format=format)
#     return cnt_scc_id
#

def new_departure_procedure(batch_departure, node_2_status, tmp_components,
                            final_components, scc_dag,
                            format, cnt_scc_id, id_wcc,
                            predecessor_in_dag_tmp, predecessor_in_dag_final):
    '''

    :param batch_departure:
    :param node_2_status:
    :param tmp_components:
    :param final_components:
    :param scc_dag:
    :param format:
    :param cnt_scc_id:
    :param id_wcc:
    :param predecessor_in_dag_tmp:
    :param predecessor_in_dag_final:
    :return:
    '''
    id_comp_to_split = set()
    id_comp_to_close = set()
    t1 = None
    for l in batch_departure:
        t1, u, v = l
        node_2_status[u][0] -= 1
        node_2_status[v][0] -= 1
        n_comp = node_2_status[u][1]
        comp = tmp_components[n_comp]
        comp.remove_link((u,v))
        if node_2_status[u][0] == 0 or node_2_status[v][0] == 0:
            id_comp_to_close.add(n_comp)
            if node_2_status[u][0] == 0:
                del node_2_status[u]
            if node_2_status[v][0] == 0:
                del node_2_status[v]
            if len(comp.set_links) == 0:
                cnt_scc_id = close_component(comp, n_comp, t1, final_components, cnt_scc_id, scc_dag,
                                             predecessor_in_dag_tmp,
                                             predecessor_in_dag_final,
                                             id_wcc=id_wcc,
                                             format=format)
                tmp_components[n_comp] = None
                id_comp_to_split.discard(n_comp)
                id_comp_to_close.discard(n_comp)
        else:
            id_comp_to_split.add(n_comp)

    for n_comp in id_comp_to_split:
        comp = tmp_components[n_comp]
        # print("before split :",comp.links)
        R = comp.split()
        if R :
            id_comp_to_close.discard(n_comp)
            # We close the current component :)
            cnt_scc_id = close_component(comp, n_comp, t1, final_components, cnt_scc_id, scc_dag,
                                         predecessor_in_dag_tmp,
                                         predecessor_in_dag_final,
                                         id_wcc=id_wcc,
                                         format=format)
            tmp_components[n_comp] = None
            for C in R:
                # New components
                comp_nodes = C.get_nodes()
                # assert is_connected(C)
                C.set_begin_time(t1)  # set new begin time
                for n in comp_nodes:
                    node_2_status[n][1] = len(tmp_components)
                    # predecessor_in_dag_tmp[len(tmp_components)] += [cnt_scc_id - 1]  # previous closed component
                tmp_components.append(C)  # to the antecedent of news comp

    # Id comp to close is necessary.
    for n_comp in id_comp_to_close:
        comp = tmp_components[n_comp]
        cnt_scc_id = close_component(comp, n_comp, t1, final_components, cnt_scc_id, scc_dag,
                                     predecessor_in_dag_tmp,
                                     predecessor_in_dag_final,
                                     id_wcc=id_wcc,
                                     format=format)
        comp.set_begin_time(t1)


    return cnt_scc_id


# def single_departure(node_to_remove, l, node_2_status, tmp_components,
#                      final_components,
#                      scc_dag,
#                      predecessor_in_dag,
#                      cnt_scc_id=None,
#                      id_wcc=None,
#                      format="cluster"):
#     # Remove the link, close the component
#     t1 = l[0]
#     n_comp_to_keep = node_2_status[node_to_remove][1]
#     # print(" SINGLE REMOVAL comp : ", n_comp_to_keep, " REMOVE : ", node_to_remove)
#     comp_to_keep = tmp_components[n_comp_to_keep]
#     comp_to_keep.remove_link((l[1], l[2]))  # Remove Link
#     # if comp_to_keep.times[0] != t1:
#     cnt_scc_id = close_component(comp_to_keep, n_comp_to_keep, t1, final_components, cnt_scc_id, scc_dag,
#                                  predecessor_in_dag,
#                                  id_wcc=id_wcc,
#                                  format=format)
#     predecessor_in_dag[n_comp_to_keep] += [cnt_scc_id - 1]
#     comp_to_keep.set_begin_time(t1)  # Set begin time
#     del node_2_status[node_to_remove]
#     return cnt_scc_id
#
#
# def double_departure(l, node_2_status, tmp_components,
#                      final_components,
#                      scc_dag,
#                      predecessor_in_dag,
#                      cnt_scc_id=None,
#                      id_wcc=None,
#                      format="cluster"):
#     t1, u, v = l
#     n_current_comp = node_2_status[u][1]
#     # print(" DOUBLE REMOVAL comp : ", n_current_comp, " u :", u, " v :", v)
#     current_comp = tmp_components[n_current_comp]
#     current_comp.remove_link((u, v))
#     if len(current_comp.set_links) == 0:  # If there is less than 2 nodes left, clear the component
#         # if current_comp.times[0] != t1:
#         cnt_scc_id = close_component(current_comp, n_current_comp, t1, final_components, cnt_scc_id, scc_dag,
#                                      predecessor_in_dag,
#                                      id_wcc=id_wcc,
#                                      format=format)
#         # End of the SCC, no successor
#         tmp_components[n_current_comp] = None
#         del node_2_status[u], node_2_status[v]
#     else:  # Else double departure means a split
#         cnt_scc_id = split_procedure(l, node_2_status, tmp_components,
#                                      scc_dag,
#                                      predecessor_in_dag,
#                                      cnt_scc_id=cnt_scc_id,
#                                      id_wcc=id_wcc,
#                                      format=format)
#
#     return cnt_scc_id
#
#
# def split_procedure(l, node_2_status, tmp_components,
#                     final_components,
#                     scc_dag,
#                     predecessor_in_dag,
#                     cnt_scc_id=None,
#                     id_wcc=None,
#                     format="cluster"):
#     t1, u, v = l
#     n_current_comp = node_2_status[u][1]
#     current_comp = tmp_components[n_current_comp]
#     current_comp.remove_link((u, v))
#     R = current_comp.split((u, v))
#     # print(" SPLIT COMP : ", n_current_comp)
#     #######################
#     # Splitting procedure #
#     #######################
#     if R:
#         # print(" SPLIT SUCCESFULL")
#         # if current_comp.times[0] != t1:
#         cnt_scc_id = close_component(current_comp, n_current_comp,
#                                      t1, final_components, cnt_scc_id, scc_dag, predecessor_in_dag,
#                                      id_wcc=id_wcc,
#                                      format=format)
#
#         u_is_present = False
#         v_is_present = False
#         for C in R:
#             C.set_begin_time(t1)  # set new begin time
#             # print("  NEW COMP : ", len(tmp_components))
#             C_nodes = C.get_nodes()
#             # print(" C LINKS :",C.set_links)
#             # print(" C NODES :", C_nodes)
#             for n in C_nodes:
#                 node_2_status[n][1] = len(tmp_components)
#             if u in C_nodes:
#                 u_is_present = True
#             if v in C_nodes:
#                 v_is_present = True
#             predecessor_in_dag[len(tmp_components)] += [cnt_scc_id - 1]  # previous closed component
#             tmp_components.append(C)  # to the antecedent of news comp
#         if not u_is_present:
#             # print("Remove u :", u)
#             del node_2_status[u]
#         if not v_is_present:
#             # print("Remove v :", v)
#             del node_2_status[v]
#         tmp_components[n_current_comp] = None
#     # else:
#     #     print(" NO SPLIT : REPLACEMENT LINK FOUND")
#     return cnt_scc_id
#
#


def close_component(comp,
                    n_comp,
                    t0,
                    final_components,
                    cnt_scc_id,
                    scc_dag,
                    predecessor_in_dag,
                    predecessor_in_dag_final,
                    id_wcc=None,
                    format="cluster"):
    '''
    Close current component
    :param comp:
    :param n_comp:
    :param t0:
    :param final_components:
    :param cnt_scc_id:
    :param scc_dag:
    :param predecessor_in_dag:
    :param predecessor_in_dag_final:
    :param id_wcc:
    :param format:
    :return:
    '''
    copy_comp = copy.copy(comp)
    copy_comp.set_end_time(t0)  # Put an end time to the previous component
    copy_comp.id = cnt_scc_id
    #
    set_nodes = set()
    for l in copy_comp.links:
        _, _, u, v = l
        set_nodes.add(u)
        set_nodes.add(v)
    copy_comp.nodes = set_nodes
    #
    if format == "object":
        final_components.append(copy_comp)
    if format == "cluster":
        c_nodes = set_nodes
        c = [(copy_comp.times[0], copy_comp.times[1], n) for n in c_nodes]
        final_components.append(c)

    # predecessor_in_dag_final[cnt_scc_id] = copy.copy(predecessor_in_dag[n_comp])
    # predecessor_in_dag[n_comp] = []
    if scc_dag:
        scc_dag.add_node(copy_comp)
        cnt_scc_id += 1
    return cnt_scc_id


class strongly_connected_component:
    def __init__(self,
                 id=None,
                 times=None,
                 nodes=None,
                 set_links=None,
                 links=None
                 ):
        '''
        A basic constructor for a connected component object
        :param id : identifier of the SCC STREAM (a string)
        :param times = [beginning time, ending time]
        :param nodes : A set of nodes present in the component
        :param set_links : A set of links present in the component (Only useful during construction)
        :param links : a list of 'segmented' links
        '''
        self.id = id
        self.times = times
        self.nodes = nodes
        self.set_links = set_links
        self.links = links

    def __repr__(self):
        rep = "Id SCC :" + str(self.id) + " time window :" + str(self.times)
        rep += "\nNodes :" + str(self.nodes)
        rep += "\nLinks :" + str(self.links)
        return rep

    def __copy__(self):
        t = copy.copy(self.times)
        l = [copy.copy(l) for l in self.links]
        return strongly_connected_component(times=t,
                                            links=l)

    def size(self):
        return len(self.nodes)

    def to_al(self):
        al = defaultdict(set)
        if self.links:
            for l in self.links:
                # print("l :",l)
                u, v = l
                al[u].add(v)
                al[v].add(u)
        return al

    def set_begin_time(self, t):
        self.times = [t, t]
        new_links = []
        for t0, t1, u, v in self.links:
            if t1 >= t and (u, v) in self.set_links:
                new_links.append([t0, t1, u, v])
        self.links = new_links

    def set_end_time(self, t):
        self.times[1] = t
        # for l in self.links:
        #     if l[1] > t:
        #         l[1] = t

    def add_link(self, link):
        u, v = link[2], link[3]
        self.set_links.add((u, v))
        self.links.append(list(link))

    def merge(self, comp):
        self.set_links |= comp.set_links
        self.links += comp.links

    def remove_link(self, link):
        self.set_links.discard(link)

    def get_nodes(self):
        return set([n for l in self.set_links for n in l])

    def split(self): #,link
        # CUSTOM BFS
        component_2_set_links = []
        node_2_component = {}
        for l in self.set_links:
            n1, n2 = l
            if n1 not in node_2_component and n2 not in node_2_component:
                node_2_component[n1] = len(component_2_set_links)
                node_2_component[n2] = len(component_2_set_links)
                component_2_set_links.append(set([l]))

            elif n1 in node_2_component and n2 not in node_2_component:
                n_comp = node_2_component[n1]
                component_2_set_links[n_comp].add(l)
                node_2_component[n2] = n_comp

            elif n1 not in node_2_component and n2 in node_2_component:
                n_comp = node_2_component[n2]
                component_2_set_links[n_comp].add(l)
                node_2_component[n1] = n_comp

            elif node_2_component[n1] != node_2_component[n2]:
                n_comp_1, n_comp_2 = node_2_component[n1], node_2_component[n2]
                component_2_set_links[n_comp_1] |= component_2_set_links[n_comp_2]
                component_2_set_links[n_comp_1].add(l)
                for e in component_2_set_links[n_comp_2]:
                    node_2_component[e[0]] = n_comp_1
                    node_2_component[e[1]] = n_comp_1
                component_2_set_links[n_comp_2] = None
            else:
                component_2_set_links[node_2_component[n1]].add(l)

        R = []
        # print("Component 2 set links :",component_2_set_links)
        if sum([1 for i in component_2_set_links if i is not None]) > 1:
            for set_links in component_2_set_links:
                if set_links is not None:
                    current_links = [l for l in self.links if (l[2], l[3]) in set_links]
                    # print("  current links :",current_links)
                    # print("  set links :",set_links)
                    R.append(
                        strongly_connected_component(set_links=set_links,
                                                     links=current_links))
        return R

    def get_stable_components(self, format="object"):
        '''
        :return: stable components stemmed from the current component
        '''

        if self.links and len(self.get_interactions_times()) > 1:
            interact_times = self.get_interactions_times()
            time_2_pos = {t: i for t, i in zip(interact_times, range(len(interact_times)))}
            inter_nodes = [set() for k in range(len(interact_times) - 1)]
            inter_links = [[] for k in range(len(interact_times) - 1)]
            for l in self.links:
                t0, t1, u, v = l
                t0 = max(t0, self.times[0])
                t1 = min(t1, self.times[1])
                for i in range(time_2_pos[t0], time_2_pos[t1]):
                    inter_nodes[i].add(u)
                    inter_nodes[i].add(v)
                    inter_links[i].append((u, v))
            stable_components = []
            if format == "object":
                for j in range(len(interact_times) - 1):
                    c = strongly_connected_component(id=(self.id, j),
                                                     times=(interact_times[j], interact_times[j + 1]),
                                                     nodes=set([u for u in inter_nodes[time_2_pos[interact_times[j]]]]),
                                                     links=[l for l in inter_links[time_2_pos[interact_times[j]]]]
                                                     )
                    stable_components.append(c)
            if format == "cluster":
                for j in range(len(interact_times) - 1):
                    c = [(interact_times[j], interact_times[j + 1], u) for u in
                         inter_nodes[time_2_pos[interact_times[j]]]]
                    stable_components.append(c)
        else:
            if format == "object":
                stable_components = [strongly_connected_component(id=self.id,
                                                                  times=self.times,
                                                                  nodes=self.nodes)]
            if format == "cluster":
                stable_components = [[(self.times[0], self.times[1], u) for u in self.nodes]]
        return stable_components

    def get_interactions_times(self):
        interact_times = set()
        for l in self.links:
            t0,t1,_,_ = l
            t0 = max(t0,self.times[0])
            t1 = min(t1,self.times[1])
            interact_times.add(t0)
            interact_times.add(t1)
        return sorted(interact_times)

    def core_number(self):
        # Get interaction times
        # Divides into slices
        # Compute cores with Batagelj (can be parallelized)
        L = defaultdict(list)
        if self.size() == 2:
            L[1] = [(self.times[0], self.times[1], list(self.nodes)[0])]
            L[1].append((self.times[0], self.times[1], list(self.nodes)[1]))
            return L
        interact_times = self.get_interactions_times()
        time_2_pos = {t: i for t, i in zip(interact_times, range(len(interact_times)))}
        inter_links = [[] for _ in range(len(interact_times) - 1)]
        for l in self.links:
            t0, t1, u, v = l
            t0 = max(t0,self.times[0])
            t1 = min(t1,self.times[1])
            for i in range(time_2_pos[t0], time_2_pos[t1]):
                inter_links[i].append((u, v))

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
        return L

    def k_core(self, k):
        L = self.core_number()
        return L[k]

    def k_clique(self, k):
        # Get interaction times
        # Divides into slices
        # Compute cores with KCList (can be parallelized)
        L = []
        if self.size() == 2:
            # L[1] = [(list(self.nodes)[0],self.times[0], self.times[1] )]
            # L[1].append((list(self.nodes)[1],self.times[0], self.times[1]))
            return L
        interact_times = self.get_interactions_times()
        time_2_pos = {t: i for t, i in zip(interact_times, range(len(interact_times)))}
        inter_links = [[] for k in range(len(interact_times) - 1)]
        for l in self.links:
            t0, t1, u, v = l
            t0 = max(t0,self.times[0])
            t1 = min(t1,self.times[1])
            for i in range(time_2_pos[t0], time_2_pos[t1]):
                inter_links[i].append((u, v))

        def para_clique(i, k):
            cliques = []
            t0, t1 = interact_times[i], interact_times[i + 1]
            current_links = inter_links[time_2_pos[t0]]
            a_l, degrees = neighborhood_and_degrees_from_links(t0, t1, current_links)
            cores, core_ordering = algo_kcores_batagelj(a_l, degrees, core_ordering=True)

            a_l = get_graph_from_ordering(t0, t1, current_links, core_ordering)

            node_label = defaultdict(int, {n: k for n in degrees})
            cliques += algo_kcliques_KCList(k, a_l, node_label, R=[])
            return cliques

        r = Parallel(n_jobs=1)(delayed(para_clique)(i, k) for i in range(len(interact_times) - 1))
        for l in r:
            L += l
        return L

    def all_cliques(self):
        # Get interaction times
        # Divides into slices
        # Compute cores with KCList (can be parallelized)
        L = defaultdict(list)
        if self.size() == 2:
            # L[1] = [(list(self.nodes)[0],self.times[0], self.times[1] )]
            # L[1].append((list(self.nodes)[1],self.times[0], self.times[1]))
            return L
        interact_times = self.get_interactions_times()
        time_2_pos = {t: i for t, i in zip(interact_times, range(len(interact_times)))}
        inter_links = [[] for k in range(len(interact_times) - 1)]
        for l in self.links:
            t0, t1, u, v = l
            t0 = max(t0,self.times[0])
            t1 = min(t1,self.times[1])
            for i in range(time_2_pos[t0], time_2_pos[t1]):
                inter_links[i].append((u, v))
        def para_clique(i):
            cliques = {}
            t0, t1 = interact_times[i], interact_times[i + 1]
            current_links = inter_links[time_2_pos[t0]]
            a_l, degrees = neighborhood_and_degrees_from_links(t0, t1, current_links)
            cores, core_ordering = algo_kcores_batagelj(a_l, degrees, core_ordering=True)
            max_core_number = max(cores.values())
            # print("\n t0,t1 :",t0,t1)
            # print("Max core number :",max_core_number)
            # print("Core ordering :",core_ordering)
            K = 3
            a_l = get_graph_from_ordering(t0, t1, current_links, core_ordering)
            # print(" a_l :",a_l)
            # print("degrees :",degrees)
            while K <= max_core_number + 1:
                node_label = defaultdict(int, {n: K for n in degrees})
                cliques[K] = algo_kcliques_KCList(K, a_l, node_label, R=[])
                K += 1
            # print("kcliques :",cliques)
            return cliques

        r = Parallel(n_jobs=1)(delayed(para_clique)(i) for i in range(len(interact_times) - 1))
        for l in r:
            for k, v in l.items():
                L[k] += v

        return L

# TODO : Faire une classe "component" dont strongly_connected_component et weakly_connected_component hériteront.
