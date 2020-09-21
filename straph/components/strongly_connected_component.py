import copy, pathlib

from joblib import Parallel, delayed

from collections import defaultdict
import msgpack

import gc
from straph.dag import condensation_dag as cdag
from straph.dag import stable_dag as sdag




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
        _, _, u, v = l
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
            new_comp = set([u, v])
            node_to_id_comp[u], node_to_id_comp[v] = len(set_comp), len(set_comp)
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
    with open(storage_file + 'comp.scf', 'rb') as input:
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


# def old_compute_stable_connected_components(S, format="object_with_links", stable_dag=False, isolated_nodes=True,
#                                           streaming_output=None, free_memory=False):
#     '''
#     Compute Stable Connected Components
#     :param S:
#     :param format:
#     :param stable_dag:
#     :return:
#     '''
#     stable_comps = []
#     if stable_dag:
#         scc, c_dag = compute_strongly_connected_components(S, format="object_with_links", condensation_dag=True)
#     else:
#         scc = compute_strongly_connected_components(S, format="object_with_links", condensation_dag=False)
#     for c in scc:
#         # print()
#         # print("comp id:",c.id," comp times :",c.times)
#         # print("comp.nodes :",c.nodes)
#         # print("comp.links :",c.links)
#         stable_comps += c.get_stable_components(format=format)
#
#     if stable_dag:
#         stable_dag = c_dag.get_stable_dag()
#         stable_dag.compute_links_inplace()
#         return stable_comps, stable_dag
#
#     return stable_comps

def compute_strongly_connected_components(S, format="object_with_links", condensation_dag=False, isolated_nodes=True,
                                          streaming_output=None, free_memory=False):
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
    if condensation_dag:
        condensation_dag = cdag.condensation_dag()
    else:
        condensation_dag = None
    #
    id_wcc = S.id

    E = S.ordered_batch_links(free_memory=free_memory)

    if streaming_output:
        opt = open(streaming_output, 'w')
    else:
        opt = None

    for batch in E:
        # print("\n Batch :",batch)
        c = batch[0][0]
        if c == 1:  # ARRIVAL
            #print("node 2 status before :",node_2_status)
            cnt_scc_id = process_batch_link_arrival(batch, node_2_status, tmp_components,
                                                    final_components,
                                                    cnt_scc_id,
                                                    condensation_dag=condensation_dag,
                                                    format=format,
                                                    streaming_output=opt)
            #print("node 2 status after :",node_2_status)

        else:  # DEPARTURE
            #print("node 2 status before :",node_2_status)

            cnt_scc_id = process_batch_link_departure(batch, node_2_status, tmp_components,
                                                      final_components,
                                                      cnt_scc_id,
                                                      condensation_dag=condensation_dag,
                                                      format=format,
                                                      streaming_output=opt)
            #print("node 2 status after :",node_2_status)


    # Add isolated Nodes
    if isolated_nodes:
        for c in S.get_isolated_nodes():
            if format == "cluster":
                final_components.append([c])
                if condensation_dag:
                    condensation_dag.add_node([c])
            elif format == "object" or format == "object_with_links":
                c = strongly_connected_component(id=cnt_scc_id,
                                                 times=[c[0], c[1]],
                                                 nodes={c[2]})
                final_components.append(c)
                if condensation_dag:
                    condensation_dag.add_node(c)
            elif format == "streaming":
                c = (c[0], c[1], 1)
                if streaming_output:
                    opt.write(str(c[0]) + ";" + str(c[1]) + ";" + str(1))
                    opt.write("\n")
                else:
                    final_components.append(c)
            cnt_scc_id += 1

    if condensation_dag:
        condensation_dag.set_id(id_wcc)
        return final_components, condensation_dag
    else:
        return final_components


def process_arrival(l, node_2_status, tmp_components, final_components, condensation_dag,
                    cnt_scc_id, id_wcc, format="cluster", streaming_output=None):
    # ARRIVAL

    u, v = l[2], l[3]
    if u not in node_2_status and v not in node_2_status:
        create_scc(l, node_2_status, tmp_components, format=format)

    elif u in node_2_status and v not in node_2_status:
        cnt_scc_id = update_scc(u, v, l, node_2_status, tmp_components,
                                final_components,
                                condensation_dag=condensation_dag,
                                cnt_scc_id=cnt_scc_id,
                                format=format, streaming_output=streaming_output
                                )

    elif u not in node_2_status and v in node_2_status:
        cnt_scc_id = update_scc(v, u, l, node_2_status, tmp_components,
                                final_components,
                                condensation_dag=condensation_dag,
                                cnt_scc_id=cnt_scc_id,
                                format=format, streaming_output=streaming_output
                                )

    elif node_2_status[u][1] != node_2_status[v][1]:
        cnt_scc_id = merge_scc(l, node_2_status, tmp_components,
                               final_components,
                               condensation_dag=condensation_dag,
                               cnt_scc_id=cnt_scc_id,
                               format=format, streaming_output=streaming_output
                               )
    else:
        node_2_status[u][0] += 1
        node_2_status[v][0] += 1
        current_comp = tmp_components[node_2_status[u][1]]
        current_comp.add_link(l)
    return cnt_scc_id


def merge_scc(l, node_2_status, tmp_components,
              final_components,
              condensation_dag=None,
              cnt_scc_id=None,
              format="cluster", streaming_output=None):
    t0, t1, u, v = l
    id_comp_u = node_2_status[u][1]
    id_comp_v = node_2_status[v][1]

    if len(tmp_components[id_comp_v].nodes) > len(tmp_components[id_comp_u].nodes):
        #  If a component is bigger than another we merge into the bigger one.
        id_comp_u, id_comp_v = id_comp_v, id_comp_u

    comp_1 = tmp_components[id_comp_u]
    comp_2 = tmp_components[id_comp_v]

    if comp_1.times[0] != t0:
        cnt_scc_id = close_component(comp_1, t0, final_components, cnt_scc_id, condensation_dag=condensation_dag,
                                     format=format, streaming_output=streaming_output)
        # predecessor_in_dag[n_comp_1] += [cnt_scc_id - 1]  # previous closed component
        comp_1.set_begin_time(t0)

    if comp_2.times[0] != t0:
        cnt_scc_id = close_component(comp_2, t0, final_components, cnt_scc_id, condensation_dag=condensation_dag,
                                     format=format, streaming_output=streaming_output)
        # predecessor_in_dag[n_comp_1] += [cnt_scc_id - 1]  # previous closed component (n_comp_1 because we merge)

    for n in comp_2.nodes:
        node_2_status[n][1] = id_comp_u  # Actualize referencement before deletion of 2nd comp

    comp_1.merge(comp_2)
    comp_1.add_link(l)  # Add the current link
    node_2_status[u][0] += 1
    node_2_status[v][0] += 1

    tmp_components[id_comp_v] = None
    # print("after merge :",comp_1.links)
    return cnt_scc_id


def update_scc(node_to_update, node_to_add, l, node_2_status, tmp_components,
               final_components,
               condensation_dag=None,
               cnt_scc_id=None,
               format="cluster", streaming_output=None
               ):
    '''

    :param node_to_update:
    :param node_to_add:
    :param l:
    :param node_2_status:
    :param tmp_components:
    :param condensation_dag:
    :param cnt_scc_id:
    :param write_to_msgpack:
    :param id_wcc:
    :return:
    '''
    t0, t1 = l[0], l[1]
    id_current_comp = node_2_status[node_to_update][1]
    current_comp = tmp_components[id_current_comp]
    if current_comp.times[0] != t0:
        cnt_scc_id = close_component(current_comp, t0, final_components, cnt_scc_id,
                                     condensation_dag=condensation_dag,
                                     format=format, streaming_output=streaming_output)
        # predecessor_in_dag[n_current_comp] += [cnt_scc_id - 1]  # previous closed component
        current_comp.set_begin_time(t0)  # Input a new begining time
    current_comp.add_node(node_to_add)  # Add the node to the comp
    current_comp.add_link(l)  # Actualize the component with the new link
    node_2_status[node_to_add] = [1, id_current_comp]
    node_2_status[node_to_update][0] += 1
    # print("after update :",current_comp.links)
    return cnt_scc_id


def create_scc(l, node_2_status, tmp_components, format="cluster"):
    '''
    Create a Strongly Connected Component from the link *l*
    :param l:
    :param node_2_status:
    :param tmp_components:
    :return:
    '''
    t0, t1, u, v = l
    new_id_comp = len(tmp_components)
    node_2_status[u] = [1, new_id_comp]
    node_2_status[v] = [1, new_id_comp]
    if format == "object_with_links":
        lks = [[t0, t1, u, v]]
    else:
        lks = None
    tmp_components.append(strongly_connected_component(times=[t0, t0],
                                                       nodes={u, v},
                                                       set_links={(u, v)},
                                                       links=lks))
    # print("creation :",{u,v})


def process_batch_link_arrival(batch, node_2_status, tmp_components, final_components,
                               cnt_scc_id, condensation_dag=None, format="cluster", streaming_output=None):
    for b in batch:
        t0, t1, u, v = b[1:]
        l = (t0, t1, u, v)
        if u not in node_2_status and v not in node_2_status:
            create_scc(l, node_2_status, tmp_components, format=format)

        elif u in node_2_status and v not in node_2_status:
            cnt_scc_id = update_scc(u, v, l, node_2_status, tmp_components,
                                    final_components,
                                    condensation_dag=condensation_dag,
                                    cnt_scc_id=cnt_scc_id,
                                    format=format, streaming_output=streaming_output
                                    )

        elif u not in node_2_status and v in node_2_status:
            cnt_scc_id = update_scc(v, u, l, node_2_status, tmp_components,
                                    final_components,
                                    condensation_dag=condensation_dag,
                                    cnt_scc_id=cnt_scc_id,
                                    format=format, streaming_output=streaming_output
                                    )

        elif node_2_status[u][1] != node_2_status[v][1]:
            cnt_scc_id = merge_scc(l, node_2_status, tmp_components,
                                   final_components,
                                   condensation_dag=condensation_dag,
                                   cnt_scc_id=cnt_scc_id,
                                   format=format, streaming_output=streaming_output
                                   )
        else:
            node_2_status[u][0] += 1
            node_2_status[v][0] += 1
            current_comp = tmp_components[node_2_status[u][1]]
            current_comp.add_link(l)

    return cnt_scc_id


def process_batch_link_departure(batch, node_2_status, tmp_components,
                                 final_components, cnt_scc_id, condensation_dag=None,
                                 format="cluster",
                                 streaming_output=None):
    '''
    :param batch:
    :param node_2_status:
    :param tmp_components:
    :param final_components:
    :param condensation_dag:
    :param format:
    :param cnt_scc_id:
    :param id_wcc:
    :return:
    '''
    id_comp_to_split = set()
    id_comp_to_close = set()
    nodes_to_remove = set()
    t1 = batch[0][1]
    for l in batch:
        u, v = l[2], l[3]
        node_2_status[u][0] -= 1
        node_2_status[v][0] -= 1
        id_comp = node_2_status[u][1]
        comp = tmp_components[id_comp]
        comp.remove_link((u, v))
        # By default we split the component
        if node_2_status[u][0] == 0 or node_2_status[v][0] == 0:
            #  If it's a node's departure, there is several cases:
            #  1. No more links in the components (it's empty)
            if not comp.set_links:
                cnt_scc_id = close_component(comp, t1, final_components, cnt_scc_id, condensation_dag=condensation_dag,
                                             format=format, streaming_output=streaming_output)
                tmp_components[id_comp] = None
                id_comp_to_split.discard(id_comp)
                id_comp_to_close.discard(id_comp)
                del node_2_status[u]
                del node_2_status[v]
            # 2. A node left but there is still some nodes inside (and other departure to come)
            else:
                if node_2_status[u][0] == 0:
                    id_comp_to_close.add(id_comp)
                    nodes_to_remove.add(u)
                    del node_2_status[u]
                if node_2_status[v][0] == 0:
                    id_comp_to_close.add(id_comp)
                    nodes_to_remove.add(v)
                    del node_2_status[v]
        else:
            id_comp_to_split.add(id_comp)

    for id_comp in id_comp_to_split:
        comp = tmp_components[id_comp]
        # print("before split :",comp.links)
        if comp.set_links:
            R = comp.split()
            if R:
                # We close the current component :)
                cnt_scc_id = close_component(comp, t1, final_components, cnt_scc_id, condensation_dag=condensation_dag,
                                             format=format, streaming_output=streaming_output)
                tmp_components[id_comp] = None
                id_comp_to_close.discard(id_comp)
                for C in R:
                    # New components
                    # assert is_connected(C)
                    C.set_begin_time(t1)  # set new begin time
                    new_id_comp = len(tmp_components)
                    for n in C.nodes:
                        node_2_status[n][1] = new_id_comp
                        # predecessor_in_dag_tmp[len(tmp_components)] += [cnt_scc_id - 1]  # previous closed component
                    tmp_components.append(C)  # to the antecedent of news comp

    # Id comp to close is necessary.
    for id_comp in id_comp_to_close:
        comp = tmp_components[id_comp]
        cnt_scc_id = close_component(comp, t1, final_components, cnt_scc_id, condensation_dag=condensation_dag,
                                     format=format, streaming_output=streaming_output)
        comp.nodes -= nodes_to_remove
        if comp.nodes:
            comp.set_begin_time(t1)  # A node left but other are still presents.
        else:
            raise ValueError("Starfoullah")

    return cnt_scc_id


def close_component(comp,
                    t,
                    final_components,
                    cnt_scc_id,
                    condensation_dag=None,
                    format="cluster",
                    streaming_output=None):
    '''
    Close current component
    :param comp:
    :param n_comp:
    :param t:
    :param final_components:
    :param cnt_scc_id:
    :param condensation_dag:
    :param id_wcc:
    :param format:
    :return:
    '''

    if format == "object" or format == "object_with_links":
        copy_comp = copy.copy(comp)
        copy_comp.set_end_time(t)  # Put an end time to the previous component
        copy_comp.id = cnt_scc_id
        final_components.append(copy_comp)
        c = copy_comp
    elif format == "cluster":
        c = [(comp.times[0], t, n) for n in comp.nodes]
        final_components.append(c)
    elif format == "streaming":
        n_nodes = len(comp.nodes)
        if streaming_output:
            streaming_output.write(str(comp.times[0]) + ";" + str(t) + ";" + str(n_nodes))
            streaming_output.write("\n")
        else:
            c = (comp.times[0], t, n_nodes)
            final_components.append(c)
    if condensation_dag:
        condensation_dag.add_node(c)
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
        rep = "Id SCC :" + str(self.id) + " time window :" + str(self.times) + "\n"
        rep += "Nodes :" + str(self.nodes) + "\n"
        rep += "Links :" + str(self.links) + "\n"
        return rep

    def __copy__(self):
        t = copy.copy(self.times)
        if self.links:
            l = [copy.copy(l) for l in self.links]
        else:
            l = None
        n = copy.copy(self.nodes)
        return strongly_connected_component(times=t,
                                            nodes=n,
                                            links=l)

    def size(self):
        if self.nodes:
            return len(self.nodes)
        else:
            return len(self.get_nodes())

    def to_adjacency_list(self):
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
        if self.links:
            self.links = [l for l in self.links if l[1] >= t and (l[2], l[3]) in self.set_links]

    def set_end_time(self, t):
        self.times[1] = t

    def add_node(self, n):
        self.nodes.add(n)

    def add_link(self, link):
        u, v = link[2], link[3]
        self.set_links.add((u, v))
        if self.links:
            self.links.append(list(link))

    def merge(self, comp):
        self.nodes |= comp.nodes
        self.set_links |= comp.set_links
        if self.links:
            self.links += comp.links

    def remove_link(self, link):
        self.set_links.discard(link)

    def remove_node(self, n):
        self.nodes.discard(n)

    def get_nodes(self):
        return {n for l in self.set_links for n in l}

    def split(self):  # ,link
        # CUSTOM BFS
        R = []

        component_2_set_links = []
        node_2_component = {}
        for l in self.set_links:
            n1, n2 = l
            if n1 not in node_2_component and n2 not in node_2_component:
                n_comp = len(component_2_set_links)
                node_2_component[n1] = n_comp
                node_2_component[n2] = n_comp
                component_2_set_links.append({l})

            elif n1 in node_2_component and n2 not in node_2_component:
                n_comp = node_2_component[n1]
                node_2_component[n2] = n_comp

                component_2_set_links[n_comp].add(l)

            elif n1 not in node_2_component and n2 in node_2_component:
                n_comp = node_2_component[n2]
                node_2_component[n1] = n_comp

                component_2_set_links[n_comp].add(l)

            elif node_2_component[n1] != node_2_component[n2]:
                n_comp_1, n_comp_2 = node_2_component[n1], node_2_component[n2]
                if len(component_2_set_links[n_comp_2]) > len(component_2_set_links[n_comp_1]):
                    #  If a component is bigger than another we merge into the bigger one.
                    n_comp_1, n_comp_2 = n_comp_2, n_comp_1

                for e in component_2_set_links[n_comp_2]:
                    node_2_component[e[0]] = n_comp_1
                    node_2_component[e[1]] = n_comp_1

                component_2_set_links[n_comp_1] |= component_2_set_links[n_comp_2]
                component_2_set_links[n_comp_1].add(l)

                component_2_set_links[n_comp_2] = None

            else:
                component_2_set_links[node_2_component[n1]].add(l)

        # print("Component 2 set links :",component_2_set_links)
        if sum([1 for el in component_2_set_links if el is not None]) > 1:
            for set_links in component_2_set_links:
                if set_links:
                    if self.links:
                        current_links = [l for l in self.links if (l[2], l[3]) in set_links]
                    else:
                        current_links = None
                    # print("  current links :",current_links)
                    # print("  set links :",set_links)
                    current_nodes = {n for l in set_links for n in l}

                    new_c = strongly_connected_component(set_links=set_links,
                                                         nodes=current_nodes,
                                                         links=current_links)
                    R.append(new_c)
        return R


    def get_stable_components(self, format="object"):
        '''
        :return: stable components stemmed from the current component
        '''

        if self.links and len(self.get_interactions_times()) > 1:
            interact_times = self.get_interactions_times()
            time_to_pos = {t: i for t, i in zip(interact_times, range(len(interact_times)))}
            inter_nodes = [set() for _ in range(len(interact_times) - 1)]
            inter_links = [[] for _ in range(len(interact_times) - 1)]
            for l in self.links:
                t0, t1, u, v = l
                t0 = max(t0, self.times[0])
                t1 = min(t1, self.times[1])
                for i in range(time_to_pos[t0], time_to_pos[t1]):
                    inter_nodes[i].add(u)
                    inter_nodes[i].add(v)
                    inter_links[i].append((u, v))
            stable_components = []
            if format == "object" or format == "object_with_links":
                for j in range(len(interact_times) - 1):
                    c = strongly_connected_component(id=(self.id, j),
                                                     times=(interact_times[j], interact_times[j + 1]),
                                                     nodes=set(
                                                         [u for u in inter_nodes[time_to_pos[interact_times[j]]]]),
                                                     links=[l for l in inter_links[time_to_pos[interact_times[j]]]]
                                                     )
                    stable_components.append(c)
            if format == "cluster":
                for j in range(len(interact_times) - 1):
                    c = [(interact_times[j], interact_times[j + 1], u) for u in
                         inter_nodes[time_to_pos[interact_times[j]]]]
                    stable_components.append(c)
        else:
            if format == "object" or format == "object_with_links":
                if self.links:
                    stable_components = [strongly_connected_component(id=self.id,
                                                                      times=self.times,
                                                                      nodes=self.nodes,
                                                                      links=[(u, v) for _, _, u, v in self.links])]
                else:
                    stable_components = [strongly_connected_component(id=self.id,
                                                                      times=self.times,
                                                                      nodes=self.nodes)]
            if format == "cluster":
                stable_components = [[(self.times[0], self.times[1], u) for u in self.nodes]]
        return stable_components

    def get_interactions_times(self):
        interact_times = set()
        for l in self.links:
            t0, t1, _, _ = l
            t0 = max(t0, self.times[0])
            t1 = min(t1, self.times[1])
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
            t0 = max(t0, self.times[0])
            t1 = min(t1, self.times[1])
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
            t0 = max(t0, self.times[0])
            t1 = min(t1, self.times[1])
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
        inter_links = [[] for _ in range(len(interact_times) - 1)]
        for l in self.links:
            t0, t1, u, v = l
            t0 = max(t0, self.times[0])
            t1 = min(t1, self.times[1])
            for i in range(time_2_pos[t0], time_2_pos[t1]):
                inter_links[i].append((u, v))

        def para_clique(i):
            cliques = {}
            t0, t1 = interact_times[i], interact_times[i + 1]
            current_links = inter_links[time_2_pos[t0]]
            if not current_links:
                # Maybe an instant component...
                #print("No current links")
                return {}
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