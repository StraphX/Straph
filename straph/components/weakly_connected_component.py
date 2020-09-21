import copy
import msgpack

from straph import stream as sg
import matplotlib.pyplot as plt


def store_wcc_to_sgf(wcc, output_file):
    '''
    Dump the weakly connected component to a msgpack_file (scf format)
    :param output_file: Writable bytes file open(...,'wb')
    :return:
    '''
    links = []
    begin_time = min([l[0] for l in wcc.links])
    end_time = max([l[1] for l in wcc.links])
    for t0, t1,u,v in wcc.links:
        links.append((1, t0, t1,u,v))  # code each link, 1 for a beginning, -1 for an ending
        links.append((-1, t1,u,v))
    links = sorted(links, key=lambda x: (x[1],-x[0]))  # Sort the links
    packer = msgpack.Packer()
    output_file.write(packer.pack((wcc.id,begin_time,end_time,links)))


class weakly_connected_component:
    def __init__(self,
                 id=None,
                 nodes=None,
                 links=None,
                 active_links = None,
                 times=None
                 ):
        '''

        :param id: identifier of the WCC
        :param nodes: set of 'segmented' nodes in WCC
        :param links: list of 'segmented' link inside WCC
        '''
        self.id = id
        self.nodes = nodes
        self.active_links = active_links
        self.links = links
        self.times = times

    def __copy__(self):
        return weakly_connected_component(nodes=copy.copy(self.nodes),
                                          active_links=copy.copy(self.active_links),
                                          links=copy.copy(self.links))

    def __repr__(self):
        rep = "Id WCC :"+str(self.id)
        rep += "\nLinks :"+str(self.links)
        return rep



###########################################
#           Streaming Version of WCC      #
###########################################


def weakly_connected_components(S):
    '''
    Compute the weakly connected components of S
    :param S: A Stream Graph
    :return: WCC of S
    '''
    comp_2_nb_links = {}  # Dictionary associating a comp to its number of current temporal links
    node_2_status = {}  # Dictionary associating a node to his current status : (current degree, id current comp)
    tmp_components = []  # List of local components
    final_components = []
    cnt_wcc_id = 0  # Id of stored wcc

    # get ordered links : edge stream
    # arrival (1,t0,t1,u,v) and departure (-1,t1,u,v)
    L = S.ordered_links()
    print("Links :",L)

    for l in L:
        c = l[0]
        if c == 1:
            u, v = l[3],l[4]
            if u not in node_2_status and v not in node_2_status:
                # print("CREATE")
                create_wcc(l, node_2_status, tmp_components, comp_2_nb_links)
                assert node_2_status[u][1] == node_2_status[v][1]

            elif u not in node_2_status and v in node_2_status:
                # print("UPDATE from :", v)
                update_wcc(v, u, l, node_2_status, tmp_components, comp_2_nb_links)
                assert node_2_status[u][1] == node_2_status[v][1]

            elif u in node_2_status and v not in node_2_status:
                # print("UPDATE from :", u)
                update_wcc(u, v, l, node_2_status, tmp_components, comp_2_nb_links)
                assert node_2_status[u][1] == node_2_status[v][1]

            elif node_2_status[u][1] != node_2_status[v][1]:
                # print("MERGE")
                merge_wccs(l, node_2_status, tmp_components, comp_2_nb_links)
                assert node_2_status[u][1] == node_2_status[v][1]
            else:
                node_2_status[u][0] += 1
                node_2_status[v][0] += 1
                n_current_comp = node_2_status[u][1]

                comp_2_nb_links[n_current_comp] += 1
                comp = tmp_components[n_current_comp]
                comp.links.append(l)
                comp.active_links.add((u,v))
                assert node_2_status[u][1] == node_2_status[v][1]
        elif c == -1:
            u, v = l[2],l[3]
            assert node_2_status[u][1] == node_2_status[v][1]
            n_current_comp = node_2_status[u][1]
            node_2_status[u][0] -= 1
            node_2_status[v][0] -= 1
            comp_2_nb_links[n_current_comp] -= 1

            comp = tmp_components[n_current_comp]
            comp.active_links.remove((u,v))
            if comp_2_nb_links[n_current_comp] == 0:
                # print("STORE COMP :", n_current_comp)
                cnt_wcc_id = end_comp(n_current_comp, tmp_components,final_components, cnt_wcc_id)

            if node_2_status[u][0] == 0:
                # print("Remove node :", u)
                del node_2_status[u]
            if node_2_status[v][0] == 0:
                # print("Remove node :", v)
                del node_2_status[v]
    # Store the rest of the wcc
    for n_comp in range(len(tmp_components)):
        if tmp_components[n_comp]:
            # print("STORE COMP")
            cnt_wcc_id = end_comp(n_comp, tmp_components,final_components, cnt_wcc_id)
    return final_components


def end_comp(n_current_comp, tmp_components, final_components, cnt_wcc_id):
    '''
    Store the component
    :param n_current_comp:
    :param tmp_components:
    :param storage_path:
    :return:
    '''
    comp_to_store = tmp_components[n_current_comp]
    comp_to_store.id = cnt_wcc_id
    final_components.append(comp_to_store)
    tmp_components[n_current_comp] = None
    cnt_wcc_id += 1
    return cnt_wcc_id


def create_wcc(l, node_2_status, tmp_components, comp_2_nb_links):
    '''
    Create a wcc component
    :param u:
    :param v:
    :param t0:
    :param t1:
    :param node_2_status:
    :param nb_comp:
    :param comp_2_nb_links:
    :return:
    '''
    u, v = l[3],l[4]
    n_comp = len(tmp_components)
    node_2_status[u] = [1, n_comp]
    node_2_status[v] = [1, n_comp]
    comp_2_nb_links[n_comp] = 1
    tmp_components.append(weakly_connected_component(links=[l], active_links=set([(u, v)])))


def update_wcc(node_to_update, node_to_add, l, node_2_status, tmp_components, comp_2_nb_links):
    '''
    Update a component with a new node
    :param node_to_update:
    :param node_to_add:
    :param t0:
    :param t1:
    :param node_2_status:
    :param tmp_components:
    :param comp_2_nb_links:
    :return:
    '''
    u,v = l[3],l[4]
    n_current_comp = node_2_status[node_to_update][1]
    comp_2_nb_links[n_current_comp] += 1
    node_2_status[node_to_update][0] += 1
    node_2_status[node_to_add] = [1, n_current_comp]

    current_comp = tmp_components[n_current_comp]
    current_comp.links.append(l)
    current_comp.active_links.add((u,v))


def merge_wccs(l, node_2_status, tmp_components, comp_2_nb_links):
    '''
    Merge two components
    :param l:
    :param node_2_status:
    :param tmp_components:
    :param comp_2_nb_links:
    :return:
    '''
    u, v = l[3],l[4]
    n_comp_1, n_comp_2 = node_2_status[u][1], node_2_status[v][1]
    comp_1 = tmp_components[n_comp_1]
    comp_2 = tmp_components[n_comp_2]

    comp_1.links += comp_2.links
    comp_1.active_links |= comp_2.active_links

    comp_1.links.append(l)
    comp_1.active_links.add((u,v))

    for h in comp_2.active_links:
        n, m = h
        if n in node_2_status and m in node_2_status:
            node_2_status[m][1] = n_comp_1
            node_2_status[n][1] = n_comp_1

    node_2_status[u][0] += 1
    node_2_status[v][0] += 1
    comp_2_nb_links[n_comp_1] += comp_2_nb_links[n_comp_2] + 1

    del comp_2_nb_links[n_comp_2]
    tmp_components[n_comp_2] = None

    node_2_status[v][1] = n_comp_1

if __name__ == '__main__':
    __directory__ = "/home/leo/Dev/Data_Stream/"
    __file__ = "sg_generated"
    S = sg.read_stream_graph(path_nodes=__directory__ + __file__ + "_nodes.sg",
                             path_links=__directory__ + __file__ + "_links.sg")
    S.plot()
    # wcc = S.weakly_connected_components()
    wcc = weakly_connected_components(S)
    for c in wcc:
        print("wcc :",c)

        # Sub = S.substream(c)
        # Sub.plot()
    # S.plot(clusters=wcc)
    plt.show()
