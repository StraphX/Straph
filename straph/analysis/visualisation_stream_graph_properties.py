import pathlib, math,gc
import msgpack
import matplotlib.pyplot as plt
import matplotlib.collections as mcol
from collections import defaultdict
import matplotlib
import numpy as np
import networkx as nx

from straph.utils import get_cmap


# ??? DATASHADER ??? useful when too big for matplotlib




def add_arrow(e, ax,direction='right', size=15, color=None,alpha=0.5):
    """
    Thanks to : https://stackoverflow.com/questions/34017866/arrow-on-a-line-plot-with-matplotlib
    add an arrow to a line.

    line:       ((a_x,_y),(b_x,b_y))
    ax:         matplotlib axes
    direction:  'left' or 'right'
    size:       size of the arrow in fontsize points
    color:      color of arrow (should be coherent with line color, or not )
    """
    print("line :",e)
    if direction == 'left':
        ax.annotate('',
                    xy=(e[0][0], e[0][1]),xycoords ='data',
                    xytext=(e[1][0], e[1][1]),textcoords = 'data',
                    arrowprops=dict(arrowstyle="<-",connectionstyle="arc3", color=color,alpha=alpha),
                    size=size,
        )
    else:
        ax.annotate('',
                    xy=(e[1][0], e[1][1]),xycoords ='data',
                    xytext=(e[0][0], e[0][1]),textcoords = 'data',
                    arrowprops=dict(arrowstyle="->",connectionstyle="arc3", color=color,alpha=alpha),
                    size=size,
        )

######################
#       PLOT WCC     #
######################
def parse_wcc_to_plot(wcc_storage_path):
    WCC = []
    ###############
    gc.disable()
    ###############
    with open(wcc_storage_path,'rb') as inpt:
        unpacker = msgpack.Unpacker(inpt, use_list=False)
        for i in unpacker:
            dict_node_end_time = defaultdict(int)
            dict_node_begin_time = defaultdict(int)
            dict_node_presence = defaultdict(list)
            present_nodes = set()
            for j in i[3]:
                if j[0] == 1:
                    t0, t1 = j[1], j[2]
                    u,v = j[3],j[4]
                    if u not in present_nodes and v not in present_nodes:
                        dict_node_begin_time[u] = t0
                        dict_node_begin_time[v] = t0
                        dict_node_end_time[u] = t1
                        dict_node_end_time[v] = t1
                        present_nodes.add(u)
                        present_nodes.add(v)
                    elif u not in present_nodes and v in present_nodes:
                        dict_node_end_time[v] = max(t1, dict_node_end_time[v])
                        dict_node_begin_time[u] = t0
                        dict_node_end_time[u] = t1
                        present_nodes.add(u)
                    elif u in present_nodes and v not in present_nodes:
                        dict_node_end_time[u] = max(t1, dict_node_end_time[u])
                        dict_node_begin_time[v] = t0
                        dict_node_end_time[v] = t1
                        present_nodes.add(v)
                    else:
                        dict_node_end_time[u] = max(t1, dict_node_end_time[u])
                        dict_node_end_time[v] = max(t1, dict_node_end_time[v])
                else:
                    t1 = j[1]
                    u,v = j[2],j[3]
                    if u in present_nodes and dict_node_end_time[u] == t1:
                        dict_node_presence[u] += [dict_node_begin_time[u], t1]
                        present_nodes.discard(u)
                        del dict_node_begin_time[u], dict_node_end_time[u]
                    if v in present_nodes and dict_node_end_time[v] == t1:
                        dict_node_presence[v] += [dict_node_begin_time[v], t1]
                        present_nodes.discard(v)
                        del dict_node_begin_time[v], dict_node_end_time[v]

            wcc = dict_node_presence
            WCC.append(wcc)
    ###############
    gc.enable()
    ###############
    return WCC


def plot_lines_wcc(wcc, label_to_id=False):
    seg = []
    for k,v in wcc.items():
        for t0,t1 in zip(v[::2],v[1::2]):
            if label_to_id:
                seg.append(((t0, label_to_id[k]), (t1, label_to_id[k])))
            else:
                seg.append(((t0, k), (t1, k)))
    return seg

def plot_wcc(wcc_storage_path, title=None,legend=None, saving_path=None, format='pdf',
             mask_biggest_comp=False,mask_big_comp=False, mask_small_comp=False,mask_medium_comp=False):
    WCC = parse_wcc_to_plot(wcc_storage_path)
    max_len = max([sum([len(v) for v in wcc.values()])/2 for wcc in WCC])
    nodes = set([i for wcc in WCC for i in wcc ])
    t_min = min([v[0] for wcc in WCC for v in wcc.values()])
    t_max = max([v[-1] for wcc in WCC for v in wcc.values()])
    print("T min :",t_min)
    print("T max :",t_max)
    print("Max len wcc :",max_len)
    colors = []
    segs = []
    color_small = "#47476b" # Light Gray/Blue
    color_medium = "#ffa31a" # Orange
    color_big = "#00b386" # Turquoise
    color_biggest = "#ff5050" # Rouge

    ####### TEMP #################
    # Transform nodes into int
    label_to_id = defaultdict(lambda :len(label_to_id))
    for wcc in WCC:
        for k,v in wcc.items():
            cnt = label_to_id[k]

    for wcc in WCC:
        len_wcc= sum([len(v) for v in wcc.values()])/2
        if len_wcc <= 2 and not mask_small_comp:
            seg = plot_lines_wcc(wcc,label_to_id)
            segs += seg
            colors+= [color_small]*len(seg)
        elif 2 < len_wcc <= 10 and not mask_medium_comp:
            seg = plot_lines_wcc(wcc,label_to_id)
            segs += seg
            colors+= [color_medium]*len(seg)
        elif len_wcc < max_len and not mask_big_comp:
            seg = plot_lines_wcc(wcc,label_to_id)
            segs += seg
            colors+= [color_big]*len(seg)
        elif not mask_biggest_comp:
            seg = plot_lines_wcc(wcc,label_to_id)
            segs += seg
            colors+= [color_biggest]*len(seg)
    print("len colors :",len(colors))
    print("len segs :",len(segs))
    print("Segs calculated")

    lnwdth = 300/len(nodes)
    if lnwdth > 0.8:
        lnwdth = 0.8
    line_coll = matplotlib.collections.LineCollection(np.array(segs),colors=colors,linewidths=lnwdth)#[1*10**(-20)])
    if mask_biggest_comp:
        line_coll.set_alpha(0.7)
    else:
        line_coll.set_alpha(0.6)
    print("Collections created")
    fig,ax = plt.subplots(1,1)
    ax.add_collection(line_coll)
    ax.set_xlim(t_min,t_max)
    ax.set_ylim(0,len(nodes)*1.05)
    print("Collections added")
    if legend :
        plt.ylabel("Nodes", fontname='Ubuntu', fontsize=12, color='#666699')
        plt.xlabel("t", fontname='Ubuntu', fontsize=12, color='#476b6b')
        list_legend =[]
        if not mask_small_comp:
            list_legend.append(matplotlib.patches.Patch(color=color_small, label="Small WCC"))
        if not mask_medium_comp:
            list_legend.append(matplotlib.patches.Patch(color=color_medium, label="Medium WCC"))
        if not mask_big_comp:
            list_legend.append(matplotlib.patches.Patch(color=color_big, label="Big WCC"))
        if not mask_biggest_comp:
            list_legend.append(matplotlib.patches.Patch(color=color_biggest, label="Biggest WCC"))
        plt.legend(handles=list_legend, loc='best')
    if title:
        ax.set_title(title, fontname='Ubuntu', fontsize=14)
    for spine in plt.gca().spines.values():
        spine.set_visible(False)
    plt.tick_params(top=False, bottom=True, right=False, left=True, labelleft=True, labelbottom=True)
    plt.tight_layout()
    if saving_path:
        fig.savefig(saving_path+"wcc_plot."+ format, format=format,dpi=2400)



############################
#       PLOT SCC
############################


def parse_scc_to_plot(scc_storage_path):
    SCC = []
    ###############
    gc.disable()
    ###############
    with open(scc_storage_path, 'rb') as inpt:
        unpacker = msgpack.Unpacker(inpt, use_list=False)
        for i in unpacker:
            dict_node_presence = defaultdict(list)
            times = i[2],i[3]
            for j in i[4]:
                if len(j)==2:
                    u, v = j
                else:
                    u,v = j[2],j[3]
                dict_node_presence[u] = times
                dict_node_presence[v] = times
            SCC.append(dict_node_presence)
    ###############
    gc.enable()
    ###############
    return SCC


def plot_lines_scc(scc,label_to_id=False):
    seg = []
    for k,v in scc.items():
        for t0,t1 in zip(v[::2],v[1::2]):
            if label_to_id:
                seg.append(((t0, label_to_id[k]), (t1, label_to_id[k])))
            else:
                seg.append(((t0, k), (t1, k)))
    return seg

def plot_scc(scc_storage_path, title=None, legend=None, saving_path=None, format='pdf',
             mask_biggest_comp=False, mask_big_comp=False, mask_small_comp=False, mask_medium_comp=False):
    SCC = parse_scc_to_plot(scc_storage_path)
    max_len = max([sum([len(v) for v in scc.values()]) / 2 for scc in SCC])
    nodes = set([i for scc in SCC for i in scc])
    t_min = min([v[0] for scc in SCC for v in scc.values()])
    t_max = max([v[-1] for scc in SCC for v in scc.values()])
    print("T min :",t_min)
    print("T max :",t_max)
    print("Max len SCC :",max_len)
    colors = []
    segs = []
    color_small = "#47476b" # Light Gray/Blue
    color_medium = "#ffa31a" # Orange
    color_big = "#00b386" # Turquoise
    color_biggest = "#ff5050" # Rouge

    ####### TEMP #################
    # Transform nodes into int
    label_to_id = defaultdict(lambda :len(label_to_id))
    for scc in SCC:
        for k,v in scc.items():
            cnt = label_to_id[k]

    for scc in SCC:
        len_wcc= sum([len(v) for v in scc.values()])/2
        if len_wcc <= 2 and not mask_small_comp:
            seg = plot_lines_scc(scc,label_to_id)
            segs += seg
            colors+= [color_small]*len(seg)
        elif 2 < len_wcc <= 10 and not mask_medium_comp:
            seg = plot_lines_scc(scc,label_to_id)
            segs += seg
            colors+= [color_medium]*len(seg)
        elif len_wcc < max_len and not mask_big_comp:
            seg = plot_lines_scc(scc,label_to_id)
            segs += seg
            colors+= [color_big]*len(seg)
        elif not mask_biggest_comp:
            seg = plot_lines_scc(scc,label_to_id)
            segs += seg
            colors+= [color_biggest]*len(seg)
    print("len colors :",len(colors))
    print("len segs :",len(segs))
    print("Segs calculated")


    lnwdth = 300/len(nodes)
    if lnwdth > 0.8:
        lnwdth = 0.8
    line_coll = matplotlib.collections.LineCollection(np.array(segs),colors=colors,
                                                      linewidths=lnwdth)#[1*10**(-20)])
    if mask_biggest_comp:
        line_coll.set_alpha(0.7)
    else:
        line_coll.set_alpha(0.6)
    print("Collections created")
    fig,ax = plt.subplots(1,1)
    ax.add_collection(line_coll)
    ax.set_xlim(t_min,t_max)
    ax.set_ylim(0,len(nodes)*1.05)
    print("Collections added")
    if legend :
        plt.ylabel("Nodes", fontname='Ubuntu', fontsize=12, color='#666699')
        plt.xlabel("t", fontname='Ubuntu', fontsize=12, color='#476b6b')
        list_legend =[]
        if not mask_small_comp:
            list_legend.append(matplotlib.patches.Patch(color=color_small, label="Small SCC"))
        if not mask_medium_comp:
            list_legend.append(matplotlib.patches.Patch(color=color_medium, label="Medium SCC"))
        if not mask_big_comp:
            list_legend.append(matplotlib.patches.Patch(color=color_big, label="Big SCC"))
        if not mask_biggest_comp:
            list_legend.append(matplotlib.patches.Patch(color=color_biggest, label="Biggest SCC"))
        plt.legend(handles=list_legend, loc='best')
    if title:
        ax.set_title(title, fontname='Ubuntu', fontsize=14)
    for spine in plt.gca().spines.values():
        spine.set_visible(False)
    plt.tick_params(top=False, bottom=True, right=False, left=True, labelleft=True, labelbottom=True)
    plt.tight_layout()
    if saving_path:
        fig.savefig(saving_path + "scc_plot." + format, format=format,dpi=2400)

####################################
#       PLOT Condensation DAG      #
####################################

def plot_condensation_dag_nx(SD_nodes, SD_links):
    g_adjacency_list = {n: [] for n in SD_nodes}
    if SD_links:
        for l in SD_links:
            n1 = l[0]
            n2 = l[1]
            g_adjacency_list[n1].append(n2)
        G = nx.from_dict_of_lists(g_adjacency_list)
        pos = nx.kamada_kawai_layout(G)
        nx.draw_networkx_nodes(G, pos, node_size=10,
                                       node_color="#339966", alpha=0.5)
        nx.draw_networkx_edges(G, pos, edge_color='#2d5986',
                                       alpha=0.5, width=1,arrows=True)

def load_times_scc(id_wcc, id_scc, scc_storage, dict_offset):
    with open(scc_storage, 'rb') as file_input:
        file_input.seek(dict_offset[(id_wcc, id_scc)])
        i = msgpack.Unpacker(file_input, use_list=False).unpack()
        times_scc = (i[2], i[3])
    return times_scc


def plot_condensation_dag_custom(id_wcc, SD_nodes, SD_links, scc_storage, offset_storage,arrow=False,label=True):
    dict_offset = msgpack.unpack(open(offset_storage, 'rb'), use_list=False)
    plt.figure()
    ax = plt.gca()
    g_adjacency_list = {n: [] for n in SD_nodes}
    pos = {n: [] for n in SD_nodes}
    y = 0
    min_x = math.inf
    max_x = -math.inf
    for id_scc in SD_nodes:     # Get pos of each node
        times_scc = load_times_scc(id_wcc, id_scc, scc_storage, dict_offset)
        t0, t1 = times_scc
        x = (t0 + t1) / 2
        min_x = min(min_x, t0)
        max_x = max(max_x, t1)
        pos[id_scc] = [x, id_scc]
    if SD_links:                # Get Adjacency list
        for l in SD_links:
            n1 = l[0]
            n2 = l[1]
            g_adjacency_list[n1].append(n2)
    # Plot nodes and links : plot nodes in increasing order
    # adjust 'y' depending on the number of neighbors
    visited = set()
    max_y = 0
    min_y = 0
    E = list()
    for n, _ in sorted(pos.items(), key=lambda x: x[1][0]):
        for u in g_adjacency_list[n]:
            if u not in visited:
                x, y = pos[u]
                max_y = max(max_y, y)
                visited.add(u)
                pos[u] = (x, y)
            E.append((pos[n], pos[u]))
    xy = np.asarray([pos[v] for v in pos])
    print("BEGIN PLOT, POSITIONS CALCULATED")
    ax.scatter(xy[:, 0], xy[:, 1],
               c="#339966",
               s=25,
               marker='o',
               alpha=0.7
               )
    if label:
        for n in pos:
            ax.annotate(n, pos[n], fontsize=10)
    print(" DRAWING NODES : DONE")
    edge_collections = mcol.LineCollection(E, colors=['#2d5986'], linewidths=2, alpha=0.5)
    ax.add_collection(edge_collections)
    if arrow == True:
        for e in E:
            add_arrow(e,ax, direction='right', color='#2d5986')
    print(" DRAWING EDGES : DONE")
    ax.set_ylim((min_y - 3, max_y + 3))
    ax.set_xlim((min_x, max_x))
    ax.set_xlabel("t", fontname='Ubuntu', fontsize=12, color='#476b6b')
    ax.tick_params(right=False, left=False, labelleft=False, labelbottom=True)

def plot_condensation_dag(dag_storage, scc_storage, offset_storage, figures_path,arrow=False):
    Unpacker = msgpack.Unpacker(open(dag_storage + "scc_dag_with_links.scf", 'rb'), use_list=False)
    format = "pdf"
    distrib_nodes_dag = []
    distrib_links_dag = []
    occurences_nodes = defaultdict(int)
    occurences_links = defaultdict(int)
    for i in Unpacker:
        id_wcc = i[0]
        nodes = i[1]
        links = i[2]
        n_nodes = len(nodes)
        distrib_nodes_dag.append(n_nodes)
        n_links = len(links)
        distrib_links_dag.append(n_links)
        for n in nodes:
            occurences_nodes[n] += 1
        for l in links:
            occurences_links[l] += 1
        # if 200 >= n_nodes >= 10 or 200 >= n_links >= 10:
        if 800 >= n_nodes >= 0 or 800 >= n_links >= 0:
            # Print big weakly component aka big graphs
            print("ID WCC : ", id_wcc, " Number of Nodes : ", n_nodes, " Number of Links : ", n_links)
            fig = plt.figure()
            plot_condensation_dag_custom(id_wcc, nodes, links, scc_storage, offset_storage,arrow=arrow)
            fig.savefig(figures_path+"condensation_dag_plot" + str(id_wcc) + "_custom." + format, format=format, dpi=2400)
            fig = plt.figure()
            plot_condensation_dag_nx(nodes, links)
            fig.savefig(figures_path+"condensation_dag_plot" + str(id_wcc) + "_nx." + format, format=format, dpi=2400)
        set_nodes = set(i[1])
        set_nodes_in_links = set()
        for l in i[2]:
            set_nodes_in_links.add(l[0])
            set_nodes_in_links.add(l[1])
        if set_nodes != set_nodes_in_links:
            # Check if there is some isolated nodes in the SCC DAG
            print("ID WCC :", i[0])
            print("len set nodes :", len(set_nodes))
            print("len set nodes in links :", len(set_nodes_in_links))
            print("number of isolated nodes :", len(set_nodes - set_nodes_in_links))
    # HISTOGRAMS :
    # ax = sns.distplot(distrib_nodes_dag, hist_kws={'alpha': 0.5, 'ec': '#510438', 'fc': '#a8508b', 'linewidth': '1'},
    #                   kde=False)
    # ax.set(yscale='log')
    # ax = sns.distplot(distrib_links_dag, hist_kws={'alpha': 0.5, 'ec': '#510438', 'fc': '#a8508b', 'linewidth': '1'},
    #                   kde=False)
    # ax.set(yscale='log')
    # ax = sns.distplot(list(occurences_nodes.values()),
    #                   hist_kws={'alpha': 0.5, 'ec': '#510438', 'fc': '#a8508b', 'linewidth': '1'}, kde=False)
    # ax.set(yscale='log')
    # ax = sns.distplot(list(occurences_links.values()),
    #                   hist_kws={'alpha': 0.5, 'ec': '#510438', 'fc': '#a8508b', 'linewidth': '1'}, kde=False)
    # ax.set(yscale='log')
    # plt.show()

############################
#       PLOT Hinge DAG
############################

def plot_hinge_dag_nx(SD_nodes, SD_links):
    g_adjacency_list = {n: [] for n in SD_nodes}
    if SD_links:
        for l in SD_links:
            n1 = l[0]
            n2 = l[1]
            g_adjacency_list[n1].append(n2)
        G = nx.from_dict_of_lists(g_adjacency_list)
        pos = nx.kamada_kawai_layout(G)
        nx.draw_networkx_nodes(G, pos, node_size=10,
                                       node_color="#339966", alpha=0.5)
        nx.draw_networkx_edges(G, pos, edge_color='#2d5986',
                                       alpha=0.5, width=1,arrows=True)

def plot_hinge_dag_custom(SD_nodes, SD_links):
    ax = plt.gca()
    min_x = min([n[0] for n in SD_nodes])
    max_x = max([n[0] for n in SD_nodes])
    min_y = min([n[1] for n in SD_nodes])
    max_y = max([n[1] for n in SD_nodes])
    xy = np.asarray(SD_nodes)
    print("BEGIN PLOT, POSITIONS CALCULATED")
    ax.scatter(xy[:, 0], xy[:, 1],
               c="#339966",
               s=10,
               marker='s',
               alpha=0.5
               )


    print(" DRAWING NODES : DONE")
    edge_collections = mcol.LineCollection(SD_links, colors=['#2d5986'], linewidths=1, alpha=0.5)
    ax.add_collection(edge_collections)
    print(" DRAWING EDGES : DONE")
    ax.set_ylim((min_y-3, max_y+3))
    ax.set_xlim((min_x-3, max_x+3))
    ax.tick_params(right=False, left=False, labelleft=False, labelbottom=True)
    ax.set_xlabel("t", fontname='Ubuntu', fontsize=12, color='#476b6b')

def plot_hinge_dag(dag_storage, scc_storage, offset_storage, figures_path):
    Unpacker = msgpack.Unpacker(open(dag_storage + "hinge_dag_with_links.scf", 'rb'), use_list=False)
    format = "pdf"
    distrib_nodes_dag = []
    distrib_links_dag = []
    occurences_nodes = defaultdict(int)
    occurences_links = defaultdict(int)
    for i in Unpacker:
        id_wcc = i[0]
        nodes = i[1]
        links = i[2]

        n_nodes = len(nodes)
        distrib_nodes_dag.append(n_nodes)
        n_links = len(links)
        distrib_links_dag.append(n_links)
        for n in nodes:
            occurences_nodes[n] += 1
        for l in links:
            occurences_links[l] += 1
        # if 200 >= n_nodes >= 0 or 200 >= n_links >= 0:
        if 800 >= n_nodes >= 50 or 800 >= n_links >= 50:
            # Print big weakly component aka big graphs
            print("ID WCC : ", id_wcc, " Number of Nodes : ", n_nodes, " Number of Links : ", n_links)
            fig = plt.figure()
            plot_hinge_dag_custom(nodes, links)
            fig.savefig(figures_path+"hinge_dag_plot" + str(id_wcc) + "_custom." + format, format=format, dpi=2400)
            fig = plt.figure()
            plot_hinge_dag_nx(nodes, links)
            fig.savefig(figures_path+"hinge_dag_plot" + str(id_wcc) + "_nx." + format, format=format, dpi=2400)
        set_nodes = set(i[1])
        set_nodes_in_links = set()
        for l in i[2]:
            set_nodes_in_links.add(l[0])
            set_nodes_in_links.add(l[1])
        if set_nodes != set_nodes_in_links:
            # Check if there is some isolated nodes in the SCC DAG
            print("ID WCC :", i[0])
            print("len set nodes :", len(set_nodes))
            print("len set nodes in links :", len(set_nodes_in_links))
            print("number of isolated nodes :", len(set_nodes - set_nodes_in_links))
    # HISTOGRAMS :
    # ax = sns.distplot(distrib_nodes_dag, hist_kws={'alpha': 0.5, 'ec': '#510438', 'fc': '#a8508b', 'linewidth': '1'},
    #                   kde=False)
    # ax.set(yscale='log')
    # ax = sns.distplot(distrib_links_dag, hist_kws={'alpha': 0.5, 'ec': '#510438', 'fc': '#a8508b', 'linewidth': '1'},
    #                   kde=False)
    # ax.set(yscale='log')
    # ax = sns.distplot(list(occurences_nodes.values()),
    #                   hist_kws={'alpha': 0.5, 'ec': '#510438', 'fc': '#a8508b', 'linewidth': '1'}, kde=False)
    # ax.set(yscale='log')
    # ax = sns.distplot(list(occurences_links.values()),
    #                   hist_kws={'alpha': 0.5, 'ec': '#510438', 'fc': '#a8508b', 'linewidth': '1'}, kde=False)
    # ax.set(yscale='log')
    # plt.show()


#############################
#       PLOT KCORES         #
#############################

def plot_lines_kcores(cores,label_to_id =False):
    seg = []
    for i in cores:
        t0, t1, n = i
        if label_to_id:
            seg.append(((t0, label_to_id[n]), (t1, label_to_id[n])))
        else:
            seg.append(((t0, n), (t1, n)))
    return seg


def plot_kcores(kcores_storage_path, title=None, saving_path=None,
                format='pdf'):
    Kcores = defaultdict(list)
    ###############
    gc.disable()
    ###############
    with open(kcores_storage_path,'rb') as ipt:
        unpacker = msgpack.Unpacker(ipt,use_list=False)
        for i in unpacker:
            k = i[0]
            for j in i[1]:
                Kcores[k].append(j)
    nodes = set([v[-1] for core in Kcores.values() for v in core])
    t_min = min([v[0] for core in Kcores.values() for v in core])
    t_max = max([v[0] for core in Kcores.values() for v in core])
    print("T min :",t_min)
    print("T max :",t_max)
    print("Max Kcores :",max(Kcores.keys()))
    c_map = get_cmap(max(Kcores.keys())+1,cmap='Blues')
    segs = []
    colors=[]

    ####### TEMP #################
    # Transform nodes into int
    label_to_id = defaultdict(lambda :len(label_to_id))
    for k,v in Kcores.items():
        for i in v:
            cnt = label_to_id[i[2]]

    for k,v in Kcores.items():
        seg = plot_lines_kcores(v,label_to_id)
        segs+=seg
        colors += [c_map(k+1)]*len(seg)
    print("len colors :",len(colors))
    print("len segs :",len(segs))
    print("Segs calculated")

    lnwdth = 300/len(nodes)
    if lnwdth > 0.8:
        lnwdth = 0.8
    line_coll = matplotlib.collections.LineCollection(np.array(segs),colors=colors,
                                                      linewidths=lnwdth)#[1*10**(-20)])
    print("Collections created and added")
    line_coll.set_alpha(0.6)
    print("Collections created")
    fig,ax = plt.subplots(1,1)
    ax.add_collection(line_coll)
    ax.set_xlim(t_min,t_max)
    ax.set_ylim(0, max(label_to_id.values()) * 1.05)
    #ax.set_ylim(0,max(nodes)*1.05)
    print("Collections added")
    if title:
        ax.set_title(title, fontname='Ubuntu', fontsize=14)
    for spine in plt.gca().spines.values():
        spine.set_visible(False)
    plt.tick_params(top=False, bottom=True, right=False, left=True, labelleft=True, labelbottom=True)
    plt.tight_layout()
    if saving_path:
        fig.savefig(saving_path + "kcores_plot." + format, format=format, dpi=2400)

###################################
#       PLOT KCLIQUES             #
###################################

def plot_lines_kcliques(cliques,label_to_id = False):
    seg = []
    for i in cliques:
        t0, t1, n = i
        if label_to_id:
            seg.append(((t0, label_to_id[n]), (t1, label_to_id[n])))
        else:
            seg.append(((t0, n), (t1, n)))
    return seg


def plot_kcliques(kcliques_storage_path, title=None, saving_path=None,
                  format='pdf'):
    Kcliques = defaultdict(set)
    ###############
    gc.disable()
    ###############
    with open(kcliques_storage_path, 'rb') as ipt:
        unpacker = msgpack.Unpacker(ipt,use_list=False)
        for i in unpacker:
            k = i[0]
            for j in i[1]:
                for n in j[2]:
                    e = (j[0],j[1],n)
                    Kcliques[k].add(e)
                    if e in Kcliques[k-1]:
                        Kcliques[k-1].remove(e)
    nodes = set([v[-1] for clique in Kcliques.values() for v in clique])
    t_min = min([v[0] for clique in Kcliques.values() for v in clique])
    t_max = max([v[0] for clique in Kcliques.values() for v in clique])
    print("T min :",t_min)
    print("T max :",t_max)
    print("Max Kcliques :",max(Kcliques.keys()))
    c_map = get_cmap(max(Kcliques.keys())+1,cmap='Purples')
    segs = []
    colors=[]

    ####### TEMP #################
    # Transform nodes into int
    label_to_id = defaultdict(lambda :len(label_to_id))
    for k,v in Kcliques.items():
        for i in v:
            cnt = label_to_id[i[2]]

    for k,v in Kcliques.items():
        seg = plot_lines_kcores(v,label_to_id)
        segs+=seg
        colors += [c_map(k+1)]*len(seg)
    print("len colors :",len(colors))
    print("len segs :",len(segs))
    print("Segs calculated")

    lnwdth = 300/len(nodes)
    if lnwdth > 0.8:
        lnwdth = 0.8
    line_coll = matplotlib.collections.LineCollection(np.array(segs),colors=colors,
                                                      linewidths=lnwdth)#[1*10**(-20)])
    print("Collections created and added")
    line_coll.set_alpha(0.6)
    print("Collections created")
    fig,ax = plt.subplots(1,1)
    ax.add_collection(line_coll)
    ax.set_xlim(t_min,t_max)
    ax.set_ylim(0, max(label_to_id.values()) * 1.05)
    #ax.set_ylim(0,max(nodes)*1.05)
    print("Collections added")
    if title:
        ax.set_title(title, fontname='Ubuntu', fontsize=14)
    for spine in plt.gca().spines.values():
        spine.set_visible(False)
    plt.tick_params(top=False, bottom=True, right=False, left=True, labelleft=True, labelbottom=True)
    plt.tight_layout()
    if saving_path:
        fig.savefig(saving_path + "kcliques_plot." + format, format=format, dpi=2400)

#########################################
#       GENERAL PLOTTEUR (#DenisBaupin) #
#########################################

def plot_properties():
    #TODO
    plot_kcliques()
    plot_kcores()
    plot_wcc()
    plot_scc()


if __name__ == '__main__':
    __directory__ = "/home/leo/Dev/Data_Stream/2018/04/"  # LIP6
    __filedate__ = "20180418"


    __directory__ = "/home/leo/Dev/Data_Stream/Socio_Patterns/Workplace/"
    __filedate__ = "workplace"

    __directory__ = "/home/leo/Dev/Data_Stream/Crawdad/Infocom/"
    __filedate__ = "infocom"
    __directory__ = "/home/leo/Dev/Data_Stream/Crawdad/Rollernet/"
    __filedate__ = "rollernet"

    __directory__ = "/home/leo/Dev/Data_Stream/Bitcoin/"
    __filedate__ = "bitcoin_otc"

    __directory__ = "/home/leo/Dev/Data_Stream/CollegeMsg/"
    __filedate__ = "CollegeMsg"
    __directory__ = "/home/leo/Dev/Data_Stream/Socio_Patterns/High_School_2013/"
    __filedate__ = "High_School_2013"

    # __directory__ = "/home/leo/Dev/Data_Stream/Enron/"
    # __file__ = "enron"
    #
    # __directory__ = "/home/leo/Dev/Data_Stream/parse_net/11September/"
    # __file__ = "Days"

    __directory__ = "/home/leo/Dev/Data_Stream/example_paper/"
    __filedate__ ="example_paper"

    __directory__ = "/home/leo/Dev/Data_Stream/ITS/"
    __filedate__ = "Suricata"


    wcc_storage_path = __directory__ + __filedate__ + "_wcc_storage/wcc.scf"
    scc_storage_path = __directory__ + __filedate__ + "_scc_storage/scc.scf"
    condensation_dag_storage_path = __directory__ + __filedate__ + "_scc_dag_storage/"
    hinge_dag_storage_path = __directory__ + __filedate__ + "_hinge_dag_storage/"

    offset_storage = condensation_dag_storage_path + "dict_offset_scc.mpck"
    kcores_storage_path = __directory__ + __filedate__ + "_kcores_storage/postprocess_kcores.scf"
    kcliques_storage_path = __directory__ + __filedate__ + "_kcliques_storage/postprocess_kcliques.scf"


    figures_path = __directory__+ "figures/"
    pathlib.Path(figures_path).mkdir(parents=True, exist_ok=True)

    dag_plot_path = figures_path+"dag_plots/"
    pathlib.Path(dag_plot_path).mkdir(parents=True, exist_ok=True)

    format = 'pdf'


    # print()
    # plot_condensation_dag(condensation_dag_storage_path, scc_storage_path, offset_storage, dag_plot_path,arrow=True)
    # plot_hinge_dag(hinge_dag_storage_path, scc_storage_path, offset_storage, dag_plot_path)
    # exit()

    # plot_wcc(wcc_storage_path, title="Weakly Connected Components", legend=False, saving_path=figures_path, format=format,
    #              mask_biggest_comp=False, mask_big_comp=False, mask_small_comp=False, mask_medium_comp=False)
    # print()
    # plot_scc(scc_storage_path, title="Strongly Connected Components", legend=False, saving_path=figures_path, format=format,
    #              mask_biggest_comp=False, mask_big_comp=False, mask_small_comp=False, mask_medium_comp=False)
    # print()
    plot_kcores(kcores_storage_path,title="K-cores",saving_path=figures_path,format=format)
    # print()
    plot_kcliques(kcliques_storage_path, title="K-cliques", saving_path=figures_path, format=format)
