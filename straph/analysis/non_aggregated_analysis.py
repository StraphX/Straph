import numpy as np

import matplotlib.pyplot as plt
import pandas as pd

import csv
import copy

import h5py


__directory__ = "/home/leo/Dev/Data_Stream/2018/04/"
__filedate__ = "20180418"


# Pour chaque t écrire le nombre de noeud présent(à priori dans nodes, subordonné à node_presence)
# De même pour les liens, écrire le nombre de liens en cours et compter le lien courant
def extract_exact_nodes(file_input,
                        file_output,
                        time_begin,
                        time_break,
                        time_bins
                        ):
    compteur = 0

    nodes = []
    node_presence = []
    links = []
    link_presence = []
    dict_node_presence = {}
    dict_link_presence = {}

    current_nodes = set()
    current_links = set()

    with open(file_input, 'r') as input, open(file_output, 'w') as output:
        writer = csv.writer(output, delimiter=';')
        writer.writerow(["time", "nb_nodes", "nb_links"])
        reader = csv.reader(input, delimiter=';')
        for i in reader:
            compteur += 1
            # Ignore header
            if compteur == 1:
                continue
            line = i
            id_time = float(line[0])
            if id_time < time_begin:
                continue
            if line[1] and line[2]:
                try:
                    id_n1 = dict_node_presence[line[1]]
                    current_nodes.add(id_n1)
                    if node_presence[id_n1][-1] >= id_time:
                        node_presence[id_n1][-1] = id_time + time_bins
                    else:
                        node_presence[id_n1] += [id_time, id_time + time_bins]
                except:
                    nodes.append(line[1])
                    node_presence.append([id_time, id_time + time_bins])
                    dict_node_presence[line[1]] = len(nodes) - 1
                    current_nodes.add(len(nodes) - 1)
                try:
                    id_n2 = dict_node_presence[line[2]]
                    current_nodes.add(id_n2)
                    if node_presence[id_n2][-1] >= id_time:
                        node_presence[id_n2][-1] = id_time + time_bins
                    else:
                        node_presence[id_n2] += [id_time, id_time + time_bins]
                except:
                    nodes.append(line[2])
                    node_presence.append([id_time, id_time + time_bins])
                    dict_node_presence[line[2]] = len(nodes) - 1
                    current_nodes.add(len(nodes) - 1)

                try:
                    id_link = dict_link_presence[(line[1], line[2])]
                    current_links.add(id_link)
                    if link_presence[id_link][-1] >= id_time:
                        link_presence[id_link][-1] = id_time + time_bins
                    else:
                        link_presence[id_link] += [id_time, id_time + time_bins]
                except:
                    links.append((line[1], line[2]))
                    link_presence.append([id_time, id_time + time_bins])
                    dict_link_presence[(line[1], line[2])] = len(links) - 1
                    current_links.add(len(links) - 1)

                #  UPDATE CURRENT NODES AND LINKS:
                for n in copy.copy(current_nodes):
                    if node_presence[n][-1] < id_time:
                        current_nodes.remove(n)
                for l in copy.copy(current_links):
                    if link_presence[l][-1] < id_time:
                        current_links.remove(l)

                writer.writerow([id_time, len(current_nodes), len(current_links)])
                if id_time + time_bins > time_break:
                    break
                    # LAST PASS :
                    # DO SOMETHING
                    # Write current_nodes, current_links
                    # END LAST PASS
    return


def is_right_on_time(l_nb_nodes):   #,l_nb_links):
    percent_previous = 0.4
    percent_mean = 0.4
    # 1st rule : nb_nodes at t > 1.4*nb_nodes at t-1 (respectively < and 0.6)
    #            nb_nodes at t > 1.25*m_nodes on t-1,...,t-50 (respectively < and 0.75)
    m_nodes = np.mean(l_nb_nodes)
    # loc_max_nodes = np.max(l_nb_nodes)
    # loc_min_nodes = np.min(l_nb_nodes)
    # l = []
    if l_nb_nodes[-1] > (1+percent_previous)*l_nb_nodes[-2] or l_nb_nodes[-1] > (1+percent_mean)*m_nodes :
        return True
    if l_nb_nodes[-1] < (1-percent_previous)*l_nb_nodes[-2] or l_nb_nodes[-1] < (1-percent_mean)*m_nodes :
        return True
    return False


def find_times_of_interest(file_input,
                        file_output,
                        time_begin,
                        time_break,
                        time_bins,
                        len_slide = 200):
    compteur = 0

    nodes = []
    node_presence = []
    links = []
    link_presence = []
    dict_node_presence = {}
    dict_link_presence = {}

    current_nodes = set()
    current_links = set()

    l_nb_nodes = []
    l_nb_nodes_interest = []

    l_nb_links = []

    time_of_interest = []
    sliding_nb_nodes = []
    sliding_nb_links = []

    l_time = []

    with open(file_input, 'r') as input, h5py.File(file_output, 'w') as output:
        reader = csv.reader(input, delimiter=';')
        for i in reader:
            compteur += 1
            # Ignore header
            if compteur == 1:
                continue
            line = i
            id_time = float(line[0])
            if id_time < time_begin:
                continue
            if line[1] and line[2]:
                try:
                    id_n1 = dict_node_presence[line[1]]
                    current_nodes.add(id_n1)
                    if node_presence[id_n1][-1] >= id_time:
                        node_presence[id_n1][-1] = id_time + time_bins
                    else:
                        node_presence[id_n1] += [id_time, id_time + time_bins]
                except:
                    nodes.append(line[1])
                    node_presence.append([id_time, id_time + time_bins])
                    dict_node_presence[line[1]] = len(nodes) - 1
                    current_nodes.add(len(nodes) - 1)
                try:
                    id_n2 = dict_node_presence[line[2]]
                    current_nodes.add(id_n2)
                    if node_presence[id_n2][-1] >= id_time:
                        node_presence[id_n2][-1] = id_time + time_bins
                    else:
                        node_presence[id_n2] += [id_time, id_time + time_bins]
                except:
                    nodes.append(line[2])
                    node_presence.append([id_time, id_time + time_bins])
                    dict_node_presence[line[2]] = len(nodes) - 1
                    current_nodes.add(len(nodes) - 1)

                try:
                    id_link = dict_link_presence[(line[1], line[2])]
                    current_links.add(id_link)
                    if link_presence[id_link][-1] >= id_time:
                        link_presence[id_link][-1] = id_time + time_bins
                    else:
                        link_presence[id_link] += [id_time, id_time + time_bins]
                except:
                    links.append((line[1], line[2]))
                    link_presence.append([id_time, id_time + time_bins])
                    dict_link_presence[(line[1], line[2])] = len(links) - 1
                    current_links.add(len(links) - 1)

                #  UPDATE CURRENT NODES AND LINKS:
                for n in copy.copy(current_nodes):
                    if node_presence[n][-1] < id_time:
                        current_nodes.remove(n)
                for l in copy.copy(current_links):
                    if link_presence[l][-1] < id_time:
                        current_links.remove(l)

                # SLIDING SHIIT
                sliding_nb_nodes.append(len(current_nodes))
                if len(sliding_nb_nodes) > len_slide:
                    sliding_nb_nodes.remove(sliding_nb_nodes[0])
                sliding_nb_links.append(len(current_links))
                if len(sliding_nb_links) > len_slide:
                    sliding_nb_links.remove(sliding_nb_links[0])
                    # Test if the time IS RIGHT m*f*
                    if is_right_on_time(sliding_nb_nodes):
                        time_of_interest.append(id_time)
                        l_nb_nodes_interest.append(len(current_nodes))
                    elif len(time_of_interest) > 2:
                        if time_of_interest[-1] - time_of_interest[-2] > 0.1:
                            time_of_interest.append(id_time)
                            l_nb_nodes_interest.append(len(current_nodes))
                l_nb_nodes.append(len(current_nodes))
                l_nb_links.append(len(current_links))
                l_time.append(id_time)
                if id_time + time_bins > time_break:
                    break
                    # LAST PASS :
                    # DO SOMETHING
                    # Write current_nodes, current_links
                    # END LAST PASS
        time_of_interest.sort()

        output["time"] = np.array(l_time)
        output["ToI"] = np.array(time_of_interest)
        output["nb_nodes"] = np.array(l_nb_nodes)
        output["nb_links"] = np.array(l_nb_links)
        output["nb_nodes_interest"] = np.array(l_nb_nodes_interest)

    return

def plot_distributions(file,protocol="",path_save=None):
    df = pd.read_csv(file, delimiter=";")#,nrows= 50000000)
    df.set_index("time", inplace=True)
    df.drop("time", inplace=True)
    # df =df.iloc[50000000:,:]
    path_save_1 = __directory__+"/Results_exact_analysis/hist_"
    path_save_2 = __directory__+"/Results_exact_analysis/plot_time_"
    path_save_3 = __directory__+"/Results_exact_analysis/cdf_"
    for c in df.columns:
        if c == "nb_nodes":
            col = "darkviolet"
        else:
            col = "slateblue"
        ax = df[c].plot(kind="hist",bins=30,color=col)
        ax.set_xlabel(c)
        ax.set_ylabel("# Occurrences")
        plt.tight_layout()
        if path_save:
            plt.savefig(path_save_1+c+"_"+protocol+".pdf", format='pdf', dpi=200)
            plt.clf()

        ax = df[c].plot(color=col)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel(c)
        plt.tight_layout()
        if path_save:
            plt.savefig(path_save_2+c+"_"+protocol+".pdf", format='pdf', dpi=200)
            plt.clf()

        ax = df[c].plot(kind="hist",cumulative=True,normed=True,bins=200, color=col)
        ax.set_ylabel("P("+c+" < x)")
        ax.set_xlabel("x")
        plt.tight_layout()
        if path_save:
            plt.savefig(path_save_3 + c + "_" + protocol + ".pdf", format='pdf', dpi=200)
            plt.clf()
    del df


def plot_time_of_interest(storage,protocol,path_save = None):
    time_of_interest = storage["ToI"]
    id_time = storage["time"]
    nb_nodes = storage["nb_nodes"]
    nb_nodes_of_interest = storage["nb_nodes_interest"]
    print("nb nodes :",len(nb_nodes),"nodes of interest :",len(nb_nodes_of_interest))
    plt.plot(id_time, nb_nodes, color="slateblue",alpha=0.4,linewidth = 1)
    plt.plot(time_of_interest,nb_nodes_of_interest,color="firebrick",alpha=0.6,linewidth=1)
    plt.tight_layout()
    if path_save:
        plt.savefig(__directory__ + "/Results_exact_analysis/ToI_" + protocol + ".pdf", format='pdf', dpi=400)
    else:
        plt.show()
    plt.clf()


if __name__ == '__main__':
    protocols = ["UDP", "ICMP", "TCP", "GRE"]
    # protocols = ["UDP","GRE"]
    ###########################################
    #   EXTRACT STREAM GRAPH                  #
    ###########################################
    for protocol in protocols:
        print("Current protocol :",protocol)
        path_csv = __directory__ + __filedate__ + "_" + protocol + ".csv"
        time_begin = 0
        time_bins = 0.01
        time_break = 10000
        time_step = 30
        path_exact_stats = __directory__ + __filedate__ + "_exact_count_" + protocol + ".csv"
        # extract_exact_nodes(path_csv,path_exact_stats,time_begin,time_break,time_bins,protocol)
        # plot_distributions(path_exact_stats,protocol,path_save=True)
        path_time_of_interest = __directory__ + __filedate__ + "_" + protocol + ".hdf5"
        find_times_of_interest(path_csv,
                               path_time_of_interest,
                               time_begin,
                               time_break,
                               time_bins)
        plot_time_of_interest(h5py.File(path_time_of_interest,'r'),protocol,path_save=True)
