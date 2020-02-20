import re, pathlib, subprocess,time,gc
import matplotlib.pyplot as plt
import msgpack
from collections import defaultdict,Counter

import pandas as pd
import seaborn as sns
import numpy as np

import straph.dag.condensation_dag as sd
import straph.stream as ss

from straph.utils import dataframe_to_latex_landscape



###########################################
#       WCC Properties                    #
###########################################


def preprocess_wcc_into_dataframe(path_wcc, output_path):
    '''
    Collect properties of Weakly Connected Component into a pandas DataFrame
    :param path_wcc: Storage path of weakly connected components
    :param output_path: Path to store the DataFrame
    :return: A DataFrame containing WCC properties
    '''
    names = []
    # Properties to collect must be specified and computed manually
    dict_data = {'Nb_of_nodes': [],
                 'Nb_of_links': [],
                 'Nb_of_segmented_links':[],
                 'Duration': [],
                 'Surface': []}
    ###############
    gc.disable()
    ###############
    with open(path_wcc,'rb') as input:
        unpacker = msgpack.Unpacker(input,use_list=False)
        for i in unpacker:
            id_wcc = i[0]
            names.append(id_wcc)
            nodes = set()
            link_presence = defaultdict(list)
            times = i[1],i[2]
            for j in i[3]:
                if j[0] == 1:
                    t0, t1 = j[1], j[2]
                    u, v = j[3],j[4]
                    link_presence[(u, v)] += [t0, t1]
                    nodes.add(u)
                    nodes.add(v)
            dict_data['Nb_of_nodes'].append(len(nodes))
            dict_data['Nb_of_links'].append(len(link_presence))
            dict_data['Nb_of_segmented_links'].append(sum([len(v)/2 for v in link_presence.values()]))
            dict_data['Duration'].append(times[1] - times[0])
            dict_data['Surface'].append(
                sum([t1 - t0 for v in link_presence.values() for t0, t1 in zip(v[::2], v[1::2])]))
    Index = pd.Index(names)
    D = pd.DataFrame(dict_data, index=Index)
    if output_path:
        pathlib.Path(output_path).mkdir(parents=True, exist_ok=True)
        D.to_msgpack(output_path + "df_wcc.mspk")
    ###############
    gc.enable()
    ###############
    return D

def stats_wcc(path_figure, path_stats, format='pdf',
              storage_path=None,with_anomalies = False,
              ):
    pathlib.Path(path_figure).mkdir(parents=True, exist_ok=True)
    if storage_path:
        D = pd.read_msgpack(storage_path + "df_wcc.mspk")
        if with_anomalies:
            with open(storage_path + "anomalies_in_wcc.mspk", 'rb') as input:
                wcc_with_anomalies = set(msgpack.load(input,use_list=False))
            print("WCC with anomalies :", wcc_with_anomalies)
    dataframe_to_latex_landscape(D, path_stats,title= "Weakly Connected Components")
    D['id'] = [i for i in D.index]
    if with_anomalies:
        D_without_anomalies = D.loc[list(set(D.index) - wcc_with_anomalies)]
        D_with_anomalies = D.loc[list(wcc_with_anomalies)]
        print("D with anomalies :",D_with_anomalies.shape)
        print("D without anomalies :",D_without_anomalies.shape)

    # columns_with_xlogscale = set(['Nb_of_links','Nb_of_nodes','Surface','Nb_of_segmented_links'])
    columns_with_xlogscale = set()
    for c in D.columns:
        if c != 'id':
            print("Plotting ", c)
            fig, ax = plt.subplots(1, 1)
            if c in columns_with_xlogscale:
                custom_bins = np.logspace(np.log(np.min(D[c])),np.log(np.max(D[c])),50)
                # custom_bins = np.logspace(np.min(D[c]),np.max(D[c]),50)
                ax = sns.distplot(D[c], bins=custom_bins, ax=ax,
                                  hist_kws={'alpha': 0.5, 'ec': '#510438', 'fc': '#a8508b', 'linewidth': '1'}, kde=False)
                ax.set(xscale='log')
            else:
                ax = sns.distplot(D[c], bins=50, ax=ax,
                                  hist_kws={'alpha': 0.5, 'ec': '#510438', 'fc': '#a8508b', 'linewidth': '1'}, kde=False)
            ax.set_xlabel(c)
            ax.set(yscale='log')
            ax.set_title("Histogram of the " + c)
            fig.savefig(path_figure + "/histogram_" + c + "_wcc.pdf", format=format)
            plt.close(fig)

            if with_anomalies:
                fig, ax = plt.subplots(1, 1)
                ax = D_without_anomalies.plot(x='id', y=c, kind='scatter', ax=ax,s=10, alpha=0.3, color='#a8508b', linewidth='0.5')
                ax = D_with_anomalies.plot(x='id', y=c, kind='scatter', ax=ax,s=20, alpha=0.7, color='#cc2900', linewidth='0.5')
                ax.set_xlabel(c)
                ax.set(yscale="log")
                ax.set_title("Scatter plot of the " + c)
                fig.savefig(path_figure + "/scatter_plot_" + c + "_wcc.pdf", format=format)
                plt.close(fig)
    D.drop(columns=['id'], inplace=True)



###########################################
#       SCC Properties                    #
###########################################


def preprocess_scc_into_dataframe(path_scc, output_path):
    '''
    Collect properties of Strongly Connected Component into a pandas DataFrame
    :param path_wcc: Storage path of strongly connected components
    :param output_path: Path to store the DataFrame
    :return: A DataFrame containing SCC properties
    '''
    names = []
    # Properties to collect must be specified and computed manually
    dict_data = {'Nb_of_nodes': [],
                 'Nb_of_links': [],
                 'Nb_of_segmented_links':[],
                 'Duration': [],
                 'Surface': []}
    t_begin = time.time()
    ###############
    gc.disable()
    ###############
    with open(path_scc, 'rb') as input:
        unpacker = msgpack.Unpacker(input,use_list=False)
        for i in unpacker:
            id_wcc,id_scc = i[0],i[1]
            times = i[2],i[3]
            if id_wcc % 100000 == 0 or (id_scc % 100000 == 0 and id_scc !=0):
                print("id wcc :",id_wcc,"id scc :",id_scc)
                print("times :",times,"size scc :",len(i[4]), " t :", time.time() - t_begin)
            names.append((id_wcc,id_scc))
            nodes = set()
            link_presence = defaultdict(list)
            for j in i[4]:
                if len(j) == 2:
                    u, v = j
                    t0, t1 = times
                    link_presence[(u, v)] += [t0, t1]
                    nodes.add(u)
                    nodes.add(v)
                else:
                    t0, t1 = j[0], j[1]
                    u, v = j[2],j[3]
                    link_presence[(u, v)] += [t0, t1]
                    nodes.add(u)
                    nodes.add(v)
            dict_data['Nb_of_nodes'].append(len(nodes))
            dict_data['Nb_of_links'].append(len(link_presence))
            dict_data['Nb_of_segmented_links'].append(sum([len(v)/2 for v in link_presence.values()]))
            dict_data['Duration'].append(times[1] - times[0])
            dict_data['Surface'].append(
                sum([t1 - t0 for v in link_presence.values() for t0, t1 in zip(v[::2], v[1::2])]))
    Index = pd.Index(names)
    D = pd.DataFrame(dict_data, index=Index)
    if output_path:
        pathlib.Path(output_path).mkdir(parents=True, exist_ok=True)
        D.to_msgpack(output_path + "df_scc.mspk")
    ###############
    gc.enable()
    ###############
    return D


def stats_scc(path_figure, path_stats, format='pdf',
              storage_path=None, with_anomalies = False):
    pathlib.Path(path_figure).mkdir(parents=True, exist_ok=True)
    if storage_path:
        D = pd.read_msgpack(storage_path + "df_scc.mspk")
        if with_anomalies:
            with open(storage_path + "anomalies_in_scc.mspk", 'rb') as input:
                scc_with_anomalies = set(msgpack.load(input,use_list=False))
            # print("SCC with anomalies :", scc_with_anomalies)
    describe_to_latex_landscape(D, path_stats,title= "Strongly Connected Components")
    D['id'] = [i for i in D.index]
    if with_anomalies:
        D_without_anomalies = D.loc[list(set(D.index) - scc_with_anomalies)]
        D_with_anomalies = D.loc[list(scc_with_anomalies)]
        print("D with anomalies :",D_with_anomalies.shape)
        print("D without anomalies :",D_without_anomalies.shape)

    columns_with_xlogscale = set()

    for c in D.columns:
        if c != 'id':
            print("Plotting ", c)
            fig, ax = plt.subplots(1, 1)
            if c in columns_with_xlogscale:
                custom_bins = np.logspace(np.log(np.min(D[c])),np.log(np.max(D[c])),50)
                # custom_bins = np.logspace(np.min(D[c]),np.max(D[c]),50)
                ax = sns.distplot(D[c], bins=custom_bins, ax=ax,
                                  hist_kws={'alpha': 0.5, 'ec': '#510438', 'fc': '#a8508b', 'linewidth': '1'}, kde=False)
                ax.set(xscale='log')
            else:
                ax = sns.distplot(D[c], bins=50, ax=ax,
                                  hist_kws={'alpha': 0.5, 'ec': '#510438', 'fc': '#a8508b', 'linewidth': '1'}, kde=False)
            ax.set_xlabel(c)
            ax.set(yscale="log")
            ax.set_title("Histogram of the " + c)
            fig.savefig(path_figure + "/histogram_" + c + "_scc.pdf", format=format)
            plt.close(fig)

            # fig, ax = plt.subplots(1, 1)
            # # ax = D_without_anomalies.plot(x='id', y=c, kind='scatter', ax=ax,s=10, alpha=0.3, color='#a8508b', linewidth='0.5')
            # ax.scatter(D_without_anomalies['id'],D_without_anomalies[c],s=10, alpha=0.3, color='#a8508b', linewidth='0.5')
            # # ax = D_with_anomalies.plot(x='id', y=c, kind='scatter', ax=ax,s=20, alpha=0.7, color='#cc2900', linewidth='0.5')
            # ax.scatter(D_with_anomalies['id'],D_with_anomalies[c],s=20, alpha=0.7, color='#cc2900', linewidth='0.5')
            # ax.set_xlabel(c)
            # ax.set(yscale="log")
            # ax.set_title("Scatter plot of the " + c)
            # fig.savefig(path_figure + "/scatter_plot_" + c + "_scc.pdf", format=format)
            # plt.close(fig)

    D.drop(columns=['id'], inplace=True)


###########################################
#       K Cores Properties                #
###########################################

def get_core_number(storage_path_kcores):
    '''
    Preprocess kcores in order to return core numbers per nodes
    :param storage_path_kcores:
    :return:
    '''
    nodes_cores_stats = {}
    ###############
    gc.disable()
    ###############
    with open(storage_path_kcores, 'rb') as ipt:
        unpacker = msgpack.Unpacker(ipt, use_list=False)
        for i in unpacker:
            print("i :",i)
            k = i[0]
            for j in i[1]:
                t0, t1, u = j
                if u not in nodes_cores_stats:
                    nodes_cores_stats[u] = {k: [t0, t1]}
                elif k not in nodes_cores_stats[u]:
                    nodes_cores_stats[u][k] = [t0, t1]
                else:
                    nodes_cores_stats[u][k] += [t0, t1]
    ###############
    gc.disable()
    ###############
    return nodes_cores_stats



def preprocess_kcores_into_dataframe(storage_path_kcores, output_path):
    '''
    Create a DataFrame with node's core numbers
    columns : (k, frequency in kcores, surface per cores, intervals inter core ?, interval in core ?)
    :param storage_path_kcores:
    :param output_path:
    :return:
    '''
    nodes_cores_stats = get_core_number(storage_path_kcores)
    names = []
    dict_data = {'2-core': [],
                 '3-core': [],
                 '4-core': [],
                 'frequency 2-core': [],
                 'frequency 3-core': [],
                 'frequency 4-core': [],
                 'frequency': [],
                 'surface 2-core': [],
                 'surface 3-core': [],
                 'surface 4-core': [],
                 'surface': [],
                 # 'time between core':[] ?
                 }
    for u in nodes_cores_stats:
        names.append(u)
        if 2 in nodes_cores_stats[u]:
            dict_data['2-core'].append(1)
            dict_data['frequency 2-core'].append(len(nodes_cores_stats[u][2]) / 2)
            dict_data['surface 2-core'].append(
                sum([t1 - t0 for t0, t1 in zip(nodes_cores_stats[u][2][::2], nodes_cores_stats[u][2][1::2])]))
        else:
            dict_data['2-core'].append(0)
            dict_data['frequency 2-core'].append(0)
            dict_data['surface 2-core'].append(0)

        if 3 in nodes_cores_stats[u]:
            dict_data['3-core'].append(1)
            dict_data['frequency 3-core'].append(len(nodes_cores_stats[u][3]) / 2)
            dict_data['surface 3-core'].append(
                sum([t1 - t0 for t0, t1 in zip(nodes_cores_stats[u][3][::2], nodes_cores_stats[u][3][1::2])]))

        else:
            dict_data['3-core'].append(0)
            dict_data['frequency 3-core'].append(0)
            dict_data['surface 3-core'].append(0)

        if 4 in nodes_cores_stats[u]:
            dict_data['4-core'].append(1)
            dict_data['frequency 4-core'].append(len(nodes_cores_stats[u][4]) / 2)
            dict_data['surface 4-core'].append(
                sum([t1 - t0 for t0, t1 in zip(nodes_cores_stats[u][4][::2], nodes_cores_stats[u][4][1::2])]))
        else:
            dict_data['4-core'].append(0)
            dict_data['frequency 4-core'].append(0)
            dict_data['surface 4-core'].append(0)

        dict_data['frequency'].append(
            dict_data['frequency 2-core'][-1] + dict_data['frequency 3-core'][-1] + dict_data['frequency 4-core'][-1])
        dict_data['surface'].append(
            dict_data['surface 2-core'][-1] + dict_data['surface 3-core'][-1] + dict_data['surface 4-core'][-1])
    Index = pd.Index(names)
    D = pd.DataFrame(dict_data, index=Index)
    if output_path:
        pathlib.Path(output_path).mkdir(parents=True, exist_ok=True)
        D.to_msgpack(output_path + "df_kcores.mspk")
    return D


def stats_kcores(path_figure, path_stats, format='pdf',
              storage_path=None, with_anomalies = False):
    '''

    :param path_figure:
    :param path_stats:
    :param format:
    :param storage_path:
    :param with_anomalies:
    :return:
    '''
    pathlib.Path(path_figure).mkdir(parents=True, exist_ok=True)
    if storage_path:
        D = pd.read_msgpack(storage_path + "df_kcores.mspk")
        if with_anomalies:
            with open(storage_path + "anomalies_in_kcores.mspk", 'rb') as input:
                kcores_with_anomalies = set(msgpack.load(input,use_list=False))
            # print("SCC with anomalies :", scc_with_anomalies)
    describe_to_latex_landscape(D, path_stats,title= "K Cores")
    D['id'] = [i for i in D.index]
    if with_anomalies:
        D_without_anomalies = D.loc[list(set(D.index) - kcores_with_anomalies)]
        D_with_anomalies = D.loc[list(kcores_with_anomalies)]
        print("D with anomalies :",D_with_anomalies.shape)
        print("D without anomalies :",D_without_anomalies.shape)

    columns_with_xlogscale = set()

    for c in D.columns:
        if c != 'id':
            print("Plotting ", c)
            fig, ax = plt.subplots(1, 1)
            if c in columns_with_xlogscale:
                custom_bins = np.logspace(np.log(np.min(D[c])),np.log(np.max(D[c])),50)
                # custom_bins = np.logspace(np.min(D[c]),np.max(D[c]),50)
                ax = sns.distplot(D[c], bins=custom_bins, ax=ax,
                                  hist_kws={'alpha': 0.5, 'ec': '#510438', 'fc': '#a8508b', 'linewidth': '1'}, kde=False)
                ax.set(xscale='log')
            else:
                ax = sns.distplot(D[c], bins=50, ax=ax,
                                  hist_kws={'alpha': 0.5, 'ec': '#510438', 'fc': '#a8508b', 'linewidth': '1'}, kde=False)
            ax.set_xlabel(c)
            ax.set(yscale="log")
            ax.set_title("Histogram of the " + c)
            fig.savefig(path_figure + "/histogram_" + c + "_kcores.pdf", format=format)
            plt.close(fig)

            # fig, ax = plt.subplots(1, 1)
            # # ax = D_without_anomalies.plot(x='id', y=c, kind='scatter', ax=ax,s=10, alpha=0.3, color='#a8508b', linewidth='0.5')
            # ax.scatter(D_without_anomalies['id'],D_without_anomalies[c],s=10, alpha=0.3, color='#a8508b', linewidth='0.5')
            # # ax = D_with_anomalies.plot(x='id', y=c, kind='scatter', ax=ax,s=20, alpha=0.7, color='#cc2900', linewidth='0.5')
            # ax.scatter(D_with_anomalies['id'],D_with_anomalies[c],s=20, alpha=0.7, color='#cc2900', linewidth='0.5')
            # ax.set_xlabel(c)
            # ax.set(yscale="log")
            # ax.set_title("Scatter plot of the " + c)
            # fig.savefig(path_figure + "/scatter_plot_" + c + "_kcores.pdf", format=format)
            # plt.close(fig)

    D.drop(columns=['id'], inplace=True)

###########################################
#       KCliques Properties               #
###########################################


def preprocess_kcliques(storage_kcliques):

    # TODO : MUST ADAPT TO KCLIQUES
    '''

    :param storage_kcliques:
    :return:
    '''
    ###############
    gc.disable()
    ###############
    nodes_cores_stats = {}
    with open(storage_kcliques + "postprocess_kcliques.scf", 'rb') as ipt:
        unpacker = msgpack.Unpacker(ipt, use_list=False)
        for i in unpacker:
            k = i[0]
            for j in i[1]:
                t0, t1, u = j
                if u not in nodes_cores_stats:
                    nodes_cores_stats[u] = {k: [t0, t1]}
                elif k not in nodes_cores_stats[u]:
                    nodes_cores_stats[u][k] = [t0, t1]
                else:
                    nodes_cores_stats[u][k] += [t0, t1]
    ###############
    gc.enable()
    ###############
    return nodes_cores_stats



def preprocess_kcliques_into_dataframe(nodes_cores_stats):
    # TODO : ADAPT FOR KCLIQUES
    # Create a DataFrame with node's index (anomalous and none anomalous),
    # columns : (k, frequency in kcores, surface per cores, intervals inter core ?, interval in core ?)
    names = []
    dict_data = {'2-core': [],
                 '3-core': [],
                 '4-core': [],
                 'frequency 2-core': [],
                 'frequency 3-core': [],
                 'frequency 4-core': [],
                 'frequency': [],
                 'surface 2-core': [],
                 'surface 3-core': [],
                 'surface 4-core': [],
                 'surface': [],
                 # 'time between core':[] ?
                 }
    for u in nodes_cores_stats:
        names.append(u)
        if 2 in nodes_cores_stats[u]:
            dict_data['2-core'].append(1)
            dict_data['frequency 2-core'].append(len(nodes_cores_stats[u][2]) / 2)
            dict_data['surface 2-core'].append(
                sum([t1 - t0 for t0, t1 in zip(nodes_cores_stats[u][2][::2], nodes_cores_stats[u][2][1::2])]))
        else:
            dict_data['2-core'].append(0)
            dict_data['frequency 2-core'].append(0)
            dict_data['surface 2-core'].append(0)

        if 3 in nodes_cores_stats[u]:
            dict_data['3-core'].append(1)
            dict_data['frequency 3-core'].append(len(nodes_cores_stats[u][3]) / 2)
            dict_data['surface 3-core'].append(
                sum([t1 - t0 for t0, t1 in zip(nodes_cores_stats[u][3][::2], nodes_cores_stats[u][3][1::2])]))

        else:
            dict_data['3-core'].append(0)
            dict_data['frequency 3-core'].append(0)
            dict_data['surface 3-core'].append(0)

        if 4 in nodes_cores_stats[u]:
            dict_data['4-core'].append(1)
            dict_data['frequency 4-core'].append(len(nodes_cores_stats[u][4]) / 2)
            dict_data['surface 4-core'].append(
                sum([t1 - t0 for t0, t1 in zip(nodes_cores_stats[u][4][::2], nodes_cores_stats[u][4][1::2])]))
        else:
            dict_data['4-core'].append(0)
            dict_data['frequency 4-core'].append(0)
            dict_data['surface 4-core'].append(0)

        dict_data['frequency'].append(
            dict_data['frequency 2-core'][-1] + dict_data['frequency 3-core'][-1] + dict_data['frequency 4-core'][-1])
        dict_data['surface'].append(
            dict_data['surface 2-core'][-1] + dict_data['surface 3-core'][-1] + dict_data['surface 4-core'][-1])
    Index = pd.Index(names)
    D = pd.DataFrame(dict_data, index=Index)
    return D


########################################
#       GENERAL PROPERTIES ANALYSIS    #
########################################

def stats_condensation_graph(figures_path, stats_path,stream_nodes_path,stream_links_path,dag_path):
    S = ss.read_stream_graph(stream_links_path,stream_nodes_path)
    print("|V| = n : ",len(S.nodes))
    print("|W| = N : ",sum([len(v)/2 for v in S.node_presence]))
    print("|Ê| = m : ",len(S.links))
    print("|E| = M : ",sum([len(v)/2 for v in S.link_presence]))
    S_degrees = Counter()
    for l in S.links:
        S_degrees[l[0]]+=1
        S_degrees[l[1]]+=1
    print("S mean degree :",sum(S_degrees.values())/len(S.nodes))
    G = sd.read_global_dag(dag_path)
    print("|n_c| : ",len(G.c_nodes))
    print("|m_c| : ",len(G.c_links))
    G_degrees = Counter()
    for l in G.c_links:
        G_degrees[l[0]]+=1
        G_degrees[l[1]]+=1
    print("Gc Mean degree : ",sum(G_degrees.values())/len(G.c_nodes))

    return


if __name__ == '__main__':

    # __directory__ = "/home/leo/Dev/Data_Stream/2018/06/"
    # __file__ = "20180604"
    # __directory__ = "/home/leo/Dev/Data_Stream/2017/07/"
    # __file__ = "20170728"

    # __directory__ = "/home/leo/Dev/Data_Stream/2018/06/"  # LIP6
    # __file__ = "20180604"

    __directory__ = "/home/leo/Dev/Data_Stream/2018/04/"  # LIP6
    __filedate__ = "20180418"

    # __directory__ = "/home/leo/Dev/Data_Stream/Enron/"
    # __file__ = "enron"

    __directory__ = "/home/leo/Dev/Data_Stream/Socio_Patterns/High_School_2013/"
    __filedate__ = "High_School_2013"
    # __directory__ = "/home/leo/Dev/Data_Stream/Socio_Patterns/Workplace/"
    # __file__ = "Workplace"

    # __directory__ = "/home/leo/Dev/Data_Stream/Crawdad/Infocom/"
    # __file__ = "infocom"
    # __directory__ = "/home/leo/Dev/Data_Stream/Crawdad/Rollernet/"
    # __file__ = "rollernet"


    # __directory__ = "/home/leo/Dev/Data_Stream/Bitcoin/bitcoin_otc/"
    # __file__ = "bitcoin_otc"

    # __directory__ = "/home/leo/Dev/Data_Stream/Bitcoin/bitcoin_alpha/"
    # __file__ = "bitcoin_alpha"

    # __directory__ = "/home/leo/Dev/Data_Stream/CollegeMsg/"
    # __file__ = "CollegeMsg"

    # __directory__ = "/home/leo/Dev/Data_Stream/parse_net/11September/"
    # __file__ = "Days"

    # __directory__ = "/home/leo/Dev/Data_Stream/example_paper/"
    # __file__ = "example_paper"

    __directory__ = "/home/leo/Dev/Data_Stream/ITS/"
    __filedate__ = "Suricata"

    wcc_storage_path = __directory__ + __filedate__ + "_wcc_storage/wcc.scf"
    scc_storage_path = __directory__ + __filedate__ + "_scc_storage/scc.scf"
    kcores_storage_path = __directory__ + __filedate__ + "_kcores_storage/postprocess_kcores.scf"
    kcliques_storage_path = __directory__ + __filedate__ + "_kcliques_storage/postprocess_kcliques.scf"
    stream_nodes_path = __directory__ + __filedate__ +"_nodes.sg"
    stream_links_path = __directory__+__filedate__+"_links.sg"
    dag_path = __directory__+__filedate__+"_scc_dag_storage/scc_dag_with_links.scf"
    figures_path = __directory__ + "figures/"
    stats_path = __directory__ + "stats/"
    preprocess_data_path = __directory__ + "preprocess_data/"


    pathlib.Path(stats_path).mkdir(parents=True, exist_ok=True)
    pathlib.Path(preprocess_data_path).mkdir(parents=True, exist_ok=True)

    D = preprocess_wcc_into_dataframe(wcc_storage_path, preprocess_data_path)
    D = preprocess_scc_into_dataframe(scc_storage_path, preprocess_data_path)
    D = preprocess_kcores_into_dataframe(kcores_storage_path, preprocess_data_path)
    # D = preprocess_kcliques_into_dataframe(kcliques_storage_path, preprocess_data_path)


    stats_wcc(figures_path, stats_path, storage_path=preprocess_data_path)
    stats_scc(figures_path, stats_path, storage_path=preprocess_data_path)
    stats_kcores(figures_path,stats_path,storage_path=preprocess_data_path)
    # stats_condensation_graph(figures_path, stats_path,stream_nodes_path,stream_links_path,dag_path)




    # if preprocess_data_path:
    #     D = pd.read_msgpack(preprocess_data_path + "df_wcc.mspk")
    #     with open(preprocess_data_path + "anomalies_in_wcc.mspk", 'rb') as input:
    #         wcc_with_anomalies = set([i.decode("utf-8") for i in msgpack.load(input)])
    #
    #     D['id'] = [int(i) for i in D.index]
    #     D_without_anomalies = D.loc[list(set(D.index) - wcc_with_anomalies)]
    #     D_with_anomalies = D.loc[list(wcc_with_anomalies)]
    #     D['cat'] = ['anomalous' if i in wcc_with_anomalies else 'normal' for i in D.index]
    #     D['cat'] = D['cat'].astype('category')
    #
    # D_wcc = pd.read_msgpack(preprocess_data + "df_wcc.mspk")
    # D_scc = pd.read_msgpack(preprocess_data + "df_scc.mspk")
    # g = sns.heatmap(D_wcc.corr(), cmap=cmap,
    #                 square=True, linewidths=.5)