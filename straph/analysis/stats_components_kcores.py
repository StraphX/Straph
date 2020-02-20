import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import json
import re
from collections import Counter


def analysis_components(S,components):

    ########################################
    # Nb de liens:noeuds par composante
    # Surfaces composantes
    ########################################
    # nb_isolated_nodes = [len(cc) for cc in components if len(cc) ==1]
    # print("Nb de noeud isolés :",sum(nb_isolated_nodes))
    # size_components = [len(cc) for cc in components if len(cc) >1]
    # print("Nb de composantes :",len(size_components))
    # sns.countplot(size_components)
    # plt.title("Size of Connected Components")
    # plt.show()
    # degrees = S.degrees()
    # sns.countplot(list(degrees.values()))
    # plt.title("Degrees")
    # plt.show()
    surfaces =[]
    for cc in components:
        surfaces.append(analysis_intra_component(S, cc))

    return

def analysis_intra_component(S, component):
    surface = sum([n[2]-n[1] for n in component])
    print("Surface component :",surface)


    nodes = set([n[0] for n in component])
    ########################################
    # columns:
    #   redondance noeud
    #   Nb de liens par noeud
    #
    ########################################
    if len(nodes) <9:
        return 0

    cnt_redon = Counter()
    cnt_surface_per_node = Counter()
    for n in component:
        cnt_redon[n[0]]+=1
        cnt_surface_per_node[n[0]] += (n[2]-n[1])
    print("redondance :",cnt_redon)
    plt.figure()
    sns.countplot(list(cnt_redon.values()))
    plt.title("Redondance nodes in component")
    plt.figure()
    sns.distplot(list(cnt_surface_per_node.values()))
    plt.title("Surface per node in component")
    plt.show()
    return surface


def analysis_kcores(S,k_cores):
    #############################################
    # Implication des noeuds dans chaque core
    # Nb de liens/noeuds par k_cores
    # Surfaces kcores
    #############################################

    # size_kcores = [len(k_cores[k]) for k in k_cores]
    # print("Nb de composantes :",len(k_cores))
    # sns.countplot(size_kcores)
    # plt.title("Size of K_cores")
    # plt.show()

    surfaces = []
    for k in k_cores:
        print("current k :",k)
        surfaces.append(analysis_intra_kcore(S, k_cores[k],k))
    return

def analysis_intra_kcore(S, core,k):
    surface = sum([n[2]-n[1] for n in core])
    print("Surface core :",surface)
    nodes = set([n[0] for n in core])
    if len(nodes) <9:
        return 0
    ####################################
    # columns:
    #   redondance noeud
    #   Nb de liens par noeud
    # Evolution en fonction du temps du kcore (global)
    ####################################
    cnt_redon = Counter()
    cnt_surface_per_node = Counter()
    times_k_core = []
    for n in core:
        cnt_redon[n[0]]+=1
        cnt_surface_per_node[n[0]] += (n[2]-n[1])
        times_k_core += [n[1],n[2]]
    plt.figure()
    sns.distplot(list(cnt_redon.values()))
    plt.title("Redondance nodes in "+str(k)+"core")
    plt.figure()
    sns.distplot(list(cnt_surface_per_node.values()))
    plt.title("Surface per node in "+str(k)+"core")

    plt.figure()
    sns.distplot(times_k_core)
    plt.title("Surface "+str(k)+"core")

    plt.show()

    return surface

def stats_nodes_links(path_directory):
    dict_data= {'nb_nodes' : [],
                'nb_links':[],
                'n_V':[],
                'n_E':[]}
    tw = []
    # for filename in glob.glob(__directory__):
    filename = path_directory+"SG.json"
    f = open(filename,'r')
    for i in f:
        js = json.loads(i)
        tw.append(float(js["time_window"][0]))
        dict_data['nb_nodes'].append(int(js['nb_nodes']))
        dict_data['nb_links'].append(int(js['nb_links']))
        dict_data['n_V'].append(float(js['n_V']))
        dict_data['n_E'].append(float(js['n_E']))
    Time_Index = pd.Index(tw)
    D = pd.DataFrame(dict_data,index=Time_Index)
    print(D.describe())
    for c in D.columns:
        pd.to_numeric(D[c])
    ax = D.plot()
    fig = ax.get_figure()
    fig.savefig(path_directory+"/figures/plot_nodes_links.eps",format='eps')


def stats_SCC(path_directory):
    dict_data= {'n_scc' : [],
                'max_size_scc':[],
                'mean_size_scc':[]}
    tw = []
    # for filename in glob.glob(__directory__):
    filename = path_directory+"SG.json"
    f = open(filename,'r')
    for i in f:
        js = json.loads(i)
        tw.append(float(js["time_window"][0]))
        dict_data['n_scc'].append(int(js["n_scc"]))
        dict_data['max_size_scc'].append(float(js["max_size_scc"]))
        dict_data['mean_size_scc'].append(float(js["mean_size_scc"]))
    Time_Index = pd.Index(tw)
    D = pd.DataFrame(dict_data,index=Time_Index)
    print(D.describe())
    for c in D.columns:
        pd.to_numeric(D[c])
    ax = D.plot()
    fig = ax.get_figure()
    fig.savefig(path_directory+"/figures/plot_scc.eps",format='eps')
    return


def stats_WCC(path_directory):
    dict_data= {'n_wcc' : [],
                'max_size_wcc':[],
                'mean_size_wcc':[]}
    tw = []
    # for filename in glob.glob(__directory__):
    filename = path_directory+"SG.json"
    f = open(filename,'r')
    for i in f:
        js = json.loads(i)
        tw.append(float(js["time_window"][0]))
        dict_data['n_wcc'].append(int(js["n_wcc"]))
        dict_data['max_size_wcc'].append(float(js["max_size_wcc"]))
        dict_data['mean_size_wcc'].append(float(js["mean_size_wcc"]))
    Time_Index = pd.Index(tw)
    D = pd.DataFrame(dict_data,index=Time_Index)
    print(D.describe())
    for c in D.columns:
        pd.to_numeric(D[c])
    ax = D.plot()
    fig = ax.get_figure()
    fig.savefig(path_directory+"/figures/plot_wcc.eps",format='eps')

def stats_Kcores(path_directory):
    tw = []
    # for filename in glob.glob(__directory__):
    filename = path_directory+"SG.json"
    f = open(filename,'r')
    set_keys= set()
    for i in f:
        js = json.loads(i)
        for k in js.keys():
            if re.findall('core',k):
                set_keys.add(k)
    print("l_keys :",set_keys)
    dict_data = {k: [] for k in set_keys}
    f = open(filename,'r')
    for i in f:
        js = json.loads(i)
        tw.append(float(js["time_window"][0]))
        for k in set_keys:
            if k in js.keys():
                dict_data[k].append(float(js[k]))
            else:
                dict_data[k].append(0)
    Time_Index = pd.Index(tw)
    D = pd.DataFrame(dict_data,index=Time_Index)
    print(D.describe())
    for c in D.columns:
        pd.to_numeric(D[c])
    l = set(D.columns)
    k=1
    while l:
        column_to_plot = []
        for c in l:
            if re.findall(str(k),c):
                column_to_plot.append(c)
        l-=set(column_to_plot)
        if column_to_plot:
            ax = D[column_to_plot].plot()
            fig = ax.get_figure()
            fig.savefig(path_directory+"/figures/plot_"+str(k)+"cores.eps",format='eps')
        k+=1
    return


def stats_degree(path_directory):
    tw = []
    # for filename in glob.glob(__directory__):
    filename = path_directory+"SG.json"
    f = open(filename,'r')
    set_keys= set()
    for i in f:
        js = json.loads(i)
        for k in js.keys():
            if re.findall('degree',k):
                set_keys.add(k)
    print("l_keys :",set_keys)
    dict_data = {k: [] for k in set_keys}
    f = open(filename,'r')
    for i in f:
        js = json.loads(i)
        tw.append(float(js["time_window"][0]))
        for k in set_keys:
            dict_data[k].append(float(js[k]))
    Time_Index = pd.Index(tw)
    D = pd.DataFrame(dict_data,index=Time_Index)
    print(D.describe())
    for c in D.columns:
        pd.to_numeric(D[c])
    ax = D.plot()
    fig = ax.get_figure()
    fig.savefig(path_directory+"/figures/plot_degree.eps",format='eps')
    return


#####################################################################
#  Diviser l'analyse en statistiques de deux différentes natures
#       - Statistiques intra composantes
#       - Statistiques sur l'ensemble des composantes
#
# Analyse entremélée :
#       - Evaluer le nombre de k_cores dans chaque composante
#       - surface des kcores par rapport à celle de leur composante
#
#####################################################################


if __name__ == '__main__':
    path = "/data/rannou/data/2017/07/protocols/stats/TCP/"
    stats_nodes_links(path)
    stats_SCC(path)
    stats_Kcores(path)
    stats_WCC(path)
    stats_degree(path)