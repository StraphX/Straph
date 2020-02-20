import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

import straph.dag.condensation_dag as sd
import straph.stream as ss


def stats_condensation_graph(list_directory,list_file):
    dict_data = {'n':[],
                 'N':[],
                 'm':[],
                 'M':[],
                 'stream_mean_degree':[],
                 'n_c':[],
                 'm_c':[],
                 'dag_mean_degree':[]}
    for d,f in zip(list_directory,list_file):
        print('\n',d+f)
        stream_links_path = d+f+"_links.sg"
        stream_nodes_path = d+f+"_nodes.sg"
        dag_path = d+f+"_scc_dag_storage/scc_dag_with_links.scf"
        S = ss.read_stream_graph(stream_links_path,stream_nodes_path)
        S_degrees = Counter()
        for l in S.links:
            S_degrees[l[0]]+=1
            S_degrees[l[1]]+=1

        G = sd.read_global_dag(dag_path)
        G_degrees = Counter()
        for l in G.c_links:
            G_degrees[l[0]]+=1
            G_degrees[l[1]]+=1

        dict_data['n'].append(len(S.nodes))
        dict_data['N'].append(sum([len(v)/2 for v in S.node_presence]))
        dict_data['m'].append(len(S.links))
        dict_data['M'].append(sum([len(v)/2 for v in S.link_presence]))
        dict_data['stream_mean_degree'].append(sum(S_degrees.values())/len(S.nodes))
        dict_data['n_c'].append(len(G.c_nodes))
        dict_data['m_c'].append(len(G.c_links))
        dict_data['dag_mean_degree'].append(sum(G_degrees.values())/len(G.c_nodes))
        print("|V| = n : ",len(S.nodes))
        print("|W| = N : ",sum([len(v)/2 for v in S.node_presence]))
        print("|Ê| = m : ",len(S.links))
        print("|E| = M : ",sum([len(v)/2 for v in S.link_presence]))
        print("S mean degree :",sum(S_degrees.values())/len(S.nodes))
        print("|n_c| : ",len(G.c_nodes))
        print("|m_c| : ",len(G.c_links))
        print("Gc Mean degree : ",sum(G_degrees.values())/len(G.c_nodes))
    D = pd.DataFrame(dict_data,index=pd.Index(list_file))
    D.sort_values(by='M',inplace=True)
    print(D)
    D[['N','M','n_c','m_c']].plot(logy=True)
    plt.xlabel('Datasets')
    plt.tight_layout()
    plt.show()

    fig = plt.figure()
    y = D['n_c']+D['m_c']
    x = D['N']+D['M']
    print("x :",x)
    print("y :",y)
    ids = x.argsort()
    sorted_x = x[ids]
    sorted_y = y[ids]
    print("sorted x :",sorted_x)
    print("sorted_y :",sorted_y)
    plt.plot(sorted_x,sorted_y)
    plt.title('Size of Condensation (m_c+n_c) in function of size of Stream (M+N)')
    plt.xlabel('N+M')
    plt.ylabel('n_c+m_c')
    plt.xscale('log')
    plt.yscale('log')
    plt.tight_layout()
    plt.show()

    fig = plt.figure()
    x = D['n_c']
    y = D['m_c']
    ids = x.argsort()
    sorted_x = x[ids]
    sorted_y = y[ids]
    plt.plot(sorted_x,sorted_y)
    plt.title('Condensation Links (m_c) in function of Condensation Nodes (n_c)')
    plt.xlabel('n_c')
    plt.ylabel('m_c')
    plt.xscale('log')
    plt.yscale('log')
    plt.tight_layout()
    plt.show()

    return







if __name__ == '__main__':

    list_directory = ["/home/leo/Dev/Data_Stream/2018/04/",
                      "/home/leo/Dev/Data_Stream/Socio_Patterns/High_School_2012/",
                      "/home/leo/Dev/Data_Stream/Socio_Patterns/High_School_2013/",
                      "/home/leo/Dev/Data_Stream/Socio_Patterns/Workplace/",
                      "/home/leo/Dev/Data_Stream/Crawdad/Infocom/",
                      "/home/leo/Dev/Data_Stream/Crawdad/Rollernet/",
                      "/home/leo/Dev/Data_Stream/Bitcoin/bitcoin_otc/",
                      "/home/leo/Dev/Data_Stream/Bitcoin/bitcoin_alpha/",
                      "/home/leo/Dev/Data_Stream/CollegeMsg/",
                      "/home/leo/Dev/Data_Stream/parse_net/11September/",
                      "/home/leo/Dev/Data_Stream/askubuntu/",
                      # "/home/leo/Dev/Data_Stream/DBLP/",
                      ]
    file_list = ["20180418",
                 "High_School_2012",
                 "High_School_2013",
                 "Workplace",
                 "infocom",
                 "rollernet",
                 "bitcoin_otc",
                 "bitcoin_alpha",
                 "CollegeMsg",
                 "Days",
                 "askubuntu",
                 # "dblp"
                 ]

    stats_condensation_graph(list_directory,file_list)