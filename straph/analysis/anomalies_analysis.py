import re, glob, os, pathlib, subprocess,time, gc
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import msgpack
from collections import defaultdict

import matplotlib
import pandas as pd
import seaborn as sns
import numpy as np


#########################################################################################
#       THE FOLLOWING FUNCTIONS ARE DESIGNED FOR THE ANALYSIS OF THE MAWILAB DATASET    #
#########################################################################################


def parse_xml(path_xml, type=None):
    '''
    Parse XML containing anomalies of MAWILAB Dataset
    :param path_xml:
    :param type: Type of anomaly (optional)
    :return:
    '''
    tree = ET.parse(path_xml)
    anomalies = []
    root = tree.getroot()
    for child in root:
        if child.tag == "anomaly":
            anom = {}
            anom['type'] = child.attrib['type']
            if type and anom['type'] != type:
                continue
            anom['value'] = re.findall(r'[^\d\W]+', child.attrib['value'])[0]
            # anom['type_detect']=
            for c in child:
                if c.tag == "from":
                    anom['begin'] = c.attrib['sec']
                    begin = c.attrib['sec']
                if c.tag == "to":
                    anom['end'] = c.attrib['sec']
                    end = c.attrib['sec']
                    # print("length "+anom['type']+" :", float(end)-float(begin))
                for d in c:
                    if d.tag == "filter":
                        try:
                            anom['info'].append(d.attrib)
                        except:
                            anom['info'] = [d.attrib]
            anomalies.append(anom)
    return anomalies


def plot_period(anomalies, type=None):
    '''
    Plot anomalies in the MAWILAB Dataset
    :param anomalies:
    :param type: Type of anomaly (optional)
    :return:
    '''
    t_min = min([min(int(i['begin']), int(i['end'])) for i in anomalies if int(i['end']) != 0 and int(i['begin']) != 0])
    for i in anomalies:
        if type:
            if i['type'] != type:
                continue
        b = int(i['begin']) - t_min
        e = int(i['end']) - t_min
        if b != 0 and (e - b) <= 900:  # and len(i['info']) >= 2:
            print("period :", (e - b), " begin :", b, " end :", e)
            plt.bar((b + (e - b) / 2), 5, width=(e - b), color='#7a1f5c', alpha=0.2)
    plt.show()


def replace_ip_with_id(anomalies, path_dict):
    '''
    For replacing IP in anomalies by their ID
    :param anomalies:
    :param path_dict:
    :return:
    '''
    ######################
    gc.disable()
    ######################
    with open(path_dict + "_dict_ip_2_nodes_label.mspk", 'rb') as output:
        ip_2_nodes_label = msgpack.load(output)
    for anom in anomalies:
        for i in anom['info']:
            try:
                i['src_ip'] = ip_2_nodes_label[bytes(i['src_ip'], 'utf-8')]
            except:
                pass
            try:
                i['dst_ip'] = ip_2_nodes_label[bytes(i['dst_ip'], 'utf-8')]
            except:
                pass
    ######################
    gc.enable()
    ######################
    return anomalies


def preprocess_anomalies_nodes(anomalies):
    '''
    Preprocess nodes anomalies
    :param anomalies:
    :return:
    '''
    t_min = min([min(int(i['begin']), int(i['end']))
                 for i in anomalies if int(i['end']) != 0 and int(i['begin']) != 0])

    dict_anomalies = {}
    for i in anomalies:
        b = int(i['begin']) - t_min
        e = int(i['end']) - t_min
        if b == 0 or (e - b) > 900:
            b = -1
            e = -1
        for j in i['info']:
            if 'src_ip' in j and not 'dst_ip' in j:
                dict_anomalies[j['src_ip']] = (b, e)
            if 'dst_ip' in j and not 'src_ip' in j:
                dict_anomalies[j['dst_ip']] = (b, e)
    print("Anomalous nodes : ", dict_anomalies)
    return dict_anomalies


def preprocess_anomalies_links(anomalies):
    '''
    Preprocess links anomalies
    :param anomalies:
    :return:
    '''
    t_min = min([min(int(i['begin']), int(i['end']))
                 for i in anomalies if int(i['end']) != 0 and int(i['begin']) != 0])
    dict_anomalies = {}
    for i in anomalies:
        b = int(i['begin']) - t_min
        e = int(i['end']) - t_min
        if b == 0 or (e - b) > 900:
            b = -1
            e = -1
        for j in i['info']:
            if 'src_ip' in j and 'dst_ip' in j:
                src = j['src_ip']
                dst = j['dst_ip']
                dict_anomalies[(src, dst)] = (b, e)
                dict_anomalies[(dst, src)] = (b, e)
    print("Anomalous Links : ", dict_anomalies)
    return dict_anomalies



######################################################
#       Anomalies in WCC                             #
######################################################

def anomalies_in_wcc(anomalies, path_wcc, storage_path=None):
    '''
    Return WCC containing anomalies
    :param anomalies:
    :param path_wcc:
    :param storage_path:
    :return:
    '''
    anomalous_nodes = preprocess_anomalies_nodes(anomalies)
    anomalous_links = preprocess_anomalies_links(anomalies)
    wcc_with_anomalies = set()
    ######################
    gc.disable()
    ######################
    with open(path_wcc,'rb') as input:
        unpacker = msgpack.Unpacker(input,use_list=False)
        for i in unpacker:
            id_wcc =i[0]
            nodes = set()
            link_presence = defaultdict(list)
            nodes_anomalies_in_wcc = defaultdict(int)
            links_anomalies_in_wcc = defaultdict(int)
            times = i[1],i[2]
            for j in i[3]:
                if j[0] == 1:
                    t0, t1 = j[1], j[2]
                    u, v = j[3],j[4]
                    link_presence[(u, v)] += [t0, t1]
                    nodes.add(u)
                    nodes.add(v)
                else:
                    continue
                if u in anomalous_nodes:
                    b, e = anomalous_nodes[u]
                    if b == -1:
                        nodes_anomalies_in_wcc[u] += 1
                    elif e >= times[0] and times[1] >= b:
                        nodes_anomalies_in_wcc[u] += 1
                if v in anomalous_nodes:
                    b, e = anomalous_nodes[v]
                    if b == -1:
                        nodes_anomalies_in_wcc[v] += 1
                    elif e >= times[0] and times[1] >= b:
                        nodes_anomalies_in_wcc[v] += 1
                if (u, v) in anomalous_links:
                    b, e = anomalous_links[(u, v)]
                    if b == -1:
                        links_anomalies_in_wcc[(u, v)] += 1
                    elif e >= times[0] and times[1] >= b:
                        links_anomalies_in_wcc[(u, v)] += 1
            if len(nodes_anomalies_in_wcc) > 0 and len(links_anomalies_in_wcc) > 0:
                print("Current WCC :", id_wcc)
                print("Nb nodes anomalies :", len(nodes_anomalies_in_wcc))
                print("Nb links anomalies :", len(links_anomalies_in_wcc))
                wcc_with_anomalies.add(id_wcc)

    if storage_path:
        pathlib.Path(storage_path).mkdir(parents=True, exist_ok=True)
        with open(storage_path + "anomalies_in_wcc.mspk", 'wb') as output:
            msgpack.dump(list(wcc_with_anomalies), output)
    ######################
    gc.enable()
    ######################
    return wcc_with_anomalies


######################################################
#       Anomalies in SCC                             #
######################################################


def anomalies_in_scc(anomalies, path_scc, storage_path=None):
    '''
    Return SCC containing anomalies
    :param anomalies:
    :param path_scc:
    :param storage_path:
    :return:
    '''
    anomalous_nodes = preprocess_anomalies_nodes(anomalies)
    anomalous_links = preprocess_anomalies_links(anomalies)
    scc_with_anomalies = set()
    t_begin =time.time()
    ######################
    gc.disable()
    ######################
    with open(path_scc, 'rb') as input:
        unpacker = msgpack.Unpacker(input,use_list=False)
        for i in unpacker:
            id_wcc,id_scc = i[0],i[1]
            times = i[2],i[3]
            if id_wcc % 100000 == 0 or (id_scc % 100000 == 0 and id_scc !=0):
                print("id wcc :",id_wcc,"id scc :",id_scc)
                print("times :",times,"size scc :",len(i[4]), " t :", time.time() - t_begin)
            nodes = set()
            link_presence = defaultdict(list)
            nodes_anomalies_in_scc = defaultdict(int)
            links_anomalies_in_scc = defaultdict(int)
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
                if u in anomalous_nodes:
                    b, e = anomalous_nodes[u]
                    if b == -1:
                        nodes_anomalies_in_scc[u] += 1
                    elif e >= times[0] and times[1] >= b:
                        nodes_anomalies_in_scc[u] += 1
                if v in anomalous_nodes:
                    b, e = anomalous_nodes[v]
                    if b == -1:
                        nodes_anomalies_in_scc[v] += 1
                    elif e >= times[0] and times[1] >= b:
                        nodes_anomalies_in_scc[v] += 1
                if (u, v) in anomalous_links:
                    b, e = anomalous_links[(u, v)]
                    if b == -1:
                        links_anomalies_in_scc[(u, v)] += 1
                    elif e >= times[0] and times[1] >= b:
                        links_anomalies_in_scc[(u, v)] += 1
            if len(nodes_anomalies_in_scc) > 0 and len(links_anomalies_in_scc) > 0:
                print("Current SCC :", (id_wcc,id_scc))
                print("Nb nodes anomalies :", len(nodes_anomalies_in_scc))
                print("Nb links anomalies :", len(links_anomalies_in_scc))
                scc_with_anomalies.add((id_wcc,id_scc))

    if storage_path:
        pathlib.Path(storage_path).mkdir(parents=True, exist_ok=True)
        with open(storage_path + "anomalies_in_scc.mspk", 'wb') as output:
            msgpack.dump(list(scc_with_anomalies), output)
    ######################
    gc.enable()
    ######################
    return scc_with_anomalies




######################################################
#       Anomalies in Kcliques                        #
######################################################


def anomalies_in_kcliques(anomalies, storage_kcliques):
    # TODO : ADADT FOR KCLIQUES
    '''
    Return kcliques containing anomalies and kcliques considered as normal.
    (For a comparison/discriminating purpose)
    :param anomalies:
    :param storage_kcliques:
    :return:
    '''
    anomalous_nodes = preprocess_anomalies_nodes(anomalies)
    print("nb anomalous nodes :",len(anomalous_nodes))
    kcliques_with_anomalies = defaultdict(list)
    kcliques_without_anomalies = defaultdict(list)
    ######################
    gc.disable()
    ######################
    with open(storage_kcliques, 'rb') as ipt:
        unpacker = msgpack.Unpacker(ipt,use_list=False)
        for i in unpacker:
            k = i[0]
            print(" k : ",k)
            for j in i[1]:
                t0,t1,u = j
                if u in anomalous_nodes:
                    b, e = anomalous_nodes[u]
                    if b == -1:
                        kcliques_with_anomalies[k].append(j)
                    elif e >= t0 and t1 >= b:
                        kcliques_with_anomalies[k].append(j)
                else:
                    kcliques_without_anomalies[k].append(j)
            if len(kcliques_with_anomalies[k]) > 0:
                print("Current K : ",k)
                print("Nb nodes anomalies : ", len(kcliques_with_anomalies[k]))
                print("Nb nodes without anomalies : ",len(kcliques_without_anomalies[k]))
#     if storage_path:
#         pathlib.Path(storage_path).mkdir(parents=True, exist_ok=True)
#         with open(storage_path + "anomalies_in_scc.mspk", 'wb') as output:
#             msgpack.dump(list(scc_with_anomalies), output)
    ######################
    gc.enable()
    ######################
    return kcliques_with_anomalies,kcliques_without_anomalies


######################################################
#       Anomalies in Kcores                          #
######################################################

def anomalies_in_kcores(anomalies, storage_kcores):
    #TODO
    '''
    Return kcores containing anomalies and kcores considered as 'normal'.
    (For a comparison/discriminating purpose)
    :param anomalies:
    :param storage_kcores:
    :return:
    '''
    anomalous_nodes = preprocess_anomalies_nodes(anomalies)
    print("nb anomalous nodes :",len(anomalous_nodes))
    kcores_with_anomalies = defaultdict(list)
    kcores_without_anomalies = defaultdict(list)
    ######################
    gc.disable()
    ######################
    with open(storage_kcores, 'rb') as ipt:
        unpacker = msgpack.Unpacker(ipt,use_list=False)
        for i in unpacker:
            k = i[0]
            print(" k : ",k)
            for j in i[1]:
                t0,t1,u = j
                if u in anomalous_nodes:
                    b, e = anomalous_nodes[u]
                    if b == -1:
                        kcores_with_anomalies[k].append(j)
                    elif e >= t0 and t1 >= b:
                        kcores_with_anomalies[k].append(j)
                else:
                    kcores_without_anomalies[k].append(j)
            if len(kcores_with_anomalies[k]) > 0:
                print("Current K : ",k)
                print("Nb nodes anomalies : ", len(kcores_with_anomalies[k]))
                print("Nb nodes without anomalies : ",len(kcores_without_anomalies[k]))
#     if storage_path:
#         pathlib.Path(storage_path).mkdir(parents=True, exist_ok=True)
#         with open(storage_path + "anomalies_in_scc.mspk", 'wb') as output:
#             msgpack.dump(list(scc_with_anomalies), output)
    ######################
    gc.enable()
    ######################
    return kcores_with_anomalies,kcores_without_anomalies



def plot_lines_kcores(node_core_info):
    seg = []
    for u in node_core_info:
        t0, t1, n = u
        seg.append(((t0, n), (t1, n)))
    return seg


def plot_kcores(kcores_with_anomalies, kcores_without_anomalies, title=None, legend=None, saving_path=None,
                format='pdf'):
    #TODO: Keep useful code
    nodes = set()
    for v in kcores_with_anomalies.values():
        for i in v:
            nodes.add(i[2])
    for v in kcores_without_anomalies.values():
        for i in v:
            nodes.add(i[2])
    t_min = 0
    t_max = 900
    print("T min :", t_min)
    print("T max :", t_max)

    dict_color_anom = {2: "#ffcc00",  # Yellow
                       3: "#e68a00",  # Orange
                       4: "#cc0000"}  # Rouge
    dict_color_normal = {2: "#8080ff",  # Light blue
                         3: "#0033cc",  # Medium blue
                         4: "#004466"}  # Strong blue
    colors_anom = []
    segs_anom = []
    colors_normal = []
    segs_normal = []
    for k, v in kcores_with_anomalies.items():
        if k >= 2:
            seg = plot_lines_kcores(v)
            segs_anom += seg
            colors_anom += [dict_color_anom[k]] * len(seg)
    for k, v in kcores_without_anomalies.items():
        if k >= 2:
            seg = plot_lines_kcores(v)
            segs_normal += seg
            colors_normal += [dict_color_normal[k]] * len(seg)

    print("len colors normal : ", len(colors_normal))
    print("len segs normal : ", len(segs_normal))
    print("len colors anomalous : ", len(colors_anom))
    print("len segs anomalous : ", len(segs_anom))
    print("Segs calculated")
    fig, ax = plt.subplots(1, 1)
    line_coll = matplotlib.collections.LineCollection(np.array(segs_normal),
                                                      colors=colors_normal, linewidths=[5 * 10 ** (-3)])
    line_coll.set_alpha(0.7)
    ax.add_collection(line_coll)

    line_coll = matplotlib.collections.LineCollection(np.array(segs_anom),
                                                      colors=colors_anom, linewidths=[5 * 10 ** (-1)])
    line_coll.set_alpha(1)
    ax.add_collection(line_coll)
    print("Collections created and added")
    ax.set_xlim(t_min, t_max)
    ax.set_ylim(0, max(nodes))
    if legend:
        plt.ylabel("Nodes", fontname='Ubuntu', fontsize=12, color='#666699')
        plt.xlabel("t", fontname='Ubuntu', fontsize=12, color='#476b6b')
        list_legend = []
        for k, current_color in dict_color_anom.items():
            list_legend.append(matplotlib.patches.Patch(color=current_color, label=str(k) + "-Cores Anomalous"))
        for k, current_color in dict_color_normal.items():
            list_legend.append(matplotlib.patches.Patch(color=current_color, label=str(k) + "-Cores"))
        plt.legend(handles=list_legend, loc='upper left', fancybox=True, prop={'size': 6})
    if title:
        ax.set_title(title, fontname='Ubuntu', fontsize=14)
    for spine in plt.gca().spines.values():
        spine.set_visible(False)
    plt.tick_params(top=False, bottom=True, right=False, left=True, labelleft=True, labelbottom=True)
    plt.tight_layout()
    if saving_path:
        fig.savefig(saving_path + "." + format, format=format, dpi=2400)



###############################################
#       GENERAL ANOMALIES ANALYSIS            #
###############################################
# TODO : Regroup functions, anomalies_in_+"*"


if __name__ == '__main__':
    # notice = parse_xml("/home/leo/Dev/data/2017/20170121_notice.xml")
    # print("nb notices :",len(notice))
    # for i in range(len(notice)):
    #     print("notice "+str(i)+" :",notice[i])
    __directory__ = "/home/leo/Dev/Data_Stream/2018/04/"
    __filedate__ = "20180418"

    wcc_storage_path = __directory__ + __filedate__ + "_wcc_storage/wcc.scf"
    scc_storage_path = __directory__ + __filedate__ + "_scc_storage/scc.scf"
    kcores_storage_path = __directory__ + __filedate__ + "_kcores_storage/kcores.scf"
    preprocess_data_path = __directory__ + "preprocess_data/"

    anomalies = parse_xml(__directory__ + __filedate__ + "_anomalous_suspicious.xml", type='anomalous')
    anomalies = replace_ip_with_id(anomalies, __directory__ + __filedate__)
    print("nb anomalies :", len(anomalies))
    # plot_period(anomalies,type='anomalous')
    for i in range(len(anomalies)):
        print("anomalie " + str(i) + " :", anomalies[i])


    # wcc_with_anomalies = anomalies_in_wcc(anomalies, wcc_storage_path,preprocess_data_path)
    scc_with_anomalies = anomalies_in_scc(anomalies, scc_storage_path, preprocess_data_path)
    # kcores_with_anomalies = anomalies_in_kcores(anomalies,kcores_storage_path) #TODO