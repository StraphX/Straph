# Copyright (C) 2017-2021 Léo Rannou - Sorbonne Université/LIP6 - Thales
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import csv
import dateutil.parser as du
import dpkt
import json
import math
import os
import socket
import time
from collections import defaultdict
from sortedcollections import SortedSet
from tqdm import tqdm

from straph import stream as sg

# TODO : parse PCAP (adpat pcap_to_csv and shit), see pcap_reader.
# TODO : parse net, to finish (for Pajek datasets).

__nb_2_protocol__ = {0: 'IPv6_HbH',  # IPv6 Hop by Hop
                     1: 'ICMP',  # Internet Control Message
                     2: 'IGMP',  # Internet Group Management
                     3: 'GGP',  #  Gateway-to-Gateway
                     4: 'IPv4',
                     5: 'ST2',  #  Sream v2 aka "IPv5"
                     6: 'TCP',  # Transmission Control
                     7: 'CBT',
                     8: 'EGP',  # Exterior Gateway Protocol
                     17: 'UDP',  # User Datagram
                     41: 'IPv6',
                     43: 'IPv6_Route',  #  Routing header for IPv6
                     47: 'GRE',  # Generic Routing encapsulation
                     50: 'ESP',  # Encap Security Payload
                     51: 'AH',  # Authenfication Header
                     58: 'IPv6_ICMP',  # IPV6 ICMP
                     103: 'PIM',  # Protocol Independent Multicast
                     }


def inet_to_str(inet):
    """
    Convert inet object to a string

    :param inet: inet network address
    :return: str: Printable/readable IP address
    """
    # First try ipv4 and then ipv6
    try:
        return socket.inet_ntop(socket.AF_INET, inet)
    except ValueError:
        return socket.inet_ntop(socket.AF_INET6, inet)


def datetime_to_timestamp(s):
    return du.parse(s).timestamp()


def pcap_to_csv(file_input, destination, protocol=None):
    """
    Transform a pcap file to a csv

    :param file_input:
    :param destination:
    :param protocol:
    :return:
    """
    counter = 0
    dict_ip = defaultdict(lambda: len(dict_ip))
    dict_label = {}
    if protocol:
        protocol = [key for key, value in __nb_2_protocol__.items() if value == protocol][0]
    print("protocol :", protocol)
    with open(destination, 'w') as output, open(file_input, 'rb') as input:
        writer = csv.writer(output, delimiter=';')
        writer.writerow(["time", "src", "dst", "protocol", "len", "src_port", "dst_port"])
        for ts, pkt in tqdm(dpkt.pcap.Reader(input)):
            eth = dpkt.ethernet.Ethernet(pkt)
            # if counter == 150000000000000000000000000000:
            #     print("Time end dkpt :", time.time() - start)
            #     break
            if counter == 0:
                t0 = ts
            counter += 1
            if isinstance(eth.data, bytes):
                continue
            ip = eth.data
            if ip.src is not None:
                ip_src = inet_to_str(ip.src).encode()
            else:
                continue
            ip_dst = inet_to_str(ip.dst).encode()

            id_src = dict_ip[ip_src]
            id_dst = dict_ip[ip_dst]
            dict_label[id_src] = ip_src
            dict_label[id_dst] = ip_dst

            # if counter % 1000000 == 0:
            #     print("Counter dkpt:", counter, "time dkpt:", time.time() - start)
            # We ignore 'ICMP' protocols, ICMP scan useless
            if ip.p == 1:
                continue
            if protocol:
                if ip.p != protocol:
                    continue
            if isinstance(eth.data, dpkt.ip6.IP6):
                # print("IPv6 !")
                len_pckt = ip.plen
            else:
                len_pckt = ip.len
            if isinstance(ip.data, dpkt.tcp.TCP) or isinstance(ip.data, dpkt.udp.UDP):
                tcp = ip.data
                src_port = tcp.sport
                dst_port = tcp.dport
            else:
                src_port = None
                dst_port = None
            # print("ip src :", ip_src, " ip dst :", ip_dst, " ip protocol :",
            #       protocol, " ip len :", len_pckt,"src port :", tcp.sport, "dest port :", tcp.dport)
            writer.writerow([round(ts - t0, 6), id_src, id_dst, ip.p, len_pckt, src_port, dst_port])


def parse_net(input_file, output_file_nodes, output_file_links, link_duration=1):
    """
    A Stream Graph reader for dataset issued by Pajek

    Format of interactions : .net
    :param input_file:
    :param output_file_nodes:
    :param output_file_links:
    :param link_duration:
    :return:
    """
    E = defaultdict(list)
    W = defaultdict(list)
    type_node = None
    with open(input_file, 'r') as input_file:
        for line in input_file:
            l = line.strip().split()
            if l[0] == '*Vertices':
                type_node = True
                continue
            if l[0] == '*Edges':
                type_node = False
                continue
            if type_node:
                continue
            else:
                u, v = int(l[0]), int(l[1])
                e = (u, v)
                if u == v:
                    # SELF LOOP : we ignore it
                    continue
                t = l[3].strip('[').strip(']').split(',')

            for current_time in t:
                current_time = int(current_time)
                if e in E and E[e][-1] >= current_time:
                    # print("Extend Link Presence")
                    E[e][-1] = max(E[e][-1], current_time + link_duration)
                else:
                    E[e] += [current_time, current_time + link_duration]

                if u in W and W[u][-1] >= current_time:
                    # print("Extend Node Presence")
                    W[u][-1] = max(W[u][-1], current_time + link_duration)
                else:
                    W[u] += [current_time, current_time + link_duration]

                if v in W and W[v][-1] >= current_time:
                    # print("Extend Node Presence")
                    W[v][-1] = max(W[v][-1], current_time + link_duration)
                else:
                    W[v] += [current_time, current_time + link_duration]

    with open(output_file_links, 'w') as output_file:
        for k, v in E.items():
            output_file.write(str(k[0]) + " " + str(k[1]) + " ")
            for t in v:
                output_file.write(str(t) + " ")
            output_file.write("\n")
    with open(output_file_nodes, 'w') as output_file:
        for k, v in W.items():
            output_file.write(str(k) + " ")
            for t in v:
                output_file.write(str(t) + " ")
            output_file.write("\n")
    return None


def parse_csv(input_file, entry_format, **kwargs):
    """
    Reader for .csv files

    :param input_file:
    :param entry_format:
    :param kwargs:
    :return:
    """
    # Convert entry format
    t_pos, b_pos, e_pos, link_duration_pos = None, None, None, None
    if len(entry_format) == 3:
        (t_pos, u_pos, v_pos) = entry_format['t_pos'], entry_format['u_pos'], entry_format['v_pos']
    elif len(entry_format) == 4 and 'link_duration_pos' in entry_format:
        (t_pos, link_duration_pos, u_pos, v_pos) = entry_format['t_pos'], entry_format['link_duration_pos'], \
                                                   entry_format['u_pos'], entry_format['v_pos']
    elif len(entry_format) == 4 and 'b_pos' in entry_format:
        (b_pos, e_pos, u_pos, v_pos) = entry_format['b_pos'], entry_format['e_pos'], \
                                       entry_format['u_pos'], entry_format['v_pos']
    else:
        raise TypeError("Entry format is not supported, see documentation !")
    E = defaultdict(list)
    W = defaultdict(list)
    cnt_rows = 0

    nodes_to_label = {}
    label_to_id = defaultdict(lambda: len(label_to_id))

    min_t, max_t = math.inf, -math.inf

    with open(input_file, 'r') as input_file:

        reader = csv.reader(input_file, delimiter=kwargs['delimiter'])

        if kwargs['ignore_header']:
            next(reader, None)

        if kwargs['link_duration']:
            link_duration = kwargs['link_duration']

        elif 'link_duration_pos' not in entry_format and 'b_pos' not in entry_format:
            link_duration = 0
            print("[WARNING] No link_duration provided, links durations are set to 0.")

        for line in tqdm(reader, desc='Parsing CSV', total=kwargs['nrows'] - 1):
            cnt_rows += 1

            if cnt_rows > kwargs['nrows']:
                break

            if kwargs['nodes_to_label']:
                # Convert Label to int
                u_label = line[u_pos]
                v_label = line[v_pos]
                if u_label in {'', ' '} or v_label in {'', ' '}:
                    # print(" Blank node line:",cnt_rows)
                    # print(" Content:",line)
                    continue

                if u_label == v_label:
                    # SELF LOOP : we ignore it
                    continue
                # If we haven't these label before they are assigned to len(label_to_id) = new_id
                u = label_to_id[u_label]
                v = label_to_id[v_label]
                nodes_to_label[u] = u_label
                nodes_to_label[v] = v_label
            else:
                u = int(line[u_pos])
                v = int(line[v_pos])
                if u == v:
                    # SELF LOOP : we ignore it
                    continue

            if kwargs['time_is_datetime']:
                if 't_pos' in entry_format:
                    t = datetime_to_timestamp(line[t_pos])
                elif 'b_pos' in entry_format:
                    b = datetime_to_timestamp(line[b_pos])
                    e = datetime_to_timestamp(line[e_pos])
                    link_duration = e - b
                    t = b
            else:
                if 't_pos' in entry_format:
                    t = float(line[t_pos].replace(',', ''))
                elif 'b_pos' in entry_format:
                    b = float(line[b_pos].replace(',', ''))
                    e = float(line[e_pos].replace(',', ''))
                    link_duration = e - b
                    t = b

            if 'link_duration_pos' in entry_format:
                link_duration = float(line[link_duration_pos].replace(',', ''))

            min_t, max_t = min(min_t, t), max(max_t, t + link_duration)

            if kwargs['is_directed']:
                l = (u, v)
            else:
                if (v, u) in E:
                    l = (v, u)
                else:
                    l = (u, v)

            if l in E and E[l][-1] >= t:
                E[l][-1] = max(E[l][-1], t + link_duration)
            else:
                E[l] += [t, t + link_duration]
            if kwargs['is_link_stream'] is False:
                if u in W and W[u][-1] >= t:
                    W[u][-1] = max(W[u][-1], t + link_duration)
                else:
                    W[u] += [t, t + link_duration]

                if v in W and W[v][-1] >= t:
                    W[v][-1] = max(W[v][-1], t + link_duration)
                else:
                    W[v] += [t, t + link_duration]
            else:
                W[u] = [min_t, max_t]
                W[v] = [min_t, max_t]

    if kwargs['is_link_stream'] is True:
        W = {k: [min_t, max_t] for k in W.keys()}

    if kwargs['delta']:
        delta = kwargs['delta']
        chrono = time.time()
        W, E = approximate_events(W, E, delta)
        print("\t Approximate events with delta :", delta, " in ", time.time() - chrono)

    S = sg.StreamGraph(times=[min_t, max_t],
                       nodes=list(W.keys()),
                       links=list(E.keys()),
                       node_presence=[W[k] for k in W.keys()],
                       link_presence=[E[k] for k in E.keys()],
                       node_to_label=nodes_to_label,
                       node_to_id={i: i for i in W.keys()})
    return S


def approximate_events(W, E, delta):
    """
    Approximation method reducing the number of distinct event times while preserving connectivity properties
    of the original dataset.

    :param W:
    :param E:
    :param delta:
    :return:
    """
    # Seems strange but avoid float imprecision
    event_times = sorted(set([t for np in W.values() for t in np] + [t for lp in E.values() for t in lp]))
    t_old = event_times[0]
    discretized_event_times = SortedSet()
    discretized_event_times.add(t_old)
    for t in event_times:
        if t - t_old >= delta:
            discretized_event_times.add(t)
            t_old = t

    new_W = {}
    for n, np in W.items():
        new_W[n] = []
        for t0, t1 in zip(np[::2], np[1::2]):
            assert t1 - t0 >= delta
            if t0 not in discretized_event_times:
                # # Catch time after t0 in discretized event times:
                t0 = discretized_event_times[discretized_event_times.bisect(t0)]
            if t1 not in discretized_event_times:
                # #Catch time before t1 in discretize event times:
                t1 = discretized_event_times[discretized_event_times.bisect(t1) - 1]

            # # new_W[n] += [t0, t1]
            # a, b = delta * math.ceil(t0 / delta), delta * math.floor(t1 / delta)
            #
            # if math.isclose(a, t0):
            #     a = t0
            # if math.isclose(b, t1):
            #     b = t1
            new_W[n] += [t0, t1]

    new_E = {}
    for l, lp in E.items():
        new_E[l] = []
        for t0, t1 in zip(lp[::2], lp[1::2]):
            assert t1 - t0 >= delta
            if t0 not in discretized_event_times:
                # # Catch time after t0 in discretized event times:
                t0 = discretized_event_times[discretized_event_times.bisect(t0)]
            if t1 not in discretized_event_times:
                # # Catch time before t1 in discretize event times:
                t1 = discretized_event_times[discretized_event_times.bisect(t1) - 1]

            # new_E[l] += [t0, t1]
            # a, b = delta * math.ceil(t0 / delta), delta * math.floor(t1 / delta)
            # if math.isclose(a, t0):
            #     a = t0
            # if math.isclose(b, t1):
            #     b = t1
            new_E[l] += [t0, t1]
    return new_W, new_E


def parse_json(input_file, entry_format, **kwargs):
    """
    A Stream Graph reader for JSON dataset.

    :param input_file:
    :param entry_format:
    :param kwargs:
    :return:
    """
    # Convert entry format
    u_pos, v_pos, t_pos, b_pos, e_pos, link_duration_pos = None, None, None, None, None, None
    if len(entry_format) == 3:
        (t_pos, u_pos, v_pos) = entry_format['t_pos'], entry_format['u_pos'], entry_format['v_pos']
    elif len(entry_format) == 4 and 'link_duration_pos' in entry_format:
        (t_pos, link_duration_pos, u_pos, v_pos) = entry_format['t_pos'], entry_format['link_duration_pos'], \
                                                   entry_format['u_pos'], entry_format['v_pos']
    elif len(entry_format) == 4 and 'b_pos' in entry_format:
        (b_pos, e_pos, u_pos, v_pos) = entry_format['b_pos'], entry_format['e_pos'], \
                                       entry_format['u_pos'], entry_format['v_pos']
    E = defaultdict(list)
    W = defaultdict(list)
    cnt_rows = 0
    id_to_label = {}
    label_to_id = defaultdict(lambda: len(label_to_id))
    min_t, max_t = math.inf, -math.inf

    with open(input_file, 'r') as input_file:
        reader = json.load(input_file)
        for line in tqdm(reader, desc="Parsing JSON", total=kwargs['nrows']):
            cnt_rows += 1

            # if cnt_rows % 100000 == 0:
            #     print((cnt_rows / kwargs['nrows']) * 100, "% loaded")

            if cnt_rows > kwargs['nrows']:
                break
            if kwargs['nodes_to_label']:
                # Convert Label to int
                u_label = line[u_pos]
                v_label = line[v_pos]
                # If we haven't these label before they are assigned to len(label_to_id = new_id)
                u = label_to_id[u_label]
                v = label_to_id[v_label]
                id_to_label[u] = u_label
                id_to_label[v] = v_label
            else:
                u = int(line[u_pos])
                v = int(line[v_pos])

            if kwargs['time_is_datetime']:
                if 't_pos' in entry_format:
                    t = datetime_to_timestamp(line[t_pos])
                elif 'b_pos' in entry_format:
                    b = datetime_to_timestamp(line[b_pos])
                    e = datetime_to_timestamp(line[e_pos])
                    link_duration = e - b
                    t = b
            else:
                if 't_pos' in entry_format:
                    t = float(line[t_pos].replace(',', ''))
                elif 'b_pos' in entry_format:
                    b = float(line[b_pos].replace(',', ''))
                    e = float(line[e_pos].replace(',', ''))
                    link_duration = e - b
                    t = b

            min_t, max_t = min(min_t, t), max(max_t, t + link_duration)

            if kwargs['link_duration']:
                link_duration = kwargs['link_duration']
            elif 'link_duration_pos' in entry_format:
                link_duration = float(line[link_duration_pos].replace(',', ''))

            if kwargs['is_directed']:
                l = (u, v)
            else:
                if (v, u) in E:
                    l = (v, u)
                else:
                    l = (u, v)
            if u == v:
                # SELF LOOP : we ignore it
                continue

            if l in E and E[l][-1] >= t:
                E[l][-1] = max(E[l][-1], t + link_duration)
            else:
                E[l] += [t, t + link_duration]
            if kwargs['is_link_stream'] is False:
                if u in W and W[u][-1] >= t:
                    W[u][-1] = max(W[u][-1], t + link_duration)
                else:
                    W[u] += [t, t + link_duration]

                if v in W and W[v][-1] >= t:
                    W[v][-1] = max(W[v][-1], t + link_duration)
                else:
                    W[v] += [t, t + link_duration]

    S = sg.StreamGraph(times=[min_t, max_t],
                       nodes=list(W.keys()),
                       links=list(E.keys()),
                       node_presence=[W[k] for k in W.keys()],
                       link_presence=[E[k] for k in E.keys()],
                       node_to_label=id_to_label)
    return S


def parse_link_stream(input_file):
    """
    Parse link stream format:
    alpha t0
    omega t1
    b e u v
    .
    .
    .
    b e v w

    :param input_file:
    :return:
    """
    E = defaultdict(list)
    W = defaultdict(list)
    cnt_rows = 0

    nodes_to_label = {}
    label_to_id = defaultdict(lambda: len(label_to_id))
    with open(input_file, 'r') as ipt:
        size = sum(1 for _ in ipt)

    with open(input_file, 'r') as input_file:
        for line in tqdm(input_file, total=size):
            cnt_rows += 1
            # if cnt_rows % 100000 == 0:
            #     print((cnt_rows / size) * 100, "% loaded")
            l = line.strip().split()
            if len(l) == 2:
                assert l[0] in ["alpha", "omega"]
                if l[0] == "alpha":
                    alpha = float(l[1])
                else:
                    omega = float(l[1])
            else:
                assert (len(l) == 4)
                b, e, u_label, v_label = l

                u = label_to_id[u_label]
                v = label_to_id[v_label]
                nodes_to_label[u] = u_label
                nodes_to_label[v] = v_label

                b = float(b)
                e = float(e)

                l = (u, v)
                if l in E:
                    l = (v, u)
                if l in E and E[l][-1] >= b:
                    E[l][-1] = max(E[l][-1], e)
                else:
                    E[l] += [b, e]
                if u not in W:
                    W[u] = [alpha, omega]
                if v not in W:
                    W[v] = [alpha, omega]

    S = sg.StreamGraph(times=[alpha, omega],
                       nodes=list(W.keys()),
                       links=list(E.keys()),
                       node_presence=[W[k] for k in W.keys()],
                       link_presence=[E[k] for k in E.keys()],
                       node_to_label=nodes_to_label,
                       node_to_id={i: i for i in W.keys()})
    return S


def parse_pcap(file_input, entry_format, **options):
    pcap_to_csv(file_input, "tmp.csv")
    S = parse_csv("tmp.csv", entry_format, **options)
    os.remove("tmp.csv")
    return S


def parser(input_file, input_format, entry_format, output_file=None, simplify_presence=False, output_format='sg',
           **kwargs):
    """
    Straph's tunable parser. Compatible with several data formats: CSV, TSV, JSon and PCAP.

    :param simplify_presence:
    :param input_file: Input FILE (name only)
    :param input_format: Format d'entrée acceptés : JSON, CSV, PCAP
    :param entry_format: Format of each line to be readed (t,u,v) = (line[x],line[y],line[w])
    :param output_file: Output FILE (name only)
    :param output_format: Format de sortie : SG,SGF,json
    :return:
    """
    with open(input_file) as ipt:
        options = {'delimiter': ',',
                   'is_link_stream': False,
                   'is_directed': False,
                   'nrows': sum(1 for _ in ipt),
                   'link_duration': False,
                   'order_sgf': False,
                   'ignore_header': True,
                   'nodes_to_label': False,
                   'time_is_datetime': False,
                   'delta': None,
                   }
    options.update(kwargs)
    if ('t_pos' in entry_format or 'link_duration_pos' in entry_format) and \
            ('b_pos' in entry_format or 'e_pos' in entry_format):
        raise TypeError('Invalide entry format :' + str(entry_format) + ' should be of type {t_pos,u_pos,v_pos} or'
                                                                        ' {t_pos,link_duration_pos,u_pos,v_pos} or'
                                                                        '{b_pos,e_pos,u_pos,v_pos} !')
    if options['link_duration'] and ('b_pos' in entry_format or 'e_pos' in entry_format):
        raise TypeError('link_duration is incompatible with entry format : {b_pos,e_pos,u_pos,v_pos} !')

    if options['link_duration'] and ('link_duration_pos' in entry_format):
        raise TypeError('link_duration is incompatible with entry format : {t_pos,link_duration_pos,u_pos,v_pos} !')

    if input_format == 'csv':
        S = parse_csv(input_file, entry_format, **options)
    elif input_format == 'json':
        S = parse_json(input_file, entry_format, **options)
    elif input_format == 'pcap':
        S = parse_pcap(input_file, entry_format, **options)
    elif input_format == 'net':
        raise ValueError("File format 'net' not yet supported.")
        # S = parse_net(input_format, entry_format, **options)
    else:
        raise TypeError('Format not supported')

    if simplify_presence is True:
        S.node_presence = [[np[0], np[-1]] for np in
                           S.node_presence]  # Set nodes to be present from ther 1st intercations to their last

    if output_file is not None:
        if isinstance(output_format, str):
            output_format = [output_format]

        for of in output_format:
            if of == 'sg':
                S.write_to_sg(output_file)
            elif of == 'json':
                S.write_to_json(output_file)

    return S


def sort_csv(input_file, entry_format, output=None, **kwargs):
    with open(input_file) as ipt:
        options = {'delimiter': ',',
                   'is_link_stream': False,
                   'is_directed': False,
                   'nrows': sum(1 for _ in ipt),
                   'link_duration': False,
                   'order_sgf': False,
                   'ignore_header': True,
                   'nodes_to_label': False,
                   'time_is_datetime': False,
                   'delta': None,
                   }
    options.update(kwargs)

    list_lines = []
    with open(input_file, 'r') as input:
        reader = csv.reader(input, delimiter=options['delimiter'])
        if options['ignore_header']:
            next(reader, None)
        for line in tqdm(reader, desc='Reading CSV before sorting', total=options['nrows']):
            list_lines.append(line)

    if len(entry_format) == 3:
        (t_pos, u_pos, v_pos) = entry_format['t_pos'], entry_format['u_pos'], entry_format['v_pos']
        if options['time_is_datetime']:
            list_lines = sorted(list_lines, key=lambda x: datetime_to_timestamp(x[t_pos]))
        else:
            list_lines = sorted(list_lines, key=lambda x: float(x[t_pos]))
    elif len(entry_format) == 4 and 'link_duration_pos' in entry_format:
        (t_pos, link_duration_pos, u_pos, v_pos) = entry_format['t_pos'], entry_format['link_duration_pos'], \
                                                   entry_format['u_pos'], entry_format['v_pos']
        if options['time_is_datetime']:
            list_lines = sorted(list_lines, key=lambda x: datetime_to_timestamp(x[t_pos]))
        else:
            list_lines = sorted(list_lines, key=lambda x: float(x[t_pos]))
    elif len(entry_format) == 4 and 'b_pos' in entry_format:
        (b_pos, e_pos, u_pos, v_pos) = entry_format['b_pos'], entry_format['e_pos'], \
                                       entry_format['u_pos'], entry_format['v_pos']

        if options['time_is_datetime']:
            list_lines = sorted(list_lines, key=lambda x: datetime_to_timestamp(x[b_pos]))
        else:
            list_lines = sorted(list_lines, key=lambda x: float(x[b_pos]))
    if output is None:
        output = input_file
    with open(output, 'w', newline='') as output:
        writer = csv.writer(output, delimiter=options['delimiter'])
        for line in tqdm(list_lines, desc='Writing CSV'):
            writer.writerow(line)
