import csv, json
import dpkt, socket
import time
import dateutil.parser as du
import msgpack
import math

from collections import defaultdict
from sortedcollections import SortedSet

from straph import stream as sg

import matplotlib.pyplot as plt

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
    """Convert inet object to a string

        Args:
            inet (inet struct): inet network address
        Returns:
            str: Printable/readable IP address
    """
    # First try ipv4 and then ipv6
    try:
        return socket.inet_ntop(socket.AF_INET, inet)
    except ValueError:
        return socket.inet_ntop(socket.AF_INET6, inet)


def datetime_to_timestamp(s):
    return du.parse(s).timestamp()


def pcap_to_csv(file_input, destination,
                path_dict, protocol=None):
    """
    Transform a pcap file to a csv
    :param file_input:
    :param destination:
    :param path_dict:
    :param protocol:
    :return:
    """
    counter = 0
    dict_ip = defaultdict(lambda: len(dict_ip))
    dict_label = {}
    start = time.time()
    if protocol:
        protocol = [key for key, value in __nb_2_protocol__.items() if value == protocol][0]
    print("protocol :", protocol)
    with open(destination, 'w') as output, open(file_input, 'rb') as input:
        writer = csv.writer(output, delimiter=';')
        writer.writerow(["time", "src", "dst", "protocol", "len", "src_port", "dst_port"])
        for ts, pkt in dpkt.pcap.Reader(input):
            eth = dpkt.ethernet.Ethernet(pkt)
            if counter == 150000000000000000000000000000:
                print("Time end dkpt :", time.time() - start)
                break
            if counter == 0:
                t0 = ts
            counter += 1
            if isinstance(eth.data, bytes):
                continue
            ip = eth.data
            try:
                ip_src = inet_to_str(ip.src).encode()
            except:
                continue
            ip_dst = inet_to_str(ip.dst).encode()

            id_src = dict_ip[ip_src]
            id_dst = dict_ip[ip_dst]
            dict_label[id_src] = ip_src
            dict_label[id_dst] = ip_dst

            if counter % 1000000 == 0:
                print("Counter dkpt:", counter, "time dkpt:", time.time() - start)
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
    with open(path_dict + "_dict_ip_2_nodes_label.mspk", 'wb') as output:
        msgpack.dump(dict_ip, output)
    with open(path_dict + "_dict_node_label_2_ip.mspk", 'wb') as output:
        msgpack.dump(dict_label, output)


def pcap_to_sgf(file_input, destination, time_bins,
                path_dict, protocol=None, threshold_collect=5000000, with_dict=True):
    """
    Transform a pcap file to a Stream Graph File
    :return:
    """
    cnt = 0
    dict_ip = defaultdict(lambda: len(dict_ip))
    dict_label = {}
    links_to_add = []
    dict_link_to_presence = {}
    packer = msgpack.Packer()
    start = time.time()
    if protocol:
        protocol = [key for key, value in __nb_2_protocol__.items() if value == protocol][0]
    print("protocol :", protocol)
    with open(destination + '_links.sgf', 'wb') as output, open(file_input, 'rb') as input:
        for ts, pkt in dpkt.pcap.Reader(input):
            eth = dpkt.ethernet.Ethernet(pkt)
            if cnt == 150000000000000000000000000000:
                print("Time end dkpt :", time.time() - start)
                break
            if cnt == 0:
                t0 = ts
            if isinstance(eth.data, bytes):
                continue
            ip = eth.data
            try:
                ip_src = inet_to_str(ip.src).encode()
            except:
                continue
            ip_dst = inet_to_str(ip.dst).encode()
            if with_dict:
                id_src = dict_ip[ip_src]
                id_dst = dict_ip[ip_dst]
                dict_label[id_src] = ip_src
                dict_label[id_dst] = ip_dst
            # We ignore 'ICMP' protocols, ICMP scan useless
            if ip.p == 1:
                continue
            if protocol:
                if ip.p != protocol:
                    continue
            cnt += 1
            current_time = round(ts - t0, 6)
            if cnt == threshold_collect:  # Every time we depass 'threshold_collect' we collect the links
                print("Counter dkpt:", cnt, "time dkpt:", time.time() - start)
                cnt = 0
                to_del = []
                for l, (t0, t1) in dict_link_to_presence.items():
                    if t1 < current_time - time_bins:
                        links_to_add.append((l, t0, t1))
                        to_del.append(l)
                for l in to_del:
                    del dict_link_to_presence[l]

            l = (id_src, id_dst)
            if l in dict_link_to_presence:
                if dict_link_to_presence[l][1] < current_time - time_bins:
                    # Previous Interaction, ended
                    # Store terminated link
                    links_to_add.append((l, dict_link_to_presence[l][0], dict_link_to_presence[l][1]))
                    dict_link_to_presence[l] = [current_time, current_time + time_bins]
                else:
                    # Extend previous interaction
                    dict_link_to_presence[l][1] = current_time + time_bins
            else:
                dict_link_to_presence[l] = [current_time, current_time + time_bins]
        for l, (t0, t1) in dict_link_to_presence.items():
            links_to_add.append((l, t0, t1))
        del dict_link_to_presence
        write_links_to_sgf(output, packer, links_to_add)

    if with_dict:
        with open(path_dict + "_dict_ip_2_nodes_label.mspk", 'wb') as output:
            msgpack.dump(dict_ip, output)
        with open(path_dict + "_dict_node_label_2_ip.mspk", 'wb') as output:
            msgpack.dump(dict_label, output)


def write_links_to_sgf(destination, packer, links):
    links = list(links)
    links = sorted(links, key=lambda x: x[0])
    with open(destination, 'wb') as output_file:
        for l in links:
            if l[0] > l[1]:
                raise AssertionError("T0 > T1")
            output_file.write(packer.pack(l))


def csv_to_sgf(file_input, destination,
               time_bins, threshold_collect=5000000, directed_links=False):
    """
    Transform a csv file into a Stream Graph File (Links only).
    Streaming File ordered by time of arrival
    :param file_input:
    :param destination:
    :param time_bins:
    :param threshold_collect:
    :param directed_links:
    :return:
    """
    links_to_add = set()
    link_2_presence_time = {}
    ignore_header = True
    cnt = 0
    t_begin = time.time()
    with open(file_input, 'r') as input:
        reader = csv.reader(input, delimiter=';')
        for line in reader:
            if ignore_header:
                ignore_header = False
                continue
            cnt += 1
            if cnt % threshold_collect == 0:
                print(threshold_collect, " lines browsed !")
            current_time = float(line[0])
            l = (int(line[1]), int(line[2]))
            if not directed_links:
                if (l[1], l[0]) in link_2_presence_time:  # We consider undirected links
                    l = (l[1], l[0])
            if l in link_2_presence_time:
                if link_2_presence_time[l][1] < current_time:
                    links_to_add.add((link_2_presence_time[l][0], link_2_presence_time[l][1], l[0], l[1],))
                    link_2_presence_time[l] = [current_time, current_time + time_bins]
                else:
                    # Extend previous interaction
                    link_2_presence_time[l][1] = max(current_time + time_bins, link_2_presence_time[l][1])
            else:
                link_2_presence_time[l] = [current_time, current_time + time_bins]
    # LAST PASS :
    for l, (t0, t1) in link_2_presence_time.items():
        links_to_add.add((t0, t1, l[0], l[1]))
    del link_2_presence_time
    packer = msgpack.Packer(use_bin_type=True, )
    write_links_to_sgf(destination, packer, links_to_add)
    print("Nb lines :", cnt)
    print("Time extraction :", time.time() - t_begin)
    return


def order_sgf(path_input, path_output):
    packer = msgpack.Packer(use_bin_type=True)
    links = []
    with open(path_input, 'rb') as input, open(path_output, 'wb') as output:
        unpacker = msgpack.Unpacker(input, use_list=False)
        for i in unpacker:
            t0, t1, u, v = i
            links.append((1, t0, t1, u, v))  # code each link, 1 for a beginning, -1 for an ending
            links.append((-1, t1, u, v))
        links = sorted(links, key=lambda x: (x[1], -x[0]))
        for l in links:
            output.write(packer.pack(l))
    return


def parse_net(input_file, output_file_nodes, output_file_links, link_duration=1):
    """
    A Stream Graph reader for dataset issued by Pajek
    Format of interactions : .net
    :param input_file:
    :param output_file_nodes:
    :param output_file_links:
    :param link_duration:
    :param delimiter:
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
            if type_node == True:
                continue
            if type_node == False:
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



def BurningBush(input_file, output_file_nodes, output_file_links, delimiter=';', null_duration_allowed=True,
                nrows=False):
    """
    A Stream Graph reader for burningBUSH:
    :param input_file:
    :param output_file_nodes:
    :param output_file_links:
    :param link_duration:
    :param delimiter:
    :return:
    """
    E = defaultdict(list)
    W = defaultdict(list)
    cnt_row = 0
    with open(input_file, 'r') as input_file:
        reader = csv.reader(input_file, delimiter=delimiter)
        next(reader, None)  # skip header :)
        for line in reader:
            l = (line[7], line[3])
            current_time = du.parser.parse(line[16]).timestamp()
            link_duration = float(line[9])
            if null_duration_allowed is False and link_duration == 0:
                continue
            if nrows and cnt_row >= nrows:
                break
            cnt_row += 1

            if (l[1], l[0]) in E:
                l = (l[1], l[0])
            u, v = l
            if u == v:
                # SELF LOOP : we ignore it
                continue
            if l in E and E[l][-1] >= current_time - link_duration:
                # print("Extend Link Presence")
                E[l][-1] = max(E[l][-1], current_time)
            else:
                E[l] += [current_time - link_duration, current_time]

            if u in W and W[u][-1] >= current_time - link_duration:
                # print("Extend Node Presence")
                W[u][-1] = max(W[u][-1], current_time)
            else:
                W[u] += [current_time - link_duration, current_time]

            if v in W and W[v][-1] >= current_time - link_duration:
                # print("Extend Node Presence")
                W[v][-1] = max(W[v][-1], current_time)
            else:
                W[v] += [current_time - link_duration, current_time]
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



def parse_csv(input_file, entry_format, **kwargs):
    """
    Reader for .csv files

    :param input_file:
    :param entry_format:
    :param kwargs:
    :return:
    """

    # Convert entry format
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
        elif 'link_duration_pos' not in entry_format:
            link_duration = 0
            print("[WARNING] No link_duration provided, links durations are set to 0.")

        for line in reader:
            cnt_rows += 1
            if cnt_rows % 1000000 == 0:
                print((cnt_rows / kwargs['size']) * 100, "% loaded")

            if kwargs['nrows'] and cnt_rows > kwargs['nrows']:
                break

            if kwargs['nodes_to_label']:
                # Convert Label to int
                u_label = line[u_pos]
                v_label = line[v_pos]
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
        print("Approximate events with delta = ", delta)
        W, E = approximate_events(W, E, delta)

    S = sg.stream_graph(times=[min_t, max_t],
                        nodes=list(W.keys()),
                        links=list(E.keys()),
                        node_presence=[W[k] for k in W.keys()],
                        link_presence=[E[k] for k in E.keys()],
                        node_to_label=nodes_to_label,
                        node_to_id={i: i for i in W.keys()})
    return S


def approximate_events(W, E, delta):
    event_times = sorted(set([t for np in W.values() for t in np] + [t for lp in E.values() for t in lp]))
    t_old = event_times[0]
    discretized_event_times = SortedSet()
    discretized_event_times.add(t_old)
    for t in event_times[1:]:
        if t - t_old >= delta:
            discretized_event_times.add(t)
            t_old = t

    new_W = {}
    for n, np in W.items():
        new_W[n] = []
        for t0, t1 in zip(np[::2], np[1::2]):
            #  assert t1-t0 >= delta
            if t0 not in discretized_event_times:
                # Catch time after t0 in discretized event times:
                t0 = discretized_event_times[discretized_event_times.bisect(t0)]
            if t1 not in discretized_event_times:
                # Catch time before t1 in discretize event times:
                t1 = discretized_event_times[discretized_event_times.bisect(t1)-1]

            new_W[n] += [t0,t1]

    new_E = {}
    for l, lp in E.items():
        new_E[l] = []
        for t0, t1 in zip(lp[::2], lp[1::2]):
            # assert t1-t0 >= delta
            if t0 not in discretized_event_times:
                # Catch time after t0 in discretized event times:
                t0 = discretized_event_times[discretized_event_times.bisect(t0)]
            if t1 not in discretized_event_times:
                # Catch time before t1 in discretize event times:
                t1 = discretized_event_times[discretized_event_times.bisect(t1)-1]

            new_E[l] += [t0,t1]
    return new_W, new_E


def parse_json(input_file, entry_format, **kwargs):
    """
    A Stream Graph reader for JSON dataset:
    [u,v,t]

    :param input_file:
    :param entry_format:
    :param kwargs:
    :return:
    """

    # Convert entry format
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
        for line in reader:
            cnt_rows += 1

            if cnt_rows % 100000 == 0:
                print((cnt_rows / kwargs['size']) * 100, "% loaded")

            if kwargs['nrows'] and cnt_rows > kwargs['nrows']:
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

    S = sg.stream_graph(times=[min_t, max_t],
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
        for line in input_file:
            cnt_rows += 1
            if cnt_rows % 100000 == 0:
                print((cnt_rows / size) * 100, "% loaded")
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

    S = sg.stream_graph(times=[alpha, omega],
                        nodes=list(W.keys()),
                        links=list(E.keys()),
                        node_presence=[W[k] for k in W.keys()],
                        link_presence=[E[k] for k in E.keys()],
                        node_to_label=nodes_to_label,
                        node_to_id={i: i for i in W.keys()})
    return S


def parse_pcap(input_file, entry_format, **kwargs):
    return


def parser(input_file, input_format, entry_format, output_file=None, output_format='sg', **kwargs):
    """
    TODO : Weighted nodes, weighted links (then modify 'entry_format').
    :param input_file: Input FILE (name only)
    :param input_format: Format d'entrée acceptés : JSON, CSV, PCAP
    :param entry_format: Format of each line to be readed (t,u,v) = (line[x],line[y],line[w])
    :param output_file: Output FILE (name only)
    :param output_format: Format de sortie : SG,SGF,json
    :param is_link_stream: If 'True' nodes are always present
    :param is_directed: IF 'False' we checked if the inverse link is already present and avoid any duplicate.
    :param ndoes_to_label: If 'True' : creation of 'nodes_to_label' attribute
    :param delimiter: Delimiter in CSV format
    :param nrows: Number of line to read
    :param link_duration: Duration of link, if not specified
    :param order_sgf: If 'True' return an ordered SGF in increasing order
    :param delta: Force links that begin/end in the same window *[t-delta/2,t+delta/2]* to begin/end at *t*
    :return:
    """
    with open(input_file) as ipt:
        options = {'delimiter': ',',
                   'is_link_stream': False,
                   'is_directed': False,
                   'nrows': None,
                   'link_duration': False,
                   'order_sgf': False,
                   'ignore_header': True,
                   'nodes_to_label': False,
                   'time_is_datetime': False,
                   'delta': None,
                   'size': (kwargs['nrows'] if 'nrows' in kwargs else sum(1 for _ in ipt))}
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
        S = parse_net(input_format, entry_format, **options)
    else:
        raise TypeError('Format not supported')

    if output_file is not None:
        if isinstance(output_format, str):
            output_format = [output_format]

        for of in output_format:
            if of == 'sg':
                S.write_to_sg(output_file)
            elif of == 'sgf':
                S.write_to_sgf(output_file)
            elif of == 'json':
                S.write_to_json(output_file)

    return S


def sort_csv(input_file, entry_format, **kwargs):
    list_lines = []
    with open(input_file, 'r') as input:
        reader = csv.reader(input, delimiter=kwargs['delimiter'])
        if kwargs['ignore_header']:
            next(reader, None)
        for line in reader:
            list_lines.append(line)

    if len(entry_format) == 3:
        (t_pos, u_pos, v_pos) = entry_format['t_pos'], entry_format['u_pos'], entry_format['v_pos']
        list_lines = sorted(list_lines, key=lambda x: float(x[t_pos]))
    elif len(entry_format) == 4 and 'link_duration_pos' in entry_format:
        (t_pos, link_duration_pos, u_pos, v_pos) = entry_format['t_pos'], entry_format['link_duration_pos'], \
                                                   entry_format['u_pos'], entry_format['v_pos']
        list_lines = sorted(list_lines, key=lambda x: float(x[t_pos]))
    elif len(entry_format) == 4 and 'b_pos' in entry_format:
        (b_pos, e_pos, u_pos, v_pos) = entry_format['b_pos'], entry_format['e_pos'], \
                                       entry_format['u_pos'], entry_format['v_pos']
        list_lines = sorted(list_lines, key=lambda x: float(x[b_pos]))
    with open(input_file, 'w') as output:
        writer = csv.writer(output, delimiter=kwargs['delimiter'])
        for line in list_lines:
            writer.writerow(line)