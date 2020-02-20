



######################################################
#          Paths Algorithms                          #
######################################################

# Foremost paths
def FoP(S, source, destination=None, start_time=None):
    """
    If *destination* is specified: Return the time to reach *destination* from *source* in *S*.
    Otherwise: Return the time to reach every reachable node in *S* from *source*.

    :param S: Stream graph
    :param source: temporal source node
    :param destination: optional
    :param start_time : optional, if not specified assume that the start time is source[0]
    :return:

    """
    L_functions = {'initialisation': FoP_initialisation,
                   'update': FoP_update,
                   }
    ttr = None
    if type(source) == int:  #  source node
        temporal_source_list = get_temporal_sources(S, source)
        for s in temporal_source_list:
            if start_time is not None and s[0] <= start_time <= s[1]:
                R = L_Algorithm(S, s, L_functions,
                                start_time=start_time)
                if R is not None:
                    ttr = FoP_postprocess(R, s, destination,
                                          start_time=start_time,
                                          ttr=ttr)
                break
            else:
                R = L_Algorithm(S, s, L_functions,
                                start_time=start_time)
                if R is not None:
                    ttr = FoP_postprocess(R, s, destination,
                                          start_time=start_time,
                                          ttr=ttr)

    else:  #  temporal source node
        R = L_Algorithm(S, source, L_functions,
                        start_time=start_time)
        if R is not None:
            ttr = FoP_postprocess(R, source, destination,
                                  start_time=start_time)

    return ttr


# Fastest Paths
def FP(S, source, destination=None):
    L_functions = {'initialisation': FP_initialisation,
                   'update': FP_update,
                   }
    latencies = None
    if type(source) == int:  #  source node
        temporal_source_list = get_temporal_sources(S, source)
        for s in temporal_source_list:
            R = L_Algorithm(S, s, L_functions)
            if R is not None:
                latencies = FP_postprocess(R, destination,
                                           latencies=latencies)
    else:  #  temporal source node
        R = L_Algorithm(S, source, L_functions)
        if R is not None:
            latencies = FP_postprocess(R, destination,
                                       latencies=latencies)
    if destination:
        print("R :", [R[i] for i in R if i[2] == destination])
        print("lat :", latencies)
    return latencies

# Shortest Foremost Paths
def SFoP(S, source, destination=None, start_time=None):
    L_functions = {'initialisation': SFoP_initialisation,
                   'update': SFoP_update,
                   }
    ttr, lengths = None, None
    if type(source) == int:  #  source node
        temporal_source_list = get_temporal_sources(S, source)
        for s in temporal_source_list:
            if start_time is not None and s[0] <= start_time <= s[1]:
                R = L_Algorithm(S, s, L_functions,
                                start_time=start_time)
                if R is not None:
                    ttr, lengths = SFoP_postprocess(R, s, destination,
                                                    start_time=start_time,
                                                    ttr=ttr,
                                                    lengths=lengths)
                break
            else:
                R = L_Algorithm(S, s, L_functions,
                                start_time=start_time)
                if R is not None:
                    ttr, lengths = SFoP_postprocess(R, s, destination,
                                                    start_time=start_time,
                                                    ttr=ttr,
                                                    lengths=lengths)

    else:  #  temporal source node
        R = L_Algorithm(S, source, L_functions,
                        start_time=start_time)
        if R is not None:
            ttr, lengths = SFoP_postprocess(R, source, destination,
                                            start_time=start_time)

    return ttr, lengths

# Shortest Fastest Paths
def SFP(S, source, destination=None):
    L_functions = {'initialisation': SFP_initialisation,
                   'update': SFP_update,
                   }
    latencies, lengths = None, None
    if type(source) == int:  #  source node
        temporal_source_list = get_temporal_sources(S, source)
        for s in temporal_source_list:
            R = L_Algorithm(S, s, L_functions, multi_criteria=True)
            if R is not None:
                latencies, lengths = SFP_postprocess(R, destination,
                                                     latencies=latencies,
                                                     lengths=lengths)
    else:  #  temporal source node
        R = L_Algorithm(S, source, L_functions, multi_criteria=True)
        if R is not None:
            latencies, lengths = SFP_postprocess(R, destination,
                                                 latencies=latencies,
                                                 lengths=lengths)
    return latencies, lengths



def FSP(S, source, destination=None):
    L_functions = {'initialisation': FSP_initialisation,
                   'update': FSP_update,
                   }
    distances, durations = None, None
    if type(source) == int:  #  source node
        temporal_source_list = get_temporal_sources(S, source)
        for s in temporal_source_list:
            R = L_Algorithm(S, s, L_functions)
            if R is not None:
                distances, durations = FSP_postprocess(R, destination,
                                                       distances=distances,
                                                       durations=durations)
    else:  #  temporal source node
        R = L_Algorithm(S, source, L_functions)
        if R is not None:
            distances, durations = FSP_postprocess(R, destination,
                                                   distances=distances,
                                                   durations=durations)
    return distances, durations


######################################
#           Foremost Paths (FoP)     #
######################################

def FoP_initialisation(S, source):
    """
    'L' initialisation for FoP

    :param source:
    :return:
    """
    L = {}
    if type(source) == int:
        L[source] = S.node_presence[source][0]  # (a_u)
    else:
        L[source[1]] = source[0]  # (a_u)
    return L, L


def FoP_update(L, R, t0, t1, u, v):
    """
    We update 'L' according to the link (t0,t1,u,v) and FoP properties.

    :param L:
    :param t0:
    :param t1:
    :param u:
    :param v:
    :return:
    """
    new_arrival = t0
    # If v is already in L, we've reached v at t <= t0.
    if v not in L:
        L[v] = new_arrival
        return True
    return False


def FoP_postprocess(R, source, destination=None, start_time=None, ttr=None):
    """
    Post process elements in 'L' to obtain time to reach 'ttr'

    :param L:
    :param source:
    :param destination: optional
    :return:
    """
    if start_time is None:
        start_time = source[0]

    if destination is None:  # Source - Target
        if ttr is None:
            ttr = {}
        for k, v in R.items():
            n = k[2]
            if n in ttr:
                ttr[n] = max(min(ttr[n], v - start_time), 0)
            else:
                ttr[n] = max(v - start_time, 0)

    else:  #  Single source

        if ttr is None:
            ttr = math.inf
        if type(destination) == int:  # destination node
            for k, v in R.items():
                if k[2] == destination:
                    ttr = max(min(ttr, v - start_time), 0)
        else:  # temporal destination node
            for k, v in R.items():
                if k[2] == destination[2] and k[0] <= destination[0] <= destination[1] <= k[1]:
                    ttr = max(min(ttr, v - start_time), 0)

    return ttr



def FoP_postprocess_pw(W):
    ttr_pw = {}
    for n in W:
        ttr_pw[n] = FoP_postprocess(W[n], n)
    return ttr_pw


def FoP_pw(S):
    """
    Return the time to reach every reachable node in *S* from every other node.
    (infinity if not reachable)
    :param S: Stream graph
    :return:
    """
    L_functions = {'initialisation': FoP_initialisation,
                   'update': FoP_update,
                   }
    W = F_Algorithm(S, L_functions)
    if W is not None:
        times_to_reach_pw = FoP_postprocess_pw(W)
        return times_to_reach_pw


######################################
#           Shortest Foremost Paths  #
######################################

def SFoP_initialisation(S, source):
    L, R = {}, {}

    if type(source) == int:

        L[source] = (0, S.node_presence[source][0])  # (d_v,a_v)
        R[source] = (0, S.node_presence[source][0])
    else:
        x, t = source
        L[x] = (0, t)  # (d_v,a_v)
        R[x] = (0, t)
    return L, R


def SFoP_update(L, R, t0, t1, u, v):
    du, au = L[u]
    new_arrival = t0
    new_distance = du + 1
    if v in L:
        current_best_distance, current_best_arrival = R[v]
        if new_arrival == current_best_arrival and new_distance < current_best_distance:
            R[v] = (new_distance, new_arrival)

        dv, av = L[v]
        if new_distance < dv:
            L[v] = (new_distance, new_arrival)
            return True
        elif new_distance == dv:
            return True
        else:
            return False
    else:
        L[v] = (new_distance, new_arrival)
        R[v] = (new_distance, new_arrival)
    return True


def SFoP_postprocess(R, source, destination=None, start_time=None,
                     ttr=None, lengths=None):
    '''
    Post Process elements in 'L' to obtain times to reach and lengths of shortest foremost path

    :param L:
    :param source:
    :param destination:
    :param start_time:
    :param ttr:
    :param lengths:
    :return:
    '''
    if start_time is None:
        start_time = source[0]

    if destination is None:
        if ttr is None:
            ttr = {}
            lengths = {}
        for k, v in R.items():
            n = k[2]
            arrival_time = v[1]
            l = v[0]
            if n in lengths:
                if arrival_time - start_time <= ttr[n]:
                    lengths[n] = min(l, lengths[n])
                    ttr[n] = min(ttr[n], arrival_time - start_time)
            else:
                lengths[n] = l
                ttr[n] = arrival_time - start_time

    else:  #  Single source mode
        if ttr is None:
            ttr = math.inf
            lengths = math.inf
        if type(destination) == int:  # destination node
            for k, v in R.items():
                if k[2] == destination:
                    arrival_time = v[1]
                    l = v[0]
                    if arrival_time - start_time <= ttr:
                        ttr = max(min(ttr, arrival_time - start_time), 0)
                        lengths = min(l, lengths)
        else:  # temporal destination node
            for k, v in R.items():
                if k[2] == destination[2] and k[0] <= destination[0] <= destination[1] <= k[1]:
                    arrival_time = v[1]
                    l = v[0]
                    if arrival_time - start_time <= ttr:
                        ttr = max(min(ttr, arrival_time - start_time), 0)
                        lengths = min(l, lengths)
    return ttr, lengths





def SFoP_W_postprocess(W):
    ttr_pw = {}
    length_pw = {}
    for n in W:
        ttr_pw[n], length_pw[n] = SFoP_postprocess(W[n], n)
    return ttr_pw, length_pw


def SFoP_pw(S):
    L_functions = {'initialisation': SFoP_initialisation,
                   'update': SFoP_update,
                   }
    W = F_Algorithm(S, L_functions)
    if W is not None:
        ttr_pw, lengths_pw = SFoP_W_postprocess(W)
        return ttr_pw, lengths_pw


######################################
#           Shortest Paths           #
######################################







def SP_W_postprocess(W):
    distances_pw = {}
    for n in W:
        distances_pw[n] = SP_postprocess(W[n])
    return distances_pw


def SP_pw(S):
    L_functions = {'initialisation': SP_initialisation,
                   'update': SP_update,
                   }
    W = F_Algorithm(S, L_functions)
    if W is not None:
        distances_pw = SP_W_postprocess(W)
        return distances_pw


######################################
#           Fastest Shortest Paths   #
######################################

def FSP_initialisation(source):
    L, R = {}, {}
    L[source] = (0, -source[1], source[0])  # (d_u,s_u,a_u)
    R[source] = (0, source[1], source[0])
    return L, R


def FSP_update(L, R, t0, t1, u, v):
    du, su, au = L[u]
    su = -su
    new_arrival = t0
    new_distance = du + 1
    new_start = min(t1, su)
    new_duration = max(new_arrival - new_start, 0)
    if v in L:

        current_best_distance = R[v][0]
        current_best_duration = max(R[v][2] - R[v][1], 0)
        if new_distance < current_best_distance or (new_distance == current_best_distance and
                                                    new_duration < current_best_duration):
            R[v] = (new_distance, new_start, new_arrival)

        dv, sv, av = L[v]
        sv = -sv
        if new_distance < dv or (new_distance == dv and new_start > sv):
            L[v] = (new_distance, -new_start, new_arrival)
        else:
            return False
    else:
        L[v] = (new_distance, -new_start, new_arrival)
        R[v] = (new_distance, new_start, new_arrival)
    return True


def FSP_postprocess(R, destination=None,
                    distances=None,
                    durations=None):
    if destination is None:
        if distances is None:
            distances = {}
            durations = {}
        for k, v in R.items():
            n = k[2]
            length = v[0]
            duration = v[2] - v[1]
            if duration < 0:
                duration = 0
            if n in distances:
                if length <= distances[n]:
                    distances[n] = min(distances[n], length)
                    durations[n] = min(durations[n], duration)
            else:
                distances[n] = length
                durations[n] = duration
    else:
        if distances is None:
            distances = math.inf
            durations = math.inf
        if type(destination) == int:  # destination node
            for k, v in R.items():
                if k[2] == destination:
                    length = v[0]
                    duration = v[2] - v[1]
                    if duration < 0:
                        duration = 0
                    if length <= distances:
                        distances = min(distances, length)
                        durations = min(durations, duration)
        else:  # temporal destination node
            for k, v in R.items():
                if k[2] == destination[2] and k[0] <= destination[0] <= destination[1] <= k[1]:
                    length = v[0]
                    duration = v[2] - v[1]
                    if duration < 0:
                        duration = 0
                    if length <= distances:
                        distances = min(distances, length)
                        durations = min(durations, duration)
    return distances, durations





def FSP_W_postprocess(W):
    distances_pw = {}
    duration_pw = {}
    for n in W:
        distances_pw[n], duration_pw[n] = FSP_postprocess(W[n])
    return distances_pw, duration_pw


def FSP_pw(S):
    L_functions = {'initialisation': FSP_initialisation,
                   'update': FSP_update,
                   }
    W = F_Algorithm(S, L_functions)
    if W is not None:
        distances_pw, duration_pw = FSP_W_postprocess(W)
        return distances_pw, duration_pw


######################################
#           Fastest Paths            #
######################################

#  IMPORTANT : On indexe par L par -s_u au lieu de s_u pour gérer les batch arrival !!

def FP_initialisation(source):
    L, R = {}, {}
    L[source] = (-source[1], source[0])  # (s_u,a_u)
    R[source] = (source[1], source[0])
    return L, R


def FP_update(L, R, t0, t1, u, v):
    '''
    Lu is sorted by the value of au and we analyze temporal links in temporal order.
    :param L:
    :param t0:
    :param t1:
    :param u:
    :param v:
    :return:
    '''

    su, au = L[u]
    su = -su  #  Trick for the sorted queue in bfs update
    new_arrival = t0
    new_start = min(t1, su)  # Remember that the aim is to leaving the source node the latest possible (thus s_u
    #  can be > to t_1).
    if v in L:
        current_best_duration = max(R[v][1] - R[v][0], 0)
        if max(new_arrival - new_start, 0) < current_best_duration:
            R[v] = (new_start, new_arrival)

        sv, av = L[v]
        sv = -sv
        if new_start > sv:
            L[v] = (-new_start, new_arrival)
            return True
        elif new_start == sv:
            return True
        else:
            # We did not update L : don't propagate the path !
            return False
    else:
        L[v] = (-new_start, new_arrival)
        R[v] = (new_start, new_arrival)
    return True


def FP_postprocess(R,
                   destination=None,
                   latencies=None):
    if destination is None:
        if latencies is None:
            latencies = {}
        for k, v in R.items():
            n = k[2]
            duration = v[1] - v[0]
            if duration < 0:
                duration = 0
            if n in latencies:
                latencies[n] = min(latencies[n], duration)
            else:
                latencies[n] = duration

    else:
        if latencies is None:
            latencies = math.inf
        if type(destination) == int:  # destination node
            for k, v in R.items():
                if k[2] == destination:
                    duration = v[1] - v[0]
                    if duration < 0:
                        duration = 0
                    latencies = min(latencies, duration)
        else:  # temporal destination node
            for k, v in R.items():
                if k[2] == destination[2] and k[0] <= destination[0] <= destination[1] <= k[1]:
                    duration = v[1] - v[0]
                    if duration < 0:
                        duration = 0
                    latencies = min(latencies, duration)

    return latencies





def FP_W_postprocess(W):
    latencies_pw = {}
    for n in W:
        latencies_pw[n] = FP_postprocess(W[n])
    return latencies_pw


def FP_pw(S):
    L_functions = {'initialisation': FP_initialisation,
                   'update': FP_update,
                   }
    W = F_Algorithm(S, L_functions)
    if W is not None:
        latencies_pw = FP_W_postprocess(W)
        return latencies_pw


######################################
#           Shortest Fastest Paths   #
######################################
def SFP_initialisation(source):
    L = defaultdict(SortedSet)
    L[source].add((-source[1], 0, source[0]))  # (-s_u,d_u,a_u)
    R = {}  # defaultdict(SortedSet)
    R[source] = (0, -source[1], source[0])  # (d_u,-s_u,a_u
    return L, R


def SFP_update(L, R, t0, t1, u, v):
    '''
    Update Lv from (t0,t1,u,v) according to SFP constraints.

    :param L:
    :param t0:
    :param t1:
    :param u:
    :param v:
    :return:
    '''
    su, du, au = L[u][0]
    su = -su
    best_duration_u = max(au - su, 0)
    Lbis = [(-su, du, au)]  #  We get elements that can be SFP or extended to be SFP
    du_old = du
    to_remove = []
    for cnt, e in enumerate(L[u][1:]):
        su, du, au = e
        su = -su
        if du >= du_old:
            to_remove.append(e)
        elif max(au - su, 0) == best_duration_u:
            Lbis.append((-su, du, au))
        else:  #  Impossible to extend to obtain a SFP
            to_remove += L[u][cnt + 1:]
            break
        du_old = du
    # for e in to_remove:
    #     L[u].remove(e)
    L[u] -= to_remove

    if v in L:
        current_best_duration = max(R[v][2] + R[v][1], 0)
        current_best_length = R[v][0]

    updated = False
    for su, du, au in Lbis:
        su = -su
        new_arrival = t0
        new_start = min(t1, su)
        new_distance = du + 1
        new_duration = max(new_arrival - new_start, 0)
        if v in L:
            sv, dv, av = L[v][0]
            sv = -sv

            if new_duration < current_best_duration or \
                    (new_duration == current_best_duration
                     and new_distance < current_best_length):
                R[v] = (new_distance, -new_start, new_arrival)
                current_best_duration = new_duration
                current_best_length = new_distance
                updated = True

            if new_start > sv or (new_start == sv and new_distance < dv) or \
                    (new_duration == 0 and new_distance < dv):
                L[v].add((-new_start, new_distance, new_arrival))
                updated = True
                # updated = True
        else:
            L[v].add((-new_start, new_distance, new_arrival))
            R[v] = (new_distance, -new_start, new_arrival)
            current_best_duration = new_duration
            current_best_length = new_distance
            updated = True
    return updated


def SFP_postprocess(R, destination=None,
                    latencies=None,
                    lengths=None):
    if destination is None:
        if latencies is None:
            latencies = {}
            lengths = {}
        for k, v in R.items():
            n = k[2]
            min_duration = max(v[2] + v[1], 0)
            min_length = v[0]

            # min_duration = max(min([i[3] + i[2] for i in v]),0)
            # min_length = min([i[1] for i in v if max(i[3] + i[2], 0) == min_duration])

            if n in latencies:
                if min_duration <= latencies[n]:
                    latencies[n] = min(latencies[n], min_duration)
                    lengths[n] = min(lengths[n], min_length)
            else:
                latencies[n] = min_duration
                lengths[n] = min_length

    else:
        if latencies is None:
            latencies = math.inf
            lengths = math.inf
        if type(destination) == int:  # destination node
            for k, v in R.items():
                if k[2] == destination:
                    min_duration = max(v[2] + v[1], 0)
                    min_length = v[0]

                    if min_duration <= latencies:
                        latencies = min(latencies, min_duration)
                        lengths = min(lengths, min_length)
        else:  # temporal destination node
            for k, v in R.items():
                if k[2] == destination[2] and k[0] <= destination[0] <= destination[1] <= k[1]:

                    min_duration = max(v[2] + v[1], 0)
                    min_length = v[0]

                    if min_duration <= latencies:
                        latencies = min(latencies, min_duration)
                        lengths = min(lengths, min_length)
    return latencies, lengths





def SFP_W_postprocess(W):
    latencies_pw = {}
    lengths_pw = {}
    for n in W:
        latencies_pw[n], lengths_pw[n] = SFP_postprocess(W[n])
    return latencies_pw, lengths_pw


def SFP_pw(S):
    L_functions = {'initialisation': SFP_initialisation,
                   'update': SFP_update,
                   }
    W = F_Algorithm(S, L_functions)
    if W is not None:
        latencies_pw, lengths_pw = SFP_W_postprocess(W)
        return latencies_pw, lengths_pw









