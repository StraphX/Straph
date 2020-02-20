
from collections import defaultdict
from sortedcontainers import SortedSet

def get_pred_and_suc(value, times):
    '''
    Return the predecessor and the successor of the value as well as the predecessor of the predecessor of the value
    in the array times.
    :param value: float
    :param times: SortedSet
    :return:
    '''
    if len(times) == 2:
        return (times[0],times[0],times[1])
    id_before = times.bisect_left(value)
    if value == times[id_before]:
        pred = value
        suc = times[id_before+1]
        if id_before-1 <0:
            pred_pred = times[0]
        else:
            pred_pred = times[id_before-1]
    else:
        pred = times[id_before-1]
        suc = times[id_before]
        if id_before-2 <0:
            pred_pred = times[0]
        else:
            pred_pred = times[id_before-2]
    return (pred_pred, pred, suc)


# proof of correctness by induction on the size of E (sum of durations of links in it)
def dynamic_connectivity(E, SCC,times):
    while len(E) > 0:
        (b, e, u, v) = E.pop()
        assert u in times
        assert v in times
        # print("link :",(b,e,u,v))
        (tu_pred, tu, tu_suc) = get_pred_and_suc(b, times[u])    # Get predecessors and successor's times of u
        (Cu, bu, eu) = SCC[(tu, u)]                              # Get the current comp of u (Cu contains nodes of u's comp)
        (tv_pred, tv, tv_suc) = get_pred_and_suc(b, times[v])    # Get predecessors and successor's times of v
        (Cv, bv, ev) = SCC[(tv, v)]                              # Get the current comp of v (Cv contains nodes of v's comp)

        if SCC[(tu, u)] == SCC[(tv, v)]:
            # print("  1. ",u," and ",v," are in the same comp")
            # u and v are in the same component at time b
            # Manage link addition after end of component
            if e > tu_suc:
                E.add((tu_suc, e, u, v))

        elif tu < b or tv < b:
            # print("  2. split beginning")
            # SPLIT their beginning
            # Different components beginning before the link
            assert Cu != Cv
            if tu < b:
                for n in Cu:
                    SCC[(tu, n)] = (Cu, tu, b)
                    SCC[(b, n)] = (Cu, b, tu_suc)
                    times[n].add(b)
            if tv < b:
                for n in Cv:
                    SCC[(tv, n)] = (Cv, tv, b)
                    SCC[(b, n)] = (Cv, b, tv_suc)
                    times[n].add(b)
            E.add((b, e, u, v))

        elif e < tu_suc or e < tv_suc:
            # print("  3. split ending")
            # SPLIT their end
            # Different components with same beginning and the link ends before them
            assert Cu != Cv
            assert tu == b
            assert tv == b
            if e < tu_suc:
                for n in Cu:
                    SCC[(tu, n)] = (Cu, tu, e)
                    SCC[(e, n)] = (Cu, e, tu_suc)
                    times[n].add(e)
            if e < tv_suc:
                for n in Cv:
                    SCC[(tv, n)] = (Cv, tv, e)
                    SCC[(e, n)] = (Cv, e, tv_suc)
                    times[n].add(e)
            E.add((b, e, u, v))
        else:
            # MERGE COMPONENTS
            # different components with same beginning and the link ends after them
            assert Cu != Cv
            assert tu == b
            assert tv == b
            assert e >= tu_suc
            assert e >= tv_suc
            # print("  4. Merge ",Cu)
            # print("     with :",Cv)
            X = Cu.union(Cv)
            new = min(tu_suc, tv_suc)

            # Is the new component the same as the next one?
            if new < min(u[1],v[1]): # U and V cannot be in the next SCC if they aren't present
                if new == tu_suc:
                    (X_s, b_s, e_s) = SCC[(tu_suc, u)]
                if new == tv_suc:
                    (X_s, b_s, e_s) = SCC[(tv_suc, v)]
                if X == X_s:
                    # yes: move the component end
                    # print(" 4.1 next component equals current one")
                    assert (tv_suc == tu_suc)
                    for n in X:
                        SCC.pop((b_s, n))
                        times[n].remove(b_s)
                    new = e_s
            # else:
            #     print(" 4.1b No next component")

            # Is the new component the same as the previous one?
            if tu_pred < tu:
                (X_p, b_p, e_p) = SCC[(tu_pred, u)]
                assert b_p == tu_pred
                assert e_p == tu
                if X == X_p:
                    # yes: move the component beginning
                    # print(" 4.2 previous component equals current one")
                    for n in X:
                        # print("b :",b,"n :",n)
                        SCC.pop((b, n))
                        times[n].remove(b)
                    b = b_p
            # else:
            #     print(" 4.2b No previous component")

            # create the new component
            for n in X:
                SCC[(b, n)] = (X, b, new)
            # print("  5. New comp :")
            # print("  X :",X)
            # print("  begin : ",b," end : ",new)

            # add the end of the components if the link ends before them
            if tu_suc > new:
                for n in Cu:
                    SCC[(new, n)] = (Cu, new, tu_suc)
                    times[n].add(new)
            if tv_suc > new:
                for n in Cv:
                    SCC[(new, n)] = (Cv, new, tv_suc)
                    times[n].add(new)

            # continue with the rest of the link
            if e > new:
                E.add((new, e, u, v))
    return SCC


def strongly_connected_components_UF(S):
    # Initialisation SCC and times
    times = defaultdict(SortedSet)
    SCC = {}
    for n,np in zip(S.nodes,S.node_presence):
        for t0,t1 in zip(np[::2],np[1::2]):
            new_n = (t0,t1,n)
            SCC[(t0, new_n)] = ({new_n}, t0, t1)
            times[new_n].add(t0)
            times[new_n].add(t1)
    # Initialisation links
    L = S.augmented_ordered_links()
    E = set()
    for l in L:
        if l[0] == 1:
            b,e,u,v = l[1::]
            # print("link :",l[1::])
            assert u in times
            assert v in times
            E.add(tuple(l[1::]))

    # for v in S.nodes:
    #     SCC[(alpha, v)] = ({v}, alpha, omega)
    #     times[v].add(alpha)
    #     times[v].add(omega)
    # Initialisation avec le (minimum,maximum) event time for each node presence.
    SCC = dynamic_connectivity(E,SCC,times)

    # NO POSTPROCESSING FOR BENCHMARK
    #SCC = postprocess_SCC(SCC)
    return SCC



def postprocess_SCC(SCC):
    scc = []
    for k,v in SCC.items():
        (X, b, e) = v
        c = []
        # if (tuple(X),b,e) not in seen:
        #     seen.add((tuple(X),b,e))
        for n in X:
            if len(X) >1 or b!=e:
                c.append((b,e,n[2]))
        if c:
            scc.append(c)
        SCC[k] = None # Free Memory
    return scc