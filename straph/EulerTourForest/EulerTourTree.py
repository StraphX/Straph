import matplotlib.pyplot as plt
import matplotlib.collections as mcol
import msgpack,math
import copy

from collections import defaultdict
from straph.EulerTourForest.Chained_Treap import CTreapNode


def euler_tour_from_edge_list(edge_list):
    '''
    Compute the euler tour of a graph ( with edges, including self loops, instead of nodes).

     https://en.wikipedia.org/wiki/Euler_tour_technique
    :param edge_list: An edge list.
    :return: A list containing the euler tour.
    '''
    a_l = defaultdict(list)
    for l in edge_list:  # Create edge list (2*m links)
        u, v = l
        a_l[u].append(v)
        a_l[v].append(u)
    tour = []
    current_node = edge_list[0][0]
    queue = [current_node]
    while queue:
        if a_l[current_node]:
            queue.append(current_node)
            current_node = a_l[current_node].pop()
        else:
            tour.append(current_node)
            current_node = queue.pop()
    edge_tour = []
    seen = set()
    for i in range(len(tour) - 1):
        u = tour[i]
        if u not in seen:
            seen.add(u)
            edge_tour.append((u, u))
        edge_tour.append((u, tour[i + 1]))
    return edge_tour

def construct_euler_tour_tree(edge_list):
    '''
    Construct an Euler Tour Tree from an edge list
    :param edge_list: An edge list
    :return: An euler tour tree
    '''
    euler_tour = euler_tour_from_edge_list(edge_list)
    ETT = EulerTourTree()
    tree_edge_2_node = defaultdict(list)
    for i, n in enumerate(euler_tour):
        node = ETT.insert(data=n, inlast=True)
        if (n[1], n[0]) in tree_edge_2_node:
            tree_edge_2_node[(n[1], n[0])].append(node)
        else:
            tree_edge_2_node[n].append(node)
    return ETT,tree_edge_2_node


class EulerTourTree(object):
    '''
    https://en.wikipedia.org/wiki/Euler_tour_technique
    '''
    def __init__(self, root=None, first=None, last=None, weight=None,begin_time=None,
                 end_time = None):
        '''

        :param root:
        :param first:
        :param last:
        :param weight:
        '''
        self.root = root
        self.first = first
        self.last = last
        self.weight = weight  # Sum of non tree edges adjacents to edge in tree # TODO: see self.replace()
        self.begin_time = begin_time
        self.end_time = end_time


    def __iter__(self):
        yield self.first
        current = self.first.suc
        while current != self.first:
            yield current
            current = current.suc

    def __copy__(self):
        '''
        Only use to write data, only need data of root and recursively
        :return:
        '''
        return EulerTourTree(root=copy.copy(self.root),begin_time = self.begin_time,end_time=self.end_time)

    def __repr__(self):
        rep = "Euler Tour : "+repr(self.get_euler_tour())+"\n"
        rep += "Priority Order :"+str(self.get_data_in_priority_order())+"\n"
        rep += str(self.root)
        return rep


    def get_data_in_priority_order(self):
        L = []
        self._get_data_in_priority_order(self.root, L)
        return L

    def _get_data_in_priority_order(self, node, L):
        if node:
            L.append(node.data)
            if node.right:
                self._get_data_in_priority_order(node.right, L)
            if node.left:
                self._get_data_in_priority_order(node.left, L)

    def get_internal_structure(self):
        return self._get_internal_structure(self.root)

    def _get_internal_structure(self, node, N=None, E=None, x_parent=0, y_parent=0):
        '''
        Return (x_pos,y_pos==depth)
        :param node:
        :param depth:
        :return:
        '''
        if not N:
            N = []
        if not E:
            E = []
        N.append(((x_parent, y_parent), (node.data,node.size)))  # ,node.priority))) # priority
        if node.left:
            if y_parent:
                y_pos = y_parent - 1
                offset = x_parent / y_pos
                if offset > 0:
                    x_pos = x_parent - offset
                else:
                    x_pos = x_parent + offset
            else:
                x_pos = -10
                y_pos = -1

            E.append(((x_parent, y_parent), (x_pos, y_pos)))
            self._get_internal_structure(node.left, N=N, E=E, x_parent=x_pos, y_parent=y_pos)
        if node.right:
            if y_parent:
                y_pos = y_parent - 1
                offset = x_parent / y_pos
                if offset > 0:
                    x_pos = x_parent + offset
                else:
                    x_pos = x_parent - offset
            else:
                x_pos = 10
                y_pos = -1
            E.append(((x_parent, y_parent), (x_pos, y_pos)))
            self._get_internal_structure(node.right, N=N, E=E, x_parent=x_pos, y_parent=y_pos)
        return N, E


    def get_euler_tour(self):
        '''
        Return the induced euler tour representation of the Tree
        :return:
        '''
        first = self.first
        euler_tour = [first.data]
        current = first.suc
        cnt =0
        while current != first:
            euler_tour.append(current.data)
            current = current.suc
            cnt += 1
        return euler_tour

    def check_heap_invariant(self):
        '''
        Check the heap invariant of the Treap
        :return: True if invariant respected, false otherwise
        '''
        return self._check_heap_invariant(self.root)

    def _check_heap_invariant(self, node):
        if node.left:
            assert node.priority <= node.left.priority
            assert self._check_heap_invariant(node.left)
        if node.right:
            assert node.priority <= node.right.priority
            assert self._check_heap_invariant(node.right)
        return True

    def plot(self, title=None):
        N, E = self.get_internal_structure()
        fig, ax = plt.subplots()
        y_min = 0
        x_min = 0
        x_max = 0
        for pos, data in N:         # Nodes are ordered as in a DFS traversal of the tree
            label = str(data[0])+ "\n" + str(data[1])   # Data + Size
            x_max = max(x_max, pos[0])
            x_min = min(x_min, pos[0])
            y_min = min(y_min, pos[1])
            ax.text(pos[0], pos[1], label, color='#2d5986',
                    bbox=dict(facecolor='none', edgecolor='#2d5986', boxstyle='round,pad=1'))

        # Set limit
        edge_collections = mcol.LineCollection(E, colors=['#2d5986'], linewidths=2, alpha=0.5)

        ax.add_collection(edge_collections)
        ax.set_ylim(y_min - 1, 1)
        ax.set_xlim(x_min - 1, x_max + 1)
        if title:
            ax.set_title(title)


    def _balance_down(self, node):
        if node.left:
            node.left = self._balance_down(node.left)
            if node.left.priority < node.priority:
                node = node.right_rotation()
        if node.right:
            node.right = self._balance_down(node.right)
            if node.right.priority < node.priority:
                node = node.left_rotation()
        return node

    def balance_down(self):
        '''
        Balance the Treap to respect the Heap invariant, from root to leaves (full browsing)
        :Status: OK
        :return:
        '''
        self.root = self._balance_down(self.root)

    def insert(self, where=None, data=None, inlast=False,priority = None):

        if not self.root:
            node = CTreapNode(data=data, size=1)
            self.root = node
            self.first = node
            self.last = node
            return node
        elif inlast:  # It means that this node is the new last
            node = self._insert(where=self.last, data=data,priority=priority)
            self.last.suc = node
            self.last = node
            self.last.suc = self.first
            self.first.pred = node
        else:
            node = self._insert(where=where, data=data,priority=priority)
        # # UPDATE SIZE OF PARENTS
        # p = node.parent
        # while p:
        #     p.update_size()
        #     p = p.parent
        return node


    def _insert(self, where, data=None,priority=None):
        '''
        (TODO) different balance : just from the children to the root
        Insert a node in the Treap after *where* (a node of the CTreap).
        Idea : If *where* doesn't have any right children, given the fact that the current node
        come after *where*, we can add it directly, it respects the order.
        Otherwise the idea is to put it just before the node that come after *where*
        so as the left children of the successor of *where*.
        :param key:
        :param data:
        :param priority:
        :param suc:
        :param pred:
        :return:
        '''
        if not where.right:
            node = CTreapNode(data=data, size=1, parent=where, pred=where, suc=where.suc)
            where.right = node
            if where.suc:
                where.suc.pred = node
            where.suc = node
            if priority is not None:
                node.priority=priority
            self.balance_down()
            return node
        suc = where.suc
        node = CTreapNode(data=data, size=1, parent=suc, pred=where, suc=suc)
        if suc.left:
            print("Noeud left a dejà un voisin, bizarre")
            raise ValueError
        suc.left = node
        suc.pred = node
        where.suc = node
        if priority is not None:
            node.priority = priority
        self.balance_down()
        return node

    def find_root(self, node):
        '''
        Find the root of the current node
        :param node:
        :return:
        '''
        current = node
        while current.parent:
            current = current.parent
        return current

    def find_first(self):
        current = self.root
        while current.left:
            current = current.left
        return current

    def find_last(self):
        current = self.root
        while current.right:
            current = current.right
        return current



    def remove(self,node):
        '''
        Remove the node
        :param node: Node to remove
        :return:
        '''
        if node == self.first:
            self.first = node.suc
        if node == self.last:
            self.last = node.pred
        if node.suc:
            node.suc.pred = node.pred
        if node.pred:
            node.pred.suc = node.suc
        if node == self.root :
            # print("Our node is a root")
            self.root = self._remove(node)
            if self.root:
                self.root.parent = None
        elif node.parent.left == node:
            # print("Our node is a Left child")
            node.parent.left = self._remove(node)
            # # UPDATE SIZE OF PARENTS
            # p = node.parent
            # while p:
            #     p.update_size()
            #     p = p.parent
        else:
            # print("Our node is a right child")
            node.parent.right = self._remove(node)
            # # UPDATE SIZE OF PARENTS
            # p = node.parent
            # while p:
            #     p.update_size()
            #     p = p.parent



    def _remove(self,node):
        if not node.left and not node.right:
            return None
        elif not node.left:
            if node.parent:
                node.right.parent = node.parent
            return node.right
        elif not node.right:
            if node.parent:
                node.left.parent = node.parent
            return node.left
        else:
            if node.left.priority < node.right.priority:
                # print("Right rotation")
                node = node.right_rotation()    # Rotation already deals with filiation
                node.right = self._remove(node.right)
            else:
                # print("Left rotation")
                node = node.left_rotation()     # Rotation already deals with filiation
                node.left = self._remove(node.left)
        return node

    def split(self, where):
        '''
        Split the Euler Tour Tree according to the split in its Euler Tour.
        [first,...,where,after_where,...,last]
        ->
        [first,...,where] [after_where,...last]

        Note : The left subtree contains at least the node *where* whereas the right subtree can be empty.
        :param where: The split is effectuated just after *where*
        :return: Left subtree and Right subtree
        '''
        # print(" Split on :", where.data)
        after_where = where.suc
        first = self.first
        last = self.last

        s = self.insert(where,priority=0)

        T_left = s.left
        T_right = s.right

        T_left.parent = None
        T_left = EulerTourTree(root=T_left)
        T_left.first = first
        first.pred = where
        T_left.last = where
        where.suc = first

        if T_right:
            T_right.parent = None
            T_right = EulerTourTree(root=T_right)
            T_right.first = after_where
            after_where.pred = last
            T_right.last = last
            last.suc = after_where

        return T_left, T_right

    def releaf(self, where):
        L, R = self.split(where=where)
        E = union_treap(R, L)
        return E

    def cut(self, nodes):
        '''
        Remove and edge from the Euler Tour Tree.
        [first,...,node_0,after_node_0,....,node_1,after_node_1,...,last]
        ->
        [first,...,before_node_0,node_0] [after_node_0,...,node_1,last]
        ->
        [first,...,before_node_0] [after_node_0,...,before_node_1] [after_node_1,...,last]
        ->
        [first,...,before_node_0,after_node_1,...,last] [after_node_0,before_node_1]

        :param e: an edge
        :return:
        '''
        # print("  Nodes :", [i.data for i in nodes])
        # Remove the first occurence of the link :
        J, K = self.split(nodes[0])
        # print(" J :\n", J)
        # print(" K :\n", K)
        # if J.first != nodes[0]:
        J.remove(nodes[0])
        # else:
        #     # It means that it only remain nodes[0] in J
        #     J = None
        # print("\n Remove first occurence : ", nodes[0].data)
        # print(" J :\n", J)
        # print(" K :\n", K)

        # J.plot(" Left after removal of " + repr(nodes[0].data))
        # K.plot(" Right after removal of " + repr(nodes[0].data))

        #  Remove the second occurence of the link
        if K and nodes[1].find_root() == K.root:
            # nodes[1] is in K
            K, L = K.split(nodes[1])
            # print(" K :\n", K)
            # print(" L :\n", L)
            # if K.first != nodes[1]:
            K.remove(nodes[1])
            # else:
                # It means that it only remain nodes[1] in K
                # K =None
            # print("\n Remove second occurence : ", nodes[1].data)
            # print(" K :\n", K)
            # print(" L :\n", L)
            # K.plot(" Left after removal of " + repr(nodes[1].data))
            E1 = K
            E2 = union_treap(J, L)
        else:
            # nodes[1] is in J
            J, L = J.split(nodes[1])
            # print(" J :\n", J)
            # print(" L :\n", L)
            # if J.first != nodes[1]:
            J.remove(nodes[1])
                # J.plot(" Left after removal of " + repr(nodes[1].data))
            # else:
                # It means that it only remain nodes[1] in J
                # J = None
            # print("\n Remove second occurence : ", nodes[1].data)
            # print(" J :\n", J)
            # print(" L :\n", L)

            E1 = union_treap(J, K)
            E2 = L

        # print("  E1 after cut : \n", E1)
        # E1.plot("E1 after cut " + repr(nodes[0].data))
        #
        # print("  E2 : \n", E2)
        # E2.plot("E2 after cut " + repr(nodes[1].data))
        return E1, E2


def union_treap(T1, T2):
    if not T2 or not T2.root:
        return T1
    if not T1 or not T1.root:
        return T2

    first = T1.first
    last = T2.last
    T1, T2 = T1.root, T2.root
    # We get the right most leaf of T1
    rl = T1
    while rl.right:
        rl = rl.right

    # We get the left most leaf of T2
    ll = T2
    while ll.left:
        ll = ll.left

    # We concat the induced euler tour
    # Maximum of first tree with minimum of second tree
    rl.suc = ll
    ll.pred = rl
    new_root = _union_treap(T1, T2)
    new_root.parent = None

    T = EulerTourTree(root=new_root)
    T.first = first
    T.first.pred = last
    T.last = last
    T.last.suc = first
    return T


def _union_treap(T1, T2):
    if not T1:
        return T2
    elif not T2:
        return T1
    elif T1.priority < T2.priority:
        T1.right = _union_treap(T1.right, T2)
        T1.right.parent = T1
        # T1.update_size()
        return T1
    else:
        T2.left = _union_treap(T1, T2.left)
        T2.left.parent = T2
        # T2.update_size()
        return T2






