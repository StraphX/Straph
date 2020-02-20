import random,copy

class CTreapNode(object):
    '''
    A node of a Treap (priority is determined at random and used to balance the Treap)
    '''

    def __init__(self, data=None, parent=None, pred=None, suc=None, size=None,left=None,right=None):
        '''

        :param value:
        :param data:
        '''
        self.data = data
        self.priority = random.random()
        self.parent = parent
        self.left = left  # Left child
        self.right = right  # Right child
        self.pred = pred  # Predecessor (used in Euler Tour)
        self.suc = suc  # successor   (used in Euler Tour)
        self.size = size  # Used to count the number of nodes in the subtree rooted at the current node
        self.__tree_number = None # Used to catch a Tree in a forest, only accessible if the current node is root

    def __copy__(self):
        '''
        Only use to write data somewhere, only need data, left and right
        :return:
        '''
        data = copy.copy(self.data)
        left=right=None
        if self.left:
            left = copy.copy(self.left)
        if self.right:
            right = copy.copy(self.right)
        return CTreapNode(data=data,left=left,right=right)


    @property
    def tree_number(self):
        if self.parent is None:
            return self.__tree_number

    @tree_number.setter
    def tree_number(self,number):
        if self.parent is None:
            self.__tree_number = number


    def find_root(self):
        '''
        Find the root of the current node
        :param node:
        :return:
        '''
        current = self
        while current.parent:
            current = current.parent
        return current

    def left_rotation(self):
        '''
        Perform a left rotation on the Treap with the current TreapNode as the root
        https://en.wikipedia.org/wiki/Tree_rotation
        Note : This doesn't change the cyclic order, successor and predecessor unchanged
        :return: New root (aka the right child)
        '''
        root = self
        pivot = root.right

        # Change filiation
        pivot.parent = root.parent
        root.parent = pivot
        if pivot.left:
            pivot.left.parent = root

        # Change size
        root.size,pivot.size = pivot.size,root.size

        # Rotate
        root.right = pivot.left
        pivot.left = root
        root = pivot

        return root

    def right_rotation(self):
        '''
        Perform a right rotation on the Treap with the current TreapNode as the root
        https://en.wikipedia.org/wiki/Tree_rotation
        Note : This doesn't change the cyclic order, successor and predecessor unchanged
        :return: New root ( aka the left child)
        '''
        root = self
        pivot = root.left

        # Change filiation
        pivot.parent = root.parent
        root.parent = pivot
        if pivot.right:
            pivot.right.parent = root

        # Change size
        root.size,pivot.size = pivot.size,root.size

        # Rotate
        root.left = pivot.right
        pivot.right = root
        root = pivot

        return root

    def update_size(self):
        '''
        Used to keep the size
        :param node:
        :return:
        '''
        c = 1
        if self.left:
            c += self.left.size
        if self.right:
            c += self.right.size
        self.size = c

    def clear(self):
        self.parent = None
        self.left = None
        self.right = None
        self.pred = None
        self.suc = None

    def __repr__(self, depth=0, left_offset=0, right_offset=0):

        # print(" depth : ",depth," right offset : ",right_offset," left offset : ",left_offset)
        ret = "\t" * depth
        if self.parent:
            ret += " parent data : " + repr(self.parent.data)
        if self.suc:
            ret += " | suc data : " + repr(self.suc.data)
        if self.pred:
            ret += " | pred data : " + repr(self.pred.data)

        ret += " | node data : " + repr(
            self.data) + " priority : " + repr(
            self.priority) + "\n"

        if self.right:
            ret += "Right " + self.right.__repr__(depth=depth + 1, right_offset=right_offset + 1,
                                                  left_offset=left_offset)
        if self.left:
            ret += "Left " + self.left.__repr__(depth=depth + 1, left_offset=left_offset + 1,
                                                right_offset=right_offset)
        return ret




