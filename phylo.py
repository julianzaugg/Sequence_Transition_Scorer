"""
This needs to be cleaned up. If I am going to use biopython, you might as well use it in your PhyloTree...instead make
PhyloTree an easier to use interface for BioPython trees.
"""

__author__ = 'julianzaugg'


from sequence import *
from Bio import Phylo


class PhyloTree:
    """Non-binary tree for representing phylogenetic relationships"""

    def __init__(self, root, **kwargs):
        self.nodes = dict()
        self.root = root
        self.label = None
        self.alignment = None

    def __getitem__(self, key):
        return self.nodes.get(key, None)

    def __len__(self):
        return len(self.nodes)

    def __iter__(self):
        return iter(self.nodes.values())

    def __str__(self):
        def recurse(root):
            if root.children:
                children = [recurse(child) for child in root.children]
                return '(%s)%s' % (','.join(children), root.label)
            else:
                return root.label
        tree = recurse(self.nodes[self.root]) + ";"
        return tree

    def _setLabel(self, label):
        self.label = label

    def addAlignment(self, aln):
        """
        Add an alignment to the tree.
        """
        self.alignment = aln

    def add(self, node):
        """
        Add a PhyloNode to the tree
        """
        self.nodes[node.label] = node

    def getNode(self, label):
        """
        Get a node object from the tree
        @return PhyloNode
        """
        return self.nodes[label]

    def getDescendants(self, node, transitive = False):
        """
        Get all descendants for a node. If not transitive, will return only immediate descendants
        @return list of descendants nodes
        """
        if not isinstance(node, PhyloNode):
            node = self.getNode(node)
        if node:
            return node.getDescendants(transitive)
        return []

    def getAncestors(self, node, transitive = False):
        """
        Get all ancestors for a node. If not transitive, will return only immediate ancestors
        @return list of ancestor nodes
        """
        if not isinstance(node, PhyloNode):
            node = self.getNode(node)
        if node:
            myroot = self.getNode(self.root)
            found = False
            branching = []
            if node == myroot:
                found = True
            while not found and myroot != None:
                branching.append(myroot)
                for child in myroot.children:
                    if child == node:
                        found = True
                        break
                    if child.isAncestorOf(node, transitive = True):
                        myroot = child
            if found and transitive:
                return branching
            elif found and len(branching) > 0:
                return [branching[len(branching)-1]]
            return []

    def getLeaves(self):
        """
        Get the leaves of the tree
        @return set of nodes without children
        """
        return {node for node in self if not node.children}

    def getTopologicalOrdering(self):
        """
        Get the ordering of all nodes in the tree based on topological order
        @return list of nodes in topological order
        """
        ordered = []
        leaves = self.getLeaves()
        visited = set()

        def visit(node):
            visited.add(node)
            for parent in self.getAncestors(node):
                if parent not in visited:
                    visit(parent)
            ordered.append(node)
        for node in leaves:
            visit(node)
        return ordered

    def getBranch(self, node):
        """
        @return ordered list of nodes in the same branch for queried node (including query node)
        """
        return self.getAncestors(node, transitive=True) +\
               [node] +\
               self.getDescendants(node, transitive=True)

    def getSubBranch(self, ancestor_node, descendant_node):
        """
        @return a ordered list of nodes between (and including) ancestor and descendant nodes
        """
        node_branch = self.getBranch(descendant_node)
        ancestor_index = node_branch.index(ancestor_node)
        descendant_index = node_branch.index(descendant_node)
        return node_branch[ancestor_index:descendant_index + 1]

    def findPath(self, start_node, end_node, path = []):
        """
        Find the path between two nodes in the tree
        @return list of nodes in order of steps in path
        """
        path = path + [start_node]
        if start_node == end_node:
            return path
        immediates = self.getDescendants(start_node) + self.getAncestors(start_node)
        for node in immediates:
            if node not in path:
                newpath = self.findPath(node, end_node, path)
                if newpath: return newpath
        return None

    def getBranches(self):
        """
        Get all complete branches in the tree
        @return dictionary of branches for the tree
        """
        leaves = self.getLeaves()
        branch_node_sets = {}
        for leaf in leaves:
            branch_nodes = self.getAncestors(leaf, transitive = True) + [leaf]
            branch_node_sets[leaf] = branch_nodes
        return branch_node_sets

class PhyloNode:

    def __init__(self, label, **kwargs):
        self.label = label
        self.children = []
        self.info = dict()
        self.parent = None

    def getDescendants(self, transitive = False):
        """
        Get all descendants for this node. If not transitive, will return only immediate descendants
        @return list of descendant nodes
        """
        children = list(self.children)
        if not transitive:
            return children
        else:
            grandchildren = []
            for child in children:
                grandchildren.extend(child.getDescendants(transitive))
            children.extend(grandchildren)
            return children

    def assignInfo(self, **kwargs):
        """
        Assign additional annotations to this node
        """
        for k,v in kwargs.items():
            self.info[k] = v

    def addChild(self, childNode):
        if not isinstance(childNode, PhyloNode):
            raise TypeError("supplied child node %s is node of type PhyloNode" % childNode)
        self.children.append(childNode)

    def isAncestorOf(self, node, transitive = False):
        """
        Determine if this node is the ancestor of the specified node.
        If transitive is True (default), all descendants are included.
        If transitive is False, only direct descendants are included.
        @return boolean
        """
        if node in self.children:
            return True
        elif transitive:
            for child_node in self.children:
                status = child_node.isAncestorOf(node, transitive)
                if status:
                    return True
        else:
            return False

    def __str__(self):
        return self.label


#This is a VERY dirty Hack. Fix later.
def load_tree(filename):
    """
    Use Biopython to load newick tree. This is then converted into
    our own Phylo tree format.
    """
    bio_tree = Phylo.read(filename, 'newick')
    tree = PhyloTree(bio_tree.root.name)
    def thing(root):
        if len(root) == 0:
            tree.add(PhyloNode(root.name))
        else:
            tree.add(PhyloNode(root.name))
            # print root.name
            tree[root.name].children = [child.name for child in root]
            # if root.name == None: print tree[root.name].children
            for rc in root:
                thing(rc)
    thing(bio_tree.root)
    for node in tree.nodes.values():
        new_children = []
        for child in node.children:
            new_children.append(tree.nodes[child])
        node.children = new_children
    return tree
