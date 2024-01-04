"""
Contains Tree and Node classes that together can read in phylogenetic trees
in Newick format, manipulate them, and print them.
"""
import copy
import pydot

class Node:
    """
    Class to keep information about nodes.
    """
    def __init__(self, name, idx, distance_to_parent=None, children=None, parent=None):
        self.name = name
        self.parent = parent
        self.children = children
        self.idx = idx

        self.distance_to_parent = distance_to_parent

        if self.children is None:
            self.children = []

class Tree:
    """
    Class to keep info about the tree. Reads in from a string of a file.
    """
    def __init__(self, newick_string="", file_name=""):
        self.newick_string = newick_string

        # read from file if preferred
        if self.newick_string == "":
            if file_name == "":
                print("ERROR: you must provide either a Newick string or a file name.")
                raise ValueError
            with open(file_name, 'r') as f:
                self.newick_string = f.readline().strip()

        self.root = None
        self.number_of_nodes = 0

        self.generate_from_string()

    def generate_from_string(self):
        """
        Function to generate tree from Newick string.
        """
        if self.newick_string[-1] != ';':
            print('ERROR: Newick strings must end with ";".')
            raise ValueError

        if self.newick_string.count('(') != self.newick_string.count(')'):
            print('ERROR: unbalanced parantheses in Newick string.')
            raise ValueError

        right_ptr = len(self.newick_string) - 1
        curr = None
        parent = None
        prev_char = ''
        idx = 0
        while right_ptr >= 0:
            char = self.newick_string[right_ptr]
            if char == ';':
                prev_char = char
                name = ''
            elif char in '(,)':
                # special things will happen
                if prev_char == ';':
                    # we read in the root
                    root = Node(name, idx, parent=parent)
                    idx += 1
                    curr = root
                elif prev_char == ')' and char in ",)":
                    # create child
                    if ':' in name:
                        n = name.split(':')[0]
                        dist = float(name.split(':')[1])
                    else:
                        n = name
                        dist = None
                    parent = curr
                    curr = Node(n, idx, parent=parent, distance_to_parent=dist)
                    idx += 1
                    parent.children.append(curr)
                elif prev_char == ',':
                    # add sibling
                    if ':' in name:
                        n = name.split(':')[0]
                        dist = float(name.split(':')[1])
                    else:
                        n = name
                        dist = None
                    curr = Node(n, idx, parent=parent, distance_to_parent=dist)
                    idx += 1
                    parent.children.append(curr)
                elif prev_char == '(':
                    # go up one level
                    curr = parent
                    parent = curr.parent
                prev_char = char
                name = ''
            else:
                name = char + name

            right_ptr -= 1

        self.root = root
        self.number_of_nodes = idx + 1

    def is_adjacent(self, node1, node2):
        """
        Given two nodes, find out if they are adjacent.
        """
        if node1 in node2.children:
            assert node1.parent == node2
            return True

        if node2 in node1.children:
            assert node2.parent == node1
            return True

        return False
    
    def add_node_to_edge(self, node_name, edge):
        """
        Given a new node name, and an edge, add a new node to the edge.
        """
        u = edge[0]
        v = edge[1]
        if isinstance(u, str):
            u = self.get_node_with_name_bfs(u)
        if isinstance(v, str):
            v = self.get_node_with_name_bfs(v)

        if not self.is_adjacent(u, v):
            raise ValueError

        node = Node(name=node_name, idx=self.number_of_nodes)

        if v.parent == u:
            u, v = v, u

        u.parent = node
        dist = u.distance_to_parent
        if dist is not None:
            u.distance_to_parent = dist / 2
            node.distance_to_parent = dist / 2
        node.parent = v
        node.children.append(u)
        v.children.remove(u)
        v.children.append(node)
        
        self.number_of_nodes += 1

    def get_node_with_name_bfs(self, target_name):
        """
        Given a node name, find the first node object that matches it
        using breadth-first search.
        
        NOTE: node names need not be unique. This function returns the first one it
        encounters in a breadth-first search.
        """
        Q = [self.root]
        visited = []

        while Q:
            node = Q.pop(0)
            if node.name == target_name:
                return node
            visited.append(node.idx)
            for child in node.children:
                if child.idx not in visited:
                    Q.append(child)

        return None

    def change_names(self, name_dicts, remove_not_in_dict=True):
        """
        Given a dictionary of old_names -> new_names, change names of nodes.
        """
        if not isinstance(name_dicts, dict):
            raise TypeError

        Q = [self.root]
        visited = []

        while Q:
            node = Q.pop(0)
            n = node.name
            if remove_not_in_dict:
                if n not in name_dicts:
                    node.name = ''
                else:
                    node.name = name_dicts[n]
            else:
                node.name = name_dicts[n]
            visited.append(node.idx)
            for child in node.children:
                if child.idx not in visited:
                    Q.append(child)

    def get_node_with_idx_bfs(self, idx):
        """
        Given a node index, use BFS to find the first node that matches it.
        """
        Q = [self.root]
        visited = []
        while Q:
            node = Q.pop(0)
            if node.idx == idx:
                return node
            visited.append(node.idx)
            for child in node.children:
                if child.idx not in visited:
                    Q.append(child)

        return None


    def re_root(self, new_root):
        """
        Function to change the root of the tree using an existing node.
        """
        if isinstance(new_root, str):
            new_root = self.get_node_with_name_bfs(new_root)
            if new_root is None:
                print(f"ERROR: node with name {new_root.name} not found.")
                raise ValueError
        elif isinstance(new_root, int):
            new_root = self.get_node_with_idx_bfs(new_root)
        elif not isinstance(new_root, Node):
            raise TypeError
        if new_root is None:
            print("ERROR: node not found.")
            raise ValueError

        old_root = self.root

        curr = new_root
        prev = None
        old_dist = None
        while curr != old_root:
            parent = curr.parent
            curr.children.append(parent)
            curr.parent = prev
            old_dist, curr.distance_to_parent = curr.distance_to_parent, old_dist
            parent.children.remove(curr)
            prev = curr
            curr = parent
        parent.parent = prev
        curr.distance_to_parent = old_dist

        self.root = new_root

    def render_graph(self, out_name, format='png'):
        """
        Function to render the graph to an image file using pydot.
        """
        graph = pydot.Dot('test', graph_type='graph', rankdir='RL')

        stack = [self.root]
        visited = []
        while stack:
            node = stack.pop()
            visited.append(node.idx)
            if node.parent is None:
                graph.add_node(pydot.Node(f"n_{node.idx}", shape='hexagon',
                                            label=node.name))
            elif len(node.children) > 0:
                graph.add_node(pydot.Node(f"n_{node.idx}", shape='diamond',
                                             label=node.name,
                                             fontcolor='blue'))
            else:
                graph.add_node(pydot.Node(f"n_{node.idx}", shape='box', label=node.name,
                                            style='filled', fillcolor='black',
                                            fontcolor='white'))

            if node.parent is not None:
                if node.distance_to_parent is None:
                    label = ""
                else:
                    label = f"{node.distance_to_parent:.1e}"
                graph.add_edge(pydot.Edge(f"n_{node.parent.idx}", f"n_{node.idx}", label=label))
            for child in node.children:
                if child.idx in visited:
                    continue

                stack.append(child)

        graph.write(out_name, format=format)

    def get_leaves_in_subtree(self, node):
        """
        Given a node, find all leaves in the subtree where node = root.
        """
        if isinstance(node, str):
            node = self.get_node_with_name_bfs(node)
        elif isinstance(node, int):
            node = self.get_node_with_idx_bfs(node)
        elif not isinstance(node, Node):
            raise TypeError
        if node is None:
            raise ValueError

        visited = []
        leaves = []
        stack = [node]
        while stack:
            node = stack.pop()
            visited.append(node.idx)

            if len(node.children) == 0:
                leaves.append(node)

            for child in node.children:
                if child.idx not in visited:
                    stack.append(child)

        return leaves

    def prune(self, min_support=80, split_char='/'):
        """
        Function that returns number of nodes to remove so that all nodes in the present tree will have
        at least min_support support.

        (Assuming support is contained in node names)
        """
        leaves = self.get_leaves_in_subtree(self.root)

        to_remove = []
        visited = []
        Q = copy.copy(leaves)

        while Q:
            node = Q.pop(0)

            if node.idx in visited:
                continue

            visited.append(node.idx)
            parent = node.parent
            if parent == self.root or parent is None or parent.name == '':
                continue

            support = min([float(s) for s in parent.name.split('/')])
            if support >= min_support:
                Q.append(parent)
                leaves.append(parent)
                continue

            siblings = [n for n in parent.children if n != node]
            leaf_siblings = [n for n in siblings if n in leaves]
            rem = [node] + leaf_siblings
            visited += [n.idx for n in leaf_siblings]

            to_remove.append((len(siblings), rem))
            if len(leaf_siblings) == len(siblings):
                Q.append(parent)
                leaves.append(parent)

        return to_remove

    def get_leaves_to_remove(self, min_support=80, split_char='/'):
        """
        Function that determines what leaves to remove so that the present
        tree has min_support.
        """
        to_remove = self.prune(min_support=min_support, split_char=split_char)

        # expand branches
        new_list = []
        for el in to_remove:
            new = []
            for node in el[1]:
                if len(node.children) != 0:
                    new.append([n for n in self.get_leaves_in_subtree(node)])
                else:
                    new.append([node])

            new_list.append((el[0], new))

        # make passes to remove all you can, unambiguously
        removed = []
        added = True
        while added:
            added = False
            new = []
            for el in new_list:
                if el[0] == len(el[1]):
                    # remove them all
                    r = [item for sublist in el[1] for item in sublist]
                    for item in r:
                        removed.append(item)
                        added = True
                else:
                    # see if we have already removed some stuff here
                    new_names = []
                    for name in el[1]:
                        nn = [n for n in name if n not in removed]
                        if len(nn) > 0:
                            new_names.append(nn)
                    if len(new_names) > 0:
                        change = len(new_names) - len(el[1])
                        if el[0] + change > 0:
                            new.append((el[0],new_names))
            new_list = copy.copy(new)

        return removed, new_list

    def get_supports(self, split_char='/'):
        """
        Traverse the tree to get a list of all supports.
        """
        supports = []
        Q = [self.root]
        visited = []
        while Q:
            node = Q.pop(0)
            if len(node.children) != 0 and node.name != '':
                support = min([float(s) for s in node.name.split('/')])
                supports.append(support)
            visited.append(node.idx)
            for child in node.children:
                if child not in visited:
                    Q.append(child)
        return supports

    def write_newick(self, no_supports=False, no_lengths=False, no_names=False):
        """
        Given the current tree, which may be modified from the original
        Newick string, generate a new Newick string.
        """
        stack = [(0, self.root)]
        prev_depth = -1
        while stack:
            depth, node = stack.pop(-1)
            name = self.get_name_formatted(node, no_support=no_supports, no_length=no_lengths, no_name=no_names)
            if depth == 0:
                # root node
                new_string = name + ";"
            else:
                if depth == prev_depth:
                    new_string = name +  ',' + new_string
                elif depth > prev_depth:
                    new_string = name + ')' + new_string
                else:
                    new_string = name + ',' + '('*(prev_depth - depth) + new_string

            for child in node.children[::-1]:
                stack.append((depth+1, child))
            prev_depth = depth

        new_string = '(' * depth + new_string

        return new_string

    def get_name_formatted(self, node, no_support=False, no_length=False, no_name=False):
        """
        Given a node, produce a string that reads

        A:0.3

        where 0.3 is distance to parent, to be printed only when no_length==False,
        A is the node name, which may or may not be printed depending on no_name and no_support.
        """
        if no_name:
            name = ""
        else:
            if no_support:
                if len(node.children) > 0:
                    name = ""
                else:
                    name = node.name
            else:
                name = node.name
            
        if no_length or node.distance_to_parent is None:
            return name

        if node == self.root:
            return name

        return f"{name}:{node.distance_to_parent}"
