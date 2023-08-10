import pytest
from sequence_analysis.tree import Tree

import pdb

class TestTree:
    def test_gen(self):
        tree_str = "((A,B)C,(D,E)F,H)G;"
        tree = Tree(tree_str)

        assert tree.root is not None
        assert tree.root.name == 'G'
        assert len(tree.root.children) == 3

    def test_gen_2(self):
        tree_str = "(A,((B,C)D,E)F,G)H;"
        tree = Tree(tree_str)
        assert tree.root is not None
        assert len(tree.root.children) == 3
        

    def test_gen_3(self):
        tree_str = "(,(((((((((((((((((((,),((,),)),(,)),),),),),),),),(((((((((((,),),(,)),(,)),),),(,)),),),),)),((,),(,))),(,)),(((,),(,)),)),),),((,),)),),),);"

        tree = Tree(tree_str)

        assert tree.root is not None
        assert tree.root.name == ''
        assert len(tree.root.children) == 3

    def test_is_adjacent(self):
        tree_str = "((A,B)C,(D,E)F,H)G;"
        tree =Tree(tree_str)

        assert tree.is_adjacent(tree.root, tree.root.children[0])

        assert not tree.is_adjacent(tree.root, tree.root.children[1].children[0])

    def test_get_node_with_name_bfs(self):
        tree_str = "((A,B)C,(D,E)F,H)G;"
        tree = Tree(tree_str)
        node = tree.get_node_with_name_bfs('E')

        assert node.name == 'E'

    def test_get_node_with_idx_bfs(self):
        tree_str = "((A,B)C,(D,E)F,H)G;"
        tree = Tree(tree_str)
        node = tree.get_node_with_idx_bfs(1)

        assert node.name == 'H'

    def test_re_root(self):
        tree_str = "((A,B)C,(D,E)F,H)G;"
        tree = Tree(tree_str)

        tree.re_root(tree.root.children[0])
        assert tree.root.name == "H"
        assert len(tree.root.children) == 1
        assert 'G' in [n.name for n in tree.root.children]

    def test_re_root_2(self):
        tree = Tree("(((A,B)C,D,X)E,((F,G)Y),Z)H;")
        new_root = tree.get_node_with_name_bfs('C')
        tree.re_root(new_root)
        
        assert tree.root.name == 'C'
        assert tree.get_node_with_name_bfs('E') in tree.root.children
        assert tree.root.parent is None
        assert tree.get_node_with_name_bfs('E').parent == tree.root

    def test_add_node_to_edge(self):
        tree = Tree("(((A,B)C,D,X)E,((F,G)Y),Z)H;")
        tree.add_node_to_edge('new',['A','C'])

        new_node = tree.get_node_with_name_bfs('new')
        assert new_node.parent == tree.get_node_with_name_bfs('C')
        assert new_node in tree.get_node_with_name_bfs('C').children

    def test_write(self):
        tree = Tree("(((A,B)C,D,X)E,((F,G)Y),Z)H;")
        tree.render_graph("test.png")

        tree.add_node_to_edge('new', ['A','C'])
        tree.render_graph('test_node_added.png')

        tree.re_root('new')
        tree.render_graph('test_new_root.png')

    def test_distances(self):
        tree = Tree("(((A:0.1,B:0.5)C:0.12,D:0.55,X:0.01)E:0.4,((F:0.8,G:0.7)Y:0.3)N:0.1,Z:0.2)H;")

        tree.render_graph('test_distances.png')

        node = tree.get_node_with_name_bfs('A')
        assert node.distance_to_parent == 0.1

    def test_distances_2(self):
        tree = Tree("(((A:0.1,B:0.5)C:0.12,D:0.55,X:0.01)E:0.4,((F:0.8,G:0.7)Y:0.3)N:0.1,Z:0.2)H;")

        tree.add_node_to_edge('new', ['B','C'])
        assert tree.get_node_with_name_bfs('new').distance_to_parent == 0.25
        assert tree.get_node_with_name_bfs('B').distance_to_parent == 0.25
        assert tree.get_node_with_name_bfs('C').distance_to_parent == 0.12

    def test_distances_3(self):
        tree = Tree("(((A:0.1,B:0.5)C:0.12,D:0.55,X:0.01)E:0.4,((F:0.8,G:0.7)Y:0.3)N:0.1,Z:0.2)H;")

        tree.re_root('C')

        tree.render_graph('test_reroot_distances.png')
        assert tree.get_node_with_name_bfs('C').distance_to_parent is None
        assert tree.get_node_with_name_bfs('B').distance_to_parent == 0.5
        assert tree.get_node_with_name_bfs('H').distance_to_parent == 0.4

    def test_get_leaves_in_subtree(self):
        tree = Tree("(((A:0.1,B:0.5)C:0.12,D:0.55,X:0.01)E:0.4,((F:0.8,G:0.7)Y:0.3)N:0.1,Z:0.2)H;")

        leaves = tree.get_leaves_in_subtree('E')
        names = [n.name for n in leaves]

        assert len(leaves) == 4
        assert 'A' in names
        assert 'B' in names
        assert 'D' in names
        assert 'X' in names

    def test_read_from_file(self):
        tree = Tree(file_name="aux_files/tree.nw")
        assert tree.root is not None

    def test_prune(self):
        tree = Tree(file_name="aux_files/tree_2.nw")

        to_remove = tree.prune()
        assert len(to_remove) > 0

    def test_leaves_to_remove(self):
        tree = Tree(file_name="aux_files/tree_2.nw")

        removed, newlist = tree.get_leaves_to_remove()

        assert len(removed) == 0
        assert len(newlist) == 3

    def test_leaves_to_remove_2(self):
        tree = Tree("(ref10_m0:0.0107059334,((((ref00_m0:0.0000010214,ref03_m0:0.0116010782)85.5/89:0.0234689348,ref05_m0:0.0093238664)94.4/98:0.0495507686,ref13_m0:0.0000026163)91.1/90:0.1176653259,((ref16_m2:0.0097760630,ref17_m2:0.0000010214)90.5/95:0.3007620051,((((ref02_m1:0.0403785457,ref03_m1:0.0495078624)87.9/84:0.0584232271,ref14_m1:0.0457792986)90.8/84:0.1432352643,ref15_m1:0.0000028335)98.3/99:0.6101056409,(ref01_m2:0.0000010214,ref03_m2:0.0000010214)91.1/98:0.3456814969)82.8/77:0.2446626484)86.7/76:0.2038092310)85.8/88:0.0780401341,ref09_m0:0.0000010214);")

        to_remove = tree.prune()
        removed, newlist = tree.get_leaves_to_remove()

    def test_get_supports(self):
        tree = Tree(file_name="aux_files/tree_2.nw")
        supports = tree.get_supports()

        assert len(supports) == 7
        assert 100 in supports
        assert 48 in supports