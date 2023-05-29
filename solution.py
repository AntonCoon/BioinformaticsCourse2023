import random
from Bio import Phylo
from Bio.Phylo.Newick import Clade


class Node:
    def __init__(self, left = None, right = None, word = None):
        self.left = left
        self.right = right
        self.words = set()
        
        if self.left is not None and self.right is not None:
            intersection = self.left.words.intersection(self.right.words)
            if len(intersection) == 0:
                self.words = self.left.words.union(self.right.words)
            else:
                self.words = intersection
        elif word is not None:
            self.words = set(word)


def tree_build(node: Clade):
    if node.name is not None:
        return Node(word = node.name)
    return Node(tree_build(node.clades[0]), tree_build(node.clades[1]))


def tree_read(node: Node, word = None):
    if word is not None and word in node.words:
        node.words = set(word)
    else:
        node.words = set(random.choice(list(node.words)))
        
    if node.left is not None and node.right is not None:
        tree_read(node.left, word)
        tree_read(node.right, word)



file_path_name = 'TEST_tree_newick.ph'
tree_input = Phylo.read(file_path_name, 'newick')
tree_result = tree_build(tree_input)
tree_read(tree_result)
Phylo.draw(tree_result)
