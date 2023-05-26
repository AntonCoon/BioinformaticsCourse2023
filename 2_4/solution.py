import random
import sys
from typing import Set

import graphviz
from Bio import Phylo
from Bio.Phylo.Newick import Clade


class Node:
    def __init__(self, left=None, right=None, letter=None):
        if left is not None and right is not None:
            self.left: Node = left
            self.right: Node = right
            intersection = self.left.letters.intersection(self.right.letters)
            if len(intersection) == 0:
                self.letters = self.left.letters.union(self.right.letters)
            else:
                self.letters = intersection
        else:
            self.left: Node = None
            self.right: Node = None
            self.letters: Set[str] = set(letter)

    def get_dot(self):
        dot = graphviz.Digraph('phylo tree', comment='Phylogenetic tree')

        def _fill_dot(v, index=1):
            node_idx = f'V_{index}'
            label = "|".join(list(v.letters))
            dot.node(node_idx, label)
            if v.left is not None and v.right is not None:
                left_idx, right_idx = 2 * index, 2 * index + 1
                _fill_dot(v.left, left_idx)
                _fill_dot(v.right, right_idx)
                dot.edge(node_idx, f'V_{left_idx}')
                dot.edge(node_idx, f'V_{right_idx}')

        _fill_dot(self)

        return dot

    def view(self):
        self.get_dot().view()

    def save(self, filename):
        self.get_dot().render(filename)


def build_tree(node: Clade):
    if node.name is not None:
        return Node(letter=node.name)
    return Node(build_tree(node.clades[0]), build_tree(node.clades[1]))


def downhill(node: Node, letter=None):
    if letter is None or letter not in node.letters:
        letter = random.choice(list(node.letters))
    node.letters = set(letter)
    if node.left is not None and node.right is not None:
        downhill(node.left, letter)
        downhill(node.right, letter)


def main(input_file_name):
    bio_tree = Phylo.read(input_file_name, 'newick')
    tree = build_tree(bio_tree.root)
    downhill(tree)
    tree.view()


if __name__ == '__main__':
    main(sys.argv[1])
