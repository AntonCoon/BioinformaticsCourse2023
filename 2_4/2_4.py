from Bio import Phylo
import matplotlib.pyplot as plt

import random

def fitch(tree):
    nodes = tree.find_clades(order='postorder')
    for node in nodes:
        if not node.is_terminal():
            children = node.clades
            intersection = set.intersection(*[set(child.name) for child in children])
            if intersection:
                node.name = ''.join(intersection)
            else:
                union = set.union(*[set(child.name) for child in children])
                node.name = random.choice(list(union))
    return tree


tree = Phylo.read('tree.newick', 'newick')
tree = fitch(tree)
Phylo.write(tree, 'output.newick', 'newick')
#Phylo.draw_ascii(tree)

Phylo.draw(tree, label_func=lambda x: x.name if x.name else '')
plt.show()