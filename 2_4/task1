from io import StringIO
from Bio import Phylo
from typing import Any, Dict
import random


def fill_internal_nodes(tree: Phylo.BaseTree) -> Phylo.BaseTree:
    for clade in tree.get_terminals():
        clade.name = clade.name.upper()

    for clade in tree.find_clades(order='postorder'):
        if not clade.is_terminal():
            counts: Dict[Any, int] = {}
            for child in clade.clades:
                counts[child.name] = counts.get(child.name, 0) + 1
            max_count = max(counts.values())
            most_common_names = [name for name, count in counts.items() if count == max_count]
            clade.name = random.choice(most_common_names)

    return tree


# Пример использования
newick_tree = "(((A, A), C), (C, G))"
tree = Phylo.read(StringIO(newick_tree), "newick")
filled_tree = fill_internal_nodes(tree)
Phylo.draw(filled_tree)
