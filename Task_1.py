# Реализуйте алгоритм фитча, который принимает на вход дерево в формате newick, в листьях этого дерева записаны символы
# из генома организма на некоторой позиции. Обратите внимание, что парсинг этого формата делать не нужно, можно
# пользоваться готовым from Bio import Phylo. Результат выполнения алгоритма - дерево, с заполненными внутренними
# узлами согласно принципу максимальной парсимонии.
# Пример:
# вход: (((A, A), C), (C, G))
# выход: дерево

from Bio import Phylo
from typing import Any, Dict
from random import choice
import io

def built_phylotree(tree: Phylo.BaseTree):
    for clade in tree.get_terminals():
        clade.name = clade.name.upper()

    for clade in tree.find_clades(order='postorder'):
        if not clade.is_terminal():
            counts: Dict[Any, int] = {}
            for child in clade.clades:
                counts[child.name] = counts.get(child.name, 0) + 1
            max_count = max(counts.values())
            names = [name for name, count in counts.items() if count == max_count]
            clade.name = choice(names)

    return tree

treedata = "(((A, A), C), (C, G))"
#treedata = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5)"
tree = Phylo.read(io.StringIO(treedata), "newick")
#print(tree)
phylotree = built_phylotree(tree)
Phylo.draw(phylotree)