import numpy as np
from typing import Union
import copy


class Graph:
    def __init__(self, seqs: list[str], threshhold: int) -> None:
        self.seqs = seqs
        self.total = len(seqs)
        self.threshold = threshhold
        self.g = {i: [] for i in range(self.total)}
        self.build_graph()

    def build_graph(self) -> None:
        # строим граф перекрытий с трешхолдом
        for i in range(self.total):
            for j in range(self.total):
                w = longest_common_suf_pref(self.seqs[i], self.seqs[j])
                if i != j and w >= self.threshold:
                    self.g[i].append([j, w])

    def merge_nodes(self, node1: int, node2: int) -> None:
        # вводим новую вершину
        new_node = len(self.seqs)
        # из new_node исходят ребра которые исходили из node2
        self.g[new_node] = copy.deepcopy(self.g[node2])
        # ребра входившие в node1 теперь входят в new_node. Веса остаются те же
        for u in self.g:
            for i in range(len(self.g[u])):
                if self.g[u][i][0] == node1 and u != new_node:
                    self.g[u][i][0] = new_node
        # добавляем новую вершину в список вершин и удаляем старые
        self.seqs.append(merge_seqs(self.seqs[node1], self.seqs[node2]))
        del self.g[node1]
        del self.g[node2]

    def find_heaviest_edge(self) -> tuple[Union[int, None], Union[int, None]]:
        # находим самое тяжелое ребро в графе
        # если ребер в графе нет возвр None, None
        w = - float('inf')
        s, e = 0, 0
        for u in self.g:
            for i in range(len(self.g[u])):
                if self.g[u][i][1] >= w and self.g[u][i][0] in self.g:
                    w = self.g[u][i][1]
                    s = u
                    e = self.g[u][i][0]
        if w == - float('inf'):
            return None, None
        else:
            return s, e
        
    def find_scs(self):
        while len(self.g) > 1:
            # находим самое тяжелое ребро
            s, e = self.find_heaviest_edge()
            if s == None and e == None:
                # есди нет ребер выбираем случайно две вершины
                s, e = np.random.choice(self.g, 2, replace=False)
            # межджим вершины
            self.merge_nodes(s, e)
        return self.seqs[-1]


def merge_seqs(s1: str, s2: str) -> str:
    # мерджим строки s1 и s2
    for i in range(min(len(s1), len(s2))):
        if s1[len(s1) - i - 1] == s2[i]:
            continue
        else:
            break
    return s1 + s2[i:]

def longest_common_suf_pref(s1: str, s2: str) -> int:
    # находим вес ребра между s1 и s2
    for i in range(min(len(s1), len(s2))):
        if s1[len(s1) - i - 1] == s2[i]:
            continue
        else:
            break
    return i


if __name__ == '__main__':
    g = Graph(['AAA', 'AAT', 'ATT', 'TTA', 'TTT'], 1)
    print(g.find_scs())