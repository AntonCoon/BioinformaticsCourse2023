# Homework 3, Task 1
# Вызовы Хиршберга (от 3 до 8 баллов)
# Реализуйте алгоритм, который принимал бы на вход две строки, штраф за удаления, вставки и несовпадения, а так-же цену
# совпадений. А возвращал бы дерево рекурсивных вызовов алгоритма Хоршберга. В вершинах дерева должны храниться пары
# строк, на которых будут происходит вызовы алгоритма Хиршберга.

import numpy as np
from binarytree import Node


class hirschberg_algo_points:

    def __init__(self, delete: float, insert: float, match: float, mismatch: float):
        self.delete = delete
        self.insert = insert
        self.match = match
        self.mismatch = mismatch

    def needleman_wunsch_algo_for_hb(self, seq1: str, seq2: str):
        seq1, seq2 = f'#{seq1}', f'#{seq2}'
        len_seq1, len_seq2 = len(seq1), len(seq2)

        F = np.zeros((2, len_seq2))

        for i in range(len_seq1):
            for j in range(len_seq2):
                Match = F[0][j - 1] + (self.match if seq1[i] == seq2[j] else self.mismatch)
                Delete = F[0][j] + self.delete
                Insert = F[1][j - 1] + self.insert
                F[1][j] = max(Match, Delete, Insert, 0)

        return F[1, :]

    def built_tree(self, seq1: str, seq2: str) -> Node:
        cur_node = Node(str((seq1, seq2)))
        self.hirschberg_algo(seq1, seq2, cur_node)
        return cur_node

    def hirschberg_algo(self, seq1: str, seq2: str, cur_tree):
        len_seq1, len_seq2 = len(seq1), len(seq2)

        if len_seq1 < 2 or len_seq2 < 2:
            return

        seq1_mid = len_seq1 // 2
        seq2_split_point = (self.needleman_wunsch_algo_for_hb(seq1[:seq1_mid], seq2)
                            + self.needleman_wunsch_algo_for_hb(seq1[seq1_mid:][::-1], seq2[::-1])[::-1]).argmax()

        cur_tree.right = Node(str((seq1[seq1_mid:], seq2[seq2_split_point:])))
        self.hirschberg_algo(seq1[seq1_mid:], seq2[seq2_split_point:], cur_tree.right)

        cur_tree.left = Node(str((seq1[:seq1_mid], seq2[:seq2_split_point])))
        self.hirschberg_algo(seq1[:seq1_mid], seq2[:seq2_split_point], cur_tree.left)


seq1 = 'ACGATATACA'  # input('Enter sequence # 1: ')
seq2 = 'CATATTA'  # input('Enter sequence # 2: ')
delete = -10
insert = -5
match = 5
mismatch = -7

align_conditions = hirschberg_algo_points(delete, insert, match, mismatch)
print(align_conditions.built_tree(seq1, seq2))
