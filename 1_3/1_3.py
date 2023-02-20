import numpy as np
from binarytree import Node
from typing import NoReturn


class HTree:
    def __init__(self, deletions: float, insertions: float, mismatch: float, match: float):
        self.deletions = deletions
        self.insertions = insertions
        self.mismatch = mismatch
        self.match = match

    def tree(self, x: str, y: str) -> Node:
        # build tree
        base = Node(str((x, y)))
        self.hb(x, y, base)
        return base

    def hb(self, x: str, y: str, tree) -> NoReturn:
        # exit condition
        x_size, y_size = len(x), len(y)
        if (x_size <= 1 or
                y_size <= 1):
            return

        # initializing indexes
        middle = (x_size // 2)
        k = (self.align(x[:middle], y)
             + self.align(x[middle:][::-1], y[::-1])[::-1]).argmax()
        
        # right branch
        right_tuple = str((x[middle:], y[k:]))
        tree.right = Node(right_tuple)
        self.hb(x[middle:], y[k:], tree.right)

        #left branch
        left_tuple = str((x[:middle], y[:k]))
        tree.left = Node(left_tuple)
        self.hb(x[:middle], y[:k], tree.left)

    def align(self, x: str, y: str) -> np.ndarray:
        x, y = f'#{x}', f'#{y}'
        x_size, y_size = len(x), len(y)
        
        matrix_shape = (2, y_size)
        matrix = np.zeros(matrix_shape)
        for i in range(x_size):
            for j in range(y_size):
                matrix[1, j] = max(matrix[0, j - 1] + (self.match if x[i] == y[j] else self.mismatch),
                                   matrix[1, j - 1] + self.insertions,
                                   matrix[0, j] + self.deletions,
                                   0)
        return matrix[1, :]


if __name__ == "__main__":
    seq1 = 'AGTACGCA'
    seq2 = 'AGTACGCA'
    deletions, insertions = -3, -3
    mismatch, match = -2, 3
    ht = HTree(deletions, insertions, mismatch, match)
    print(ht.tree(seq1, seq2))
    
#output
# 1.
# seq1 = 'AGTACGCA'
# seq2 = 'AGTACGCA'
#                                             ___________________________________________('AGTACGCA', 'AGTACGCA')__________________________________________
#                                            /                                                                                                             \
#                   _________________('AGTA', 'AGTA')________________                                                               _________________('CGCA', 'CGCA')________________
#                  /                                                 \                                                             /                                                 \
#       _____('AG', 'AG')____                               _____('TA', 'TA')____                                       _____('CG', 'CG')____                               _____('CA', 'CA')____
#      /                     \                             /                     \                                     /                     \                             /                     \
# ('A', 'A')              ('G', 'G')                  ('T', 'T')              ('A', 'A')                          ('C', 'C')              ('G', 'G')                  ('C', 'C')              ('A', 'A')
# 
# 2.
# seq1 = 'AGTACGCA'
# seq2 = 'GTCGTCCT'
#                    __________________________________________('AGTACGCA', 'TATGC')__________________________________________
#                   /                                                                                                         \
#       _____('AGTA', 'TA')________________                                                            _________________('CGCA', 'TGC')_____
#      /                                   \                                                          /                                     \
# ('AG', '')                      _____('TA', 'TA')____                                    _____('CG', 'TG')____                        ('CA', 'C')
#                                /                     \                                  /                     \
#                           ('T', 'T')              ('A', 'A')                       ('C', 'T')              ('G', 'G')
#
# 3.
# seq1 = 'AGTACGCA'
# seq2 = 'TATGC'
#                     ___________________('AGTACGCA', 'GTCGTCCT')___________________
#                    /                                                              \
#        _____('AGTA', 'GT')_____                                       _____('CGCA', 'CGTCCT')__________________
#       /                        \                                     /                                         \
# ('AG', 'G')                ('TA', 'T')                          ('CG', '')                          _____('CA', 'CGTCCT')______
#                                                                                                    /                           \
#                                                                                               ('C', 'C')                  ('A', 'GTCCT')
