import numpy as np
import binarytree #pip install binarytree

# вспомогательная функция
def nw(a: str, b: str, del_pen: float, ins_pen: float, mismatch_pen: float, match: float) -> np.ndarray:
    """Needleman-Wunsch Algorithm

    Parameters
    ----------
    a : str
        first seq
    b : str
        second seq
    del_pen: float
        penalty for deletion
    ins_pen: float
        penalty for insertion
    mismatch_pen: float
        penalty for mismatch
    match: float
        matching score

    Returns
    -------
    last row in pairwise sequence alignment table
    """
    a = "-" + a 
    b = "-" + b
    f = np.zeros((2, len(b)))
    for i in range(len(a)):
        for j in range(len(b)):
            if i == 0:
                f[1, j] = ins_pen*j
            elif j == 0:
                f[1, 0] = del_pen*i
            else:
                delete = f[0, j] + del_pen
                insert = f[1, j-1] + ins_pen
                f[1, j] = max(delete, insert, f[0, j-1] + (mismatch_pen if a[i]!=b[j] else match))
        f[0, :] = f[1, :]
    return f[1, :]


def hirschberg_tree(a: str, b: str, del_pen: float, ins_pen: float, mismatch_pen: float, match: float) -> binarytree.Node:
    """Возвращает дерево рекурсивных вызовов алгоритма Хиршберга.
    В вершинах дерева хранятся строки, 
    на которых происходят вызовы алгоритма Хиршберга.

    Parameters
    ----------
    a : str
        first seq
    b : str
        second seq
    del_pen: float
        penalty for deletion
    ins_pen: float
        penalty for insertion
    mismatch_pen: float
        penalty for mismatch
    match: float
        matching score

    Returns
    -------
    binarytree.Node object
    """
    def hirschberg(tree, a, b):
        if len(a) <= 1 or len(b) <= 1:
            return
        mid = len(a)//2
        scores = nw(a[:mid], b, del_pen, ins_pen, mismatch_pen, match) + \
            nw(a[mid:][::-1], b[::-1], del_pen, ins_pen, mismatch_pen, match)[::-1]
        cut_ind = np.argmax(scores)
        tree.left = binarytree.Node(str((a[:mid], b[:cut_ind])))
        tree.right = binarytree.Node(str((a[mid:], b[cut_ind:])))
        hirschberg(tree.left, a[:mid], b[:cut_ind])
        hirschberg(tree.right, a[mid:], b[cut_ind:])
    root = binarytree.Node(str((a, b)))
    hirschberg(root, a, b)
    return root


if __name__ == "__main__":
    a = 'AGTACGCA'
    b = 'TATGC'
    del_pen = -2
    ins_pen = -2
    match = 2
    mismatch = -1
    print(hirschberg_tree(a, b, del_pen, ins_pen, mismatch, match))