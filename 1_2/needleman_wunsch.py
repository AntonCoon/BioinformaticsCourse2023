import numpy as np
from collections import deque

def needleman_wunsch(a: str, b: str, S: dict, d: float):
    """Needleman-Wunsch Algorithm

    Parameters
    ----------
    a : str
        first seq
    b : str
        second seq
    S : dict
        substitution matrix
    d: float
        gap penalty

    Returns
    -------
    tuple of (aligned seq a, aligned seq b, alignment score)
    """
    # добавляем фиктивный символ в начало каждой строки для удобства дальнейшей индексации
    a = "-" + a 
    b = "-" + b
    f = np.zeros((len(a), len(b)))
    for i in range(len(a)):
        for j in range(len(b)):
            if min(i, j) == 0:
                f[i, j] = d*(i + j)
            else:
                delete = f[i-1, j] + d
                insert = f[i, j-1] + d
                try:
                    Sij = S[(a[i], b[j])]
                except KeyError:
                    Sij = S[(b[j], a[i])]
                f[i, j] = max(delete, insert, f[i-1, j-1] + Sij)
    # восстанавливаем последовательности
    al_a = deque()
    al_b = deque()
    i = len(a) - 1
    j = len(b) - 1
    while i > 0 and j > 0:
        try:
            Sij = S[(a[i], b[j])]
        except KeyError:
            Sij = S[(b[j], a[i])]
        if f[i, j] == f[i-1, j-1] + Sij:
            al_a.appendleft(a[i])
            al_b.appendleft(b[j])
            i -= 1
            j -= 1
        elif f[i, j] == f[i, j-1] + d:
            al_a.appendleft('-')
            al_b.appendleft(b[j])
            j -= 1
        else:
            al_a.appendleft(a[i])
            al_b.appendleft('-')
            i -= 1
    while i > 0:
        al_a.appendleft(a[i])
        al_b.appendleft('-')
        i -= 1
    while j > 0:
        al_a.appendleft('-')
        al_b.appendleft(b[j])
        j -= 1
    return ''.join(al_a),''.join(al_b), f[-1, -1]


if __name__ == "__main__":
    from Bio.SubsMat import MatrixInfo as matlist
    matrix = matlist.blosum62
    print(needleman_wunsch("GATTACA", "CATTAGA", matrix, -2))