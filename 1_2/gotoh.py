import numpy as np
from collections import deque

def gotoh(a: str, b: str, S: dict, alpha: float, beta: float) -> tuple[str, str, float]:
    """Gotoh Algorithm

    Parameters
    ----------
    a : str
        first seq
    b : str
        second seq
    S : dict
        substitution matrix
    alpha: float
        gap open penalty
    beta: float
        gap extension penalty

    Returns
    -------
    tuple of (aligned seq a, aligned seq b, alignment score)
    """
    # добавляем фиктивный символ в начало каждой строки для удобства дальнейшей индексации
    a = "-" + a 
    b = "-" + b
    A = np.zeros((len(a), len(b)))
    A[0, :] = - np.inf
    B = np.zeros((len(a), len(b)))
    B[:, 0] = - np.inf
    D = np.zeros((len(a), len(b)))
    D[:, 0] = alpha + beta*np.arange(len(a))
    D[0, :] = alpha + beta*np.arange(len(b))
    D[0, 0] = 0
    for i in range(1, len(a)):
        for j in range(1, len(b)):
            A[i, j] = max(A[i-1, j] + beta, D[i-1, j] + alpha + beta)
            B[i, j] = max(B[i, j-1] + beta, D[i, j-1] + alpha + beta)
            if (a[i], b[j]) in S:
                Sij = S[(a[i], b[j])]
            else:
                Sij = S[(b[j], a[i])]
            D[i, j] = max(A[i, j], B[i, j], D[i-1, j-1] + Sij)
    # восстанавливаем последовательности
    al_a = deque()
    al_b = deque()
    i = len(a) - 1
    j = len(b) - 1
    while i > 0 and j > 0:
        if D[i, j] == B[i, j]:
            al_a.appendleft('-')
            al_b.appendleft(b[j])
            j -= 1
        elif D[i, j] == A[i, j]:
            al_a.appendleft(a[i])
            al_b.appendleft('-')
            i -= 1
        else:
            al_a.appendleft(a[i])
            al_b.appendleft(b[j])
            i -= 1
            j -= 1
    while i > 0:
        al_a.appendleft(a[i])
        al_b.appendleft('-')
        i -= 1
    while j > 0:
        al_a.appendleft('-')
        al_b.appendleft(b[j])
        j -= 1
    return ''.join(al_a),''.join(al_b), D[-1, -1]


if __name__ == "__main__":
    from Bio.SubsMat import MatrixInfo as matlist
    matrix = matlist.blosum62
    print(gotoh("GATTACA", "CATGA", matrix, -100, 0))