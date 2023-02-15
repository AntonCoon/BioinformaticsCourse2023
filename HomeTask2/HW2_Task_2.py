# Homework 2, Task 1
# Афинные гэпы (4 балла)
# Реализуйте выравнивание с афинными гэпами, алгоритм на вход принимает две строки, матрицу замен, штраф за начало гэпа
# α, и за его продолжение β. В результате возвращает выравнивание и его вес. Сложность алгоритма квадратичная по
# памяти и по времени.

import numpy as np
import Bio.Align as Align

matrix = Align.substitution_matrices.load("BLOSUM62")
#print(matrix)

def gotoh_algo(seq1: str, seq2: str, alpha: float, beta: float):
    len_seq1 = len(seq1)
    len_seq2 = len(seq2)

    F = np.zeros((len_seq1 + 1, len_seq2 + 1))
    A = np.zeros((len_seq1 + 1, len_seq2 + 1)) #delete matrix
    B = np.zeros((len_seq1 + 1, len_seq2 + 1)) #insert matrix

    for k in range(1, len_seq2 + 1):
        A[0][k] = alpha + beta * k
        B[0][k] = alpha + beta * k
        F[0][k] = alpha + beta * k
    for k in range(1, len_seq1 + 1):
        A[k][0] = alpha + beta * k
        B[k][0] = alpha + beta * k
        F[k][0] = alpha + beta * k

    for i in range(1, len_seq1 + 1):
        for j in range(1, len_seq2 + 1):
            A[i][j] = max(A[i-1][j] + beta, F[i-1][j] + alpha)
            B[i][j] = max(B[i][j-1] + beta, F[i][j-1] + alpha)
            Match = F[i - 1][j - 1] + matrix[seq1[i - 1]][seq2[j - 1]]
            F[i][j] = max(Match, A[i][j], B[i][j])

    score = F[len_seq1][len_seq2]

    i, j = len_seq1, len_seq2
    r_seq1, r_seq2 = [], []
    while i > 0 or j > 0:

        if i > 0 and F[i][j] == max(A[i][j], F[i-1][j] + alpha):
            r_seq1.append(seq1[i - 1])
            r_seq2.append('-')
            i -= 1
        elif j > 0 and F[i][j] == max(B[i][j], F[i][j-1] + alpha):
            r_seq1.append('-')
            r_seq2.append(seq2[j - 1])
            j -= 1
        else:
            r_seq1.append(seq1[i - 1])
            r_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
    # Reverse the strings.
    r_seq1 = ''.join(r_seq1)[::-1]
    r_seq2 = ''.join(r_seq2)[::-1]
    return score, '\n'.join([r_seq1, r_seq2])


seq1 = input('Enter sequence # 1: ')
seq2 = input('Enter sequence # 2: ')
gap_start_penalty = int(input('Enter the start gap penalty. Negative number! '))
gap_stop_penalty = int(input('Enter the stop gap penalty. Negative number! '))
s, align = gotoh_algo(seq1, seq2, gap_start_penalty, gap_stop_penalty)
print('Score:', s)
print(align)
