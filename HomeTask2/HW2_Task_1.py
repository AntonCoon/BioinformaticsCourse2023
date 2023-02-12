# Homework 2, Task 1
# Нидлман Вунш (3 балла)
# Реализуйте алгоритм Нидлмана Вунша для выравнивания последовательностеей. На вход принимается две строки, матрица
# замен и стоимость гэпа. В результате верните оптимальное выравнивание и его вес.
# При проверке помните, что оптимальных выравниваний может быть несколько, но вес у них должен совпадать.

import numpy as np
import Bio.Align as Align

matrix = Align.substitution_matrices.load("BLOSUM62")
#print(matrix)

def needleman_wunsch_algo(seq1: str, seq2: str, gap_penalty: int):
    len_seq1 = len(seq1)
    len_seq2 = len(seq2)

    F = np.zeros((len_seq1 + 1, len_seq2 + 1))
    F[:, 0] = np.linspace(0, len_seq1 * gap_penalty, len_seq1 + 1)
    F[0, :] = np.linspace(0, len_seq2 * gap_penalty, len_seq2 + 1)

    for i in range(1, len_seq1 + 1):
        for j in range(1, len_seq2 + 1):
            Match = F[i - 1][j - 1] + matrix[seq1[i - 1]][seq2[j - 1]]
            Delete = F[i - 1][j] + gap_penalty
            Insert = F[i][j - 1] + gap_penalty
            F[i][j] = max(Match, Delete, Insert)

    score = F[len_seq1][len_seq2]

    i, j = len_seq1, len_seq2
    r_seq1, r_seq2 = [], []
    while i > 0 or j > 0:

        if i > 0 and j > 0 and F[i][j] == F[i - 1][j - 1] + matrix[seq1[i - 1]][seq2[j - 1]]:
            r_seq1.append(seq1[i - 1])
            r_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and (F[i][j] == F[i - 1][j] + gap_penalty):
            r_seq1.append(seq1[i - 1])
            r_seq2.append('-')
            i -= 1
        else:
            r_seq1.append('-')
            r_seq2.append(seq2[j - 1])
            j -= 1
    # Reverse the strings.
    r_seq1 = ''.join(r_seq1)[::-1]
    r_seq2 = ''.join(r_seq2)[::-1]
    return score, '\n'.join([r_seq1, r_seq2])


seq1 = input('Enter sequence # 1: ')
seq2 = input('Enter sequence # 2: ')
gap = int(input('Enter the gap penalty. Negative number! '))
s, align = needleman_wunsch_algo(seq1, seq2, gap)
print('Score:', s)
print(align)
