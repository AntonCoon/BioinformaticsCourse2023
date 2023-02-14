# Нидлман Вунш (3 балла)
# Реализуйте алгоритм Нидлмана Вунша для выравнивания последовательностеей. На вход принимается две строки,
# матрица замен и стоимость гэпа.
# В результате верните оптимальное выравнивание и его вес. При проверке помните, что оптимальных выравниваний может быть
# несколько, но вес у них должен совпадать.
import Bio.SubsMat.MatrixInfo as matlist
from Bio import SeqIO
import numpy as np
import math

matrix = matlist.blosum62


def NildmanWood(s1, s2, matrix, gap_weight):
    n = len(s1)
    m = len(s2)
    D = np.zeros((m + 1, n + 1))

    for i in range(m + 1):
        D[i][0] = gap_weight * i

    for j in range(n + 1):
        D[0][j] = gap_weight * j

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if (s1[j - 1], s2[i - 1]) in matrix:
                change = D[i - 1][j - 1] + matrix[(s1[j - 1], s2[i - 1])]

            else:
                change = D[i - 1][j - 1] + matrix[tuple(reversed((s1[j - 1], s2[i - 1])))]

            delete = D[i][j] + gap_weight
            insert = D[i][j] + gap_weight
            D[i][j] = max(change, insert, delete)
    result1 = ""
    result2 = ""
    i = m
    j = n
    # Смотрим откуда пришли
    while i > 0 and j > 0:
        if (s1[j - 1], s2[i - 1]) in matrix:
            if D[i][j] == D[i - 1][j - 1] + matrix[(s1[j - 1], s2[i - 1])]:
                result1 = result1 + s1[j - 1]
                result2 = result2 + s2[i - 1]
                i -= 1
                j -= 1
        else:
            if D[i][j] == D[i - 1][j - 1] + matrix[tuple(reversed((s1[j - 1], s2[i - 1])))]:
                result1 = result1 + s1[j - 1]
                result2 = result2 + s2[i - 1]
                i -= 1
                j -= 1
        if D[i][j] == D[i - 1][j] + gap_weight:
            result1 = result1 + s1[j - 1]
            result2 = result2 + "_"
            j -= 1
        elif D[i][j] == D[i][j - 1] + gap_weight:
            result1 = result1 + "_"
            result2 = result2 + s2[i - 1]
            i -= 1

    # доходим строки до конца
    while j > 0:
        result1 = result1 + s1[j - 1]
        result2 = result2 + "_"
        j -= 1
    while i > 0:
        result1 = result1 + "_"
        result2 = result2 + s2[i - 1]
        i -= 1

    result1 = result1[::-1]
    result2 = result2[::-1]

    return result1, result2, D[m][n]


str1, str2, weight = NildmanWood("AGTA", "ATA", matrix, -1)
print(f'Вес: {weight}')
print("Выравнивание: ")
print(str1)
print(str2)


# Афинные гэпы (4 балла)
# Реализуйте выравнивание с афинными гэпами, алгоритм на вход принимает две строки, матрицу замен, штраф за начало гэпа α,
# и за его продолжение β. В результате возвращает выравнивание и его вес. Сложность алгоритма квадратичная по памяти и по времени
def g(alpha, beta, k):
    return alpha + beta * k


def GotohAlgorithm(s1, s2, matrix, alpha, beta):
    n = len(s1)
    m = len(s2)
    D = np.zeros((m + 1, n + 1))
    for i in range(m + 1):
        D[i][0] = g(alpha, beta, i) * i

    for j in range(n + 1):
        D[0][j] = g(alpha, beta, j) * j

    A = np.zeros((m + 1, n + 1))
    for j in range(n + 1):
        A[0][j] = - math.inf

    B = np.zeros((m + 1, n + 1))
    for i in range(m + 1):
        B[i][0] = - math.inf

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            A[i][j] = max(A[i - 1][j] + beta, D[i - 1][j] + g(alpha, beta, 1))
            B[i][j] = max(B[i][j - 1] + beta, D[i][j - 1] + g(alpha, beta, 1))

            if (s1[j - 1], s2[i - 1]) in matrix:
                change = D[i - 1][j - 1] + matrix[(s1[j - 1], s2[i - 1])]

            else:
                change = D[i - 1][j - 1] + matrix[tuple(reversed((s1[j - 1], s2[i - 1])))]

            delete = A[i][j]
            insert = B[i][j]
            D[i][j] = max(change, insert, delete)

    result1 = ""
    result2 = ""
    i = m
    j = n

    # Смотрим откуда пришли
    while i > 0 and j > 0:
        if (s1[j - 1], s2[i - 1]) in matrix:
            if D[i][j] == D[i - 1][j - 1] + matrix[(s1[j - 1], s2[i - 1])]:
                result1 = result1 + s1[j - 1]
                result2 = result2 + s2[i - 1]
                i -= 1
                j -= 1
        else:
            if D[i][j] == D[i - 1][j - 1] + matrix[tuple(reversed((s1[j - 1], s2[i - 1])))]:
                result1 = result1 + s1[j - 1]
                result2 = result2 + s2[i - 1]
                i -= 1
                j -= 1
        if D[i][j] == A[i][j]:
            result1 = result1 + s1[j - 1]
            result2 = result2 + "_"
            j -= 1
        elif D[i][j] == B[i][j]:
            result1 = result1 + "_"
            result2 = result2 + s2[i - 1]
            i -= 1

    # доходим строки до конца
    while j > 0:
        result1 = result1 + s1[j - 1]
        result2 = result2 + "_"
        j -= 1
    while i > 0:
        result1 = result1 + "_"
        result2 = result2 + s2[i - 1]
        i -= 1

    result1 = result1[::-1]
    result2 = result2[::-1]

    return result1, result2, D[m][n]


str1, str2, weight = GotohAlgorithm("AGTA", "ATA", matrix, -3, -2)
print(f'Вес: {weight}')
print("Выравнивание: ")
print(str1)
print(str2)

