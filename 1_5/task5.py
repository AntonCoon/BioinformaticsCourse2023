import numpy as np
from Bio import pairwise2


def multiple_alignment(strings):
    # инициализируем матрицу расстояний Левенштейна
    distance_matrix = np.zeros((len(strings), len(strings)))

    # заполняем матрицу расстояний Левенштейна между всеми парами строк
    for i in range(len(strings)):
        for j in range(i + 1, len(strings)):
            distance_matrix[i, j] = pairwise2.align.globalxx(strings[i], strings[j], score_only=True)
            distance_matrix[j, i] = distance_matrix[i, j]

    # пока не останется только одна строка
    while len(strings) > 1:
        # находим две строки с минимальным расстоянием Левенштейна
        min_distance = np.min(distance_matrix[np.nonzero(distance_matrix)])
        min_indices = np.where(distance_matrix == min_distance)
        i, j = min_indices[0][0], min_indices[1][0]

        # заменяем их консенсусной строкой
        alignment = pairwise2.align.globalxx(strings[i], strings[j], one_alignment_only=True)[0]
        consensus = ""
        for k in range(len(alignment[0])):
            if alignment[0][k] == alignment[1][k]:
                consensus += alignment[0][k]
            else:
                consensus += "-"

        # выводим результаты парного выравнивания
        print("Alignment {}: ".format(len(strings) - 2))
        print(strings[i])
        print(strings[j])
        print()

        # удаляем исходные строки и добавляем консенсусную в список строк
        strings.pop(j)
        strings.pop(i)
        strings.append(consensus)

        # пересчитываем матрицу расстояний Левенштейна
        distance_matrix = np.zeros((len(strings), len(strings)))
        for i in range(len(strings)):
            for j in range(i + 1, len(strings)):
                distance_matrix[i, j] = pairwise2.align.globalxx(strings[i], strings[j], score_only=True)
                distance_matrix[j, i] = distance_matrix[i, j]

    # выводим последнее выравнивание
    print("Alignment {}: ".format(len(strings) - 1))
    print(strings[0])
    print()

    # возвращаем конечную строку выравнивания
    return strings[0]
strings = ["ACGGACGTA", "AGGACGTG", "AGGACGTAG"]
alignment = multiple_alignment(strings)
print("Final alignment: ", alignment)
