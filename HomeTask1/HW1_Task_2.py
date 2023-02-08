# Homework 1, Task 2
# Ближайшая подстрока (1 балл)
# Даны 2 строки над алфавитом ATGC, реализуйте алгоритм, который позволит находить подстроку в более длинной строке,
# которая ближе всего по расстоянию Хэминга к боллее короткой. Верните позицию начала этой подстроки, саму подстроку
# и рассстояние Xэмминга.

from HW1_Task_1 import Hamming_distance

def nearest_substring(line_1: str, line_2: str):
    lenght_1, lenght_2 = len(line_1), len(line_2)
    if lenght_1 < lenght_2:
        line_1, line_2 = line_2, line_1
        lenght_1, lenght_2 = lenght_2, lenght_1

    all_substrings_list = [line_1[i: i+ lenght_2] for i in range(0, lenght_1 - lenght_2 + 1)]

    list_of_distances = [(i, substring, Hamming_distance(line_2, substring)) for i, substring in enumerate(all_substrings_list)]
    list_of_distances = sorted(list_of_distances, key=lambda x: (x[2], x[0]))

    return list_of_distances[0]