# Homework 1, Task 1
# Хэмминг (1 балл)
# Реализуйте вычисление расстояния Хэмминга.

def Hamming_distance(line_1: str, line_2: str):
    Ham_distance = 0
    if len(line_1) != len(line_2):
        raise ValueError('Different length of the the lines! Check data')
    for character_1, character_2 in zip(line_1, line_2):
        if character_1 != character_2:
            Ham_distance += 1
    return Ham_distance