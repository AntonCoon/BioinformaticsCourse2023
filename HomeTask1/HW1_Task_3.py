# Homework 1, Task 3
# Левенштейн (4 балла)
# Реализуйте алгоритм вычисления расстояния редактирования. O(n•m) по времени и O(min(n, m)) по памяти (обратите
# внимание на то, что возвращать способ получения из одной строки другую не нужно). Алгоритм возвращаеет одно
# число - расстояние редактирования между входными строками.

def Edit_distance(line_1, line_2):
    lenght_1, lenght_2 = len(line_1), len(line_2)
    
    if lenght_1 > lenght_2:
        line_1, line_2 = line_2, line_1
        lenght_1, lenght_2 = lenght_2, lenght_1

    # создаем массив для запоминания строк после преобразования
    rmmbr = [[0 for i in range(lenght_1 + 1)] for j in range(2)]

    # первая строка - входные условия: вторая строка пустая
    for i in range(lenght_1 + 1):
        rmmbr[0][i] = i

    # само заполнение
    # цикл по каждому символу второй строки
    for i in range(1, lenght_2 + 1):
        # цикл для сравнения символа второй строки с символом первой
        for j in range(0, lenght_1 + 1):
            if j == 0:
                rmmbr[i % 2][j] = i
            elif (line_1[j - 1] == line_2[i - 1]):
                rmmbr[i % 2][j] = rmmbr[(i - 1) % 2][j - 1]
            else:
                rmmbr[i % 2][j] = (1 + min(rmmbr[(i - 1) % 2][j], min(rmmbr[i % 2][j - 1], rmmbr[(i - 1) % 2][j - 1])))

    return rmmbr[lenght_2 % 2][lenght_1]
