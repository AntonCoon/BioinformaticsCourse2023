## Домашнее задание по первой лекции.
### Практическая часть
'''
Во всех практических заданиях ниже входные данные
будут храниться в формате [fasta](https://en.wikipedia.org/wiki/FASTA_format).
Писать парсер fasta файлов самостоятельно не рекомендуется,
пользуйтесь библиотеками. Если пишете на python, то можете
использовать, например, класс SeqIO из BioPython.

from Bio import SeqIO

Остальные функции в академических целях нужно
написать не используя библиотек. В следующих ДЗ или
для тестирования Вашего решения в этом, можно
использовать библиотеку [Distance](https://pypi.org/project/Distance/0.1/)
для python или аналогичные для других ЯП. Входные файлы,
на которых можно попробовать запустить свое решение
находятся в папке `./data`.
1. Хэмминг (1 балл)
Реализуйте вычисление расстояния Хэмминга.
2. Ближайшая подстрока (1 балл)
Даны 2 строки над алфавитом ATGC, реализуйте алгоритм,
 который позволит находить подстроку в более длинной
 строке, которая ближе всего по расстоянию Хэминга к
 боллее короткой. Верните позицию начала этой подстроки,
 саму подстроку и рассстояние хэмминга.
3. Левенштейн (4 балла)
Реализуйте алгоритм вычисления расстояния
редактирования. *O(n•m)* по времени и *O(min(n, m))* по памяти
(обратите внимание на то, что возвращать способ получения из
одной строки другую не нужно). Алгоритм возвращаеет одно число -
расстояние редактирования между входными строками.
'''

from Bio import SeqIO


def hamming_distance(string1: str, string2: str):
    if len(string1) != len(string2):
        raise Exception('Hamming distance can only be calculated for strings of equal length')
    else:
        distance = 0
        strings_length = len(string1)

        for i in range(strings_length):
            if string1[i] != string2[i]:
                distance += 1

        return distance


def nearest_substring(string1: str, string2: str):
    l1, l2 = len(string1), len(string2)

    if l1 > l2:
        string1, string2 = string2, string1
        l1, l2 = l2, l1

    substrings = [string2[i: i + l1] for i in range(0, l2 - l1 + 1)]
    return min((i, substring, hamming_distance(string1, substring)) for i, substring in enumerate(substrings))


def levenshtein_distance(string1: str, string2: str):
    l1, l2 = len(string1), len(string2)

    if l1 > l2:
        string1, string2 = string2, string1
        l1, l2 = l2, l1

    current_row = range(l1 + 1)
    for i in range(1, l2 + 1):
        previous_row, current_row = current_row, [i] + [0] * l1
        for j in range(1, l1 + 1):
            add, delete, change = previous_row[j] + 1, current_row[j - 1] + 1, previous_row[j - 1]
            if string1[j - 1] != string2[i - 1]:
                change += 1
            current_row[j] = min(add, delete, change)

    return current_row[l1]


if __name__ == '__main__':
    print("Проверка на первых тестовых данных")
    data = list(SeqIO.parse('f8.fasta', 'fasta'))
    s1 = str(data[0].seq)
    s2 = str(data[1].seq)
    print("Nearest substring:", nearest_substring(s1, s2))
    print("Levenshtein distance:", levenshtein_distance(s1, s2))

    print("Проверка на вторых тестовых данных")
    data = list(SeqIO.parse('gattaca.fasta', 'fasta'))
    s1 = str(data[0].seq)
    s2 = str(data[1].seq)
    print("Hamming distance:", hamming_distance(s1, s2))
    print("Nearest substring:", nearest_substring(s1, s2))
    print("Levenshtein distance:", levenshtein_distance(s1, s2))
