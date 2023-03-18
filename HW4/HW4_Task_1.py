# Homework 4, Task 1
# Определите частоты генерации для каждого из нуклеотидов внутри CpG островков и вне их. Посчитайте так-же частоты
# для всех упорядоченных пар нуклеотидов и сравните частоту пары CG внутри островков и снаружи. Сделайте вывод.

from Parsing_fasta import read_fasta

def count_weights_single(ids, lines):
    weight = [['A', 0], ['C', 0], ['G', 0], ['T', 0]]
    for i in range(len(ids)):
        line = lines[i]
        len_line = len(lines[i])
        for j in range(len_line):
            for k in range(4):
                if line[j] == weight[k][0]:
                    weight[k][1] += 1

    total = sum(weight[i][1] for i in range(4))
    for k in range(4):
        weight[k][1] = round(weight[k][1] /total, 3)

    return weight

def count_weights(ids, lines):
    weight = [['-', 'A', 'C', 'G', 'T'], ['A', 0, 0, 0, 0], ['C', 0, 0, 0, 0], ['G', 0, 0, 0, 0], ['T', 0, 0, 0, 0]]
    for i in range(len(ids)):
        line = lines[i]
        len_line = len(lines[i]) - 1
        for j in range(len_line):
            pair = line[j] + line[j + 1]
            for y in range(1, 5):
                for z in range(1, 5):
                    if pair == weight[y][0] + weight[0][z]:
                        weight[y][z] += 1

    for i in range(1, 5):
        total = weight[i][1] + weight[i][2] + weight[i][3] + weight[i][4]
        for j in range(1, 5):
            weight[i][j] = round(weight[i][j] / total, 3)

    return weight

nonIslands_id, nonIslands_lines = read_fasta('nonIslands.fasta')
islands_id, islands_lines = read_fasta('islands.fasta')

# single
single_weight_nonIslands, single_weight_islands = count_weights_single(nonIslands_id, nonIslands_lines), \
                                                  count_weights_single(islands_id, islands_lines)

# pairs
weight_nonIslands, weight_islands = count_weights(nonIslands_id, nonIslands_lines), \
                                    count_weights(islands_id, islands_lines)
