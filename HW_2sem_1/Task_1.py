# 1. Симуляция секвенирования (3 балла)
# Первое задание заключается в том, чтобы написать симуляцию секвенирования ридов при помощи illumina. Для этого
# необходимо сгенерировать случайно строку длинной 50 000 символов, над алфавитом {A, T, G, C}. После этого случайно
# выбирать подстроки, длинна которых распределена нормально со средним 250 и среднеквадратическим отклонением 30,
# каждую такую подстроку "секвенировать".
# Под словом "секвенировать" в данном случае имеется ввиду симуляция того процесса, который мы обсуждали на лекции.
# Вы считываете очередной символ подстроки, например, 100 раз, но в N случаях из 100 ошибаетесь и считываете любой
# нуклеотид, кроме правильного с одинаковой вероятностью. N принадлежит равномерному распределению от 0 до 60.
# По получившемуся набору вычисляете наибоее вероятный нуклеотид (тот, которого больше всего) и его качество прочтения,
# чтобы записать его в формате FASTQ.
# При симуляции важно запоминать ID рида и для каждого ID позиции, в которых нуклеотид был считан неверно и какой
# должен быть на самом деле (для этого просто можно хранить 2 числа и помнить что у вас есть исходная строка,
# откуда берутся риды). Это понадобится для выполнения следующих заданий.
# Всего ридов пусть будет 50К.

import random
import numpy as np
from math import log10
alphabet = ["A", "T", "G", "C"]

def write_fastq(seq):
    with open('my_seq.fastq', 'w') as file:
        lines =[]
        for i in range(len(seq)):
            lines.append('@ Sequense #' + str(i+1) + '\n' + str(seq[i][0])+'\n'+'+'+'\n'+ str(seq[i][1]))
        file.write("\n".join(str(line) for line in lines))

def reading_seq(sequence):
    length = len(sequence)
    new_seq, quality, errors_ID = [], [], []
    reads = [[] for _ in range(100)]

    for i in range(length):
        count = [["A", 0], ["T", 0], ["G", 0], ["C", 0]]
        n = random.randrange(61)
        for k in range(4):
            if sequence[i] == count[k][0]:
                count[k][1] = 100 - n
                cur_nucleo_in_position = [count[k][0]]*count[k][1]
        n_alphabet = alphabet.copy()
        n_alphabet.remove(sequence[i])
        for j in range(n):
            nucleotide = random.choice(n_alphabet)
            cur_nucleo_in_position.append(nucleotide)
            for k in range(4):
                if nucleotide == count[k][0]:
                    count[k][1] += 1

        random.shuffle(cur_nucleo_in_position)
        for i in range(100):
            reads[i].append(cur_nucleo_in_position[i])
        final_nucleotide = max(count, key=lambda x: x[1])

        new_seq.append(final_nucleotide[0])
        quality.append(chr(int(-10*log10(final_nucleotide[1]/100)) + 33))

        if final_nucleotide[0] != sequence[i]:
            errors_ID.append(i)

    count.clear()

    for i in range(100):
        reads[i] = ''.join(reads[i])

    return ''.join(new_seq), ''.join(quality), errors_ID, reads



line=[]
for i in range(50000):
    line.append(random.choice(alphabet))

seq = ''.join(line)

ss_lengths = []
for i in range(100):
    ss_lengths.append(int(np.random.normal(250, 30)))

ss = []
for i in range(len(ss_lengths)):
    start = random.randrange(len(seq) - ss_lengths[i])
    ss.append([start, seq[start:start+ss_lengths[i]]])

nsg = []
all_reads = []
for i in range(len(ss)):
    sseq, q, err, readings = reading_seq(ss[i][1])
    nsg.append([sseq, q, err])
    all_reads.append(readings)

write_fastq(nsg)
