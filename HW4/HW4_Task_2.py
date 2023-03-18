# Homework 4, Task 2
# Напишите марковскую модель, которая имеет открытые состояния {A, T, G, C}, и скрытые состояния {+, -}. Когда модель
# в состоянии +, то вероятность генерации некоторого символа нуклеотида соответствует его частоте внутри CpG островков,
# вычислиному в первом задании, если состояние -, то частоте вне островков. Вероятность остаться внутри островка 0.95,
# а перейти в обычный геном 0.05. Для остальной части генома соответствующие вероятности 0.995 и 0.005. Саму модель
# можно реализовать в виде итератора, определив метод next, который возвращает пару - состояние и нуклеотид, который
# в этом состоянии произведен.

# Воспользуйтесь данной моделью для того чтобы сгенерировать набор из 20 последовательностей длинной от 1 000
# до 100 000, причем к каждой последовательности должна прилагаться последовательность состояний.

import random
import HW4_Task_1

x = ['A', 'C', 'G', 'T'] # open states
p = ['+', '-'] # hidden states

seq = [[None, None] for _ in range(20)]
for i in range(20):
    seq_len = random.randint(1000, 100000)
    hs = random.choice(p)
    os = random.choice(x)
    cur_seq = [os]
    cur_states = [hs]

    for j in range(seq_len - 1):
        for _ in range(4):
            if os == 'A':
                num_line = 1
            elif os == 'C':
                num_line = 2
            elif os == 'G':
                num_line = 3
            else:
                num_line = 4
        if hs == '+':
            nucleotide = random.choices(x, HW4_Task_1.weight_islands[num_line][1:5])
        else:
            nucleotide = random.choices(x, HW4_Task_1.weight_nonIslands[num_line][1:5])

        cur_seq.append(*nucleotide)
        cur_states.append(*hs)

        if hs == '+':
            hs = random.choices(p, [0.95, 0.05])
        else:
            hs = random.choices(['-', '+'], [0.995, 0.005])

    seq[i][0] = ''.join(cur_seq)
    seq[i][1] = ''.join(cur_states)



