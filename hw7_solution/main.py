import random
import numpy as np
import matplotlib.pyplot as plt


# make sequence and flip it k times
def simulate(k: int):
    sequence = np.arange(1, 1001)
    
    for _ in range(k):
        index1 = random.randint(0, len(sequence) - 1)
        index2 = random.randint(index1, len(sequence) - 1)
        reverse_segment(index1, index2+1, sequence)
        
    return np.array(sequence)

# reverse just a segment of a sequence
def reverse_segment(a:int, b:int, sequence):
    sequence[a:b] = sequence[a:b][::-1]
    for i in range(a, b):
        sequence[i] = - sequence[i]

# greedy flip sort algorithm
def flip_sort(sequence):
    sequence = sequence.copy()
    n = len(sequence)
    counter = 0
    for i in range(n):
        if i+1 == sequence[i]:
            continue
        
        tmp = i
        for j in range(i, n):
            if abs(sequence[j]) == i+1:
                tmp = j
                break
        reverse_segment(i, tmp+1, sequence)
        counter += 1
        
        if i+1 == - sequence[i]:
            reverse_segment(i, i+1, sequence)
            counter += 1
            
    return counter


f = 10  # number of flips
result = simulate(f)
re = flip_sort(result)
print(re)  # to watch the amount of actions greedy algorithm needs to make sequence great again

k_used = []
algo_ans = []
n = 1000
for k in range(1, n*2, 20):
    k_used.append(k/n)
    res = simulate(k)
    ans = flip_sort(res)
    algo_ans.append(ans/n)
    
# If you want to see the graph, just run this code 
plt.plot(k_used, algo_ans)
plt.plot(k_used, k_used)
plt.xlabel("Given reverses quantity")
plt.ylabel("Used reverses quantity")
plt.legend(["What we really see", "What we want to see"])
plt.show()

'''
Что видно из графика и в целом из результатов для различных k:

Изначально у нас есть некоторая последовательность {1,2,...,1000}, в которой мы абсолютно рандомным образом 
совершаем k поворотов ее сегментов. А затем пользуемся жадным алгоритмом для того, чтобы вернуть
последовательность в изначальное упорядоченное состояние. И мы можем наблюдать такой факт, что 
с ростом числа поворотов k алгоритм выходит на некое "плато" и, с течением итераций, в какой-то
момент количество действий алгоритма для возвращения последовательности к исходному состоянию 
становится меньше, чем реальное k. Так как в алгоритме присутствует некая доля случайности, то
определить точку поломки совсем точно не получится, однако, можно с высокой долей уверенности сказать,
что она располагается на отрезке [1300, 1500].
Соответственно, можем сказать, что при k > 1.5n жадный алгоритм будет делать меньше операций, чем 
нам потребовалось изначально.  

Если попробовать объяснить это аналитически, то можно заметить, что надо рассмотреть два варианта:
1. Когда нам надо сделать только локальный разворот, т.е. поменять знак у числа на конкретной позиции. 
   С вероятностью 1/2 нам нужен (не нужен) этот поворот, поэтому, если посчитать грубо, то 
   в среднем это 1/2 * n операций. 
2. Еще ~n операций (в худшем случае) нам надо сделать, чтобы повернуть все n элементов на свое место. 
Итого это ~1.5n операций, что отвечает эксперементальным результатам 
'''
