import numpy as np 

alphabet = ['A', 'T', 'G', 'C']
sequence = np.random.choice(alphabet, 100)
print(sequence)
# Задаем длину подстроки
subseq_len = 15

# Открываем файл для записи
with open('D:/PYTHON/ALGORITHMS_IN_BIOINFORMATICS/10HW/reads.fastq', 'w') as f:
    # Проходим по строке с шагом subseq_len
    for i in range(0, len(sequence)-subseq_len+1, 1):
        # Выбираем подстроку длины subseq_len
        subseq = sequence[i:i+subseq_len]
        # Формируем запись в формате FASTQ
        record = '@read{}\n{}\n+\n{}\n'.format(i, ''.join(subseq), '~'*subseq_len)
        # Записываем запись в файл
        f.write(record  + '\n')
