import random
import numpy as np
import itertools
import subprocess

# generate sequence over the alphabet
alphabet = ["A", "T", "G", "C"]
sequence = ''.join([random.choice(alphabet) for _ in range(50000)])

# pick a substring with normally distributed length 
lengths = np.random.normal(loc=250, scale=30, size=50000).astype(int)
subsequences = []
for length in lengths:
    start = random.randint(0, 50000 - length)
    end = start + length
    subsequence = sequence[start:end]
    subsequences.append(subsequence)
    
# sequencing 
reads = []
for i, subsequence in enumerate(subsequences):
    errors = np.random.uniform(low=0, high=60, size=len(subsequence)).astype(int)
    read = ""
    quality = ""
    position_errors = []
    for j, nucleotide in enumerate(subsequence):
        error = random.randint(0, 99)
        if error < errors[j]:
            wrong_nucleotide = random.choice([n for n in alphabet if n != nucleotide])
            read += wrong_nucleotide
            quality += "!"
            position_errors.append((j, nucleotide))
        else:
            read += nucleotide
            quality += "H"
    read_id = f"read_{i}"
    reads.append((read_id, read, quality, position_errors))
    
# make fastQ
with open("reads.fastq", "w") as f:
    for read in reads:
        f.write(f"@{read[0]}\n{read[1]}\n+\n{read[2]}\n")
        
# Trimmomatic
# Обработка ридов Trimmomatic
subprocess.run(["java", "-jar", "trimmomatic.jar", "SE", "reads.fastq", "trimmed_reads.fastq", "ILLUMINACLIP:TruSeq3-SE.fa:2:30:10", "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"])

# Анализ результата Trimmomatic
trimmed_reads = []
with open("trimmed_reads.fastq", "r") as f:
    lines = f.readlines()
    for i in range(0, len(lines), 4):
        read_id = int(lines[i].strip()[4:])
        trimmed_seq = lines[i+1].strip()
        trimmed_errors = reads[read_id][3]
        trimmed_errors_removed = [i for i in trimmed_errors if trimmed_seq[i] != reads[read_id][1][i]]
        trimmed_correct_removed = [i for i in range(len(trimmed_seq)) if trimmed_seq[i] != reads[read_id][1][i] and i not in trimmed_errors_removed]
        trimmed_reads.append((read_id, trimmed_seq, trimmed_errors_removed, trimmed_correct_removed))

correct_removed = 0
error_removed = 0
for read in trimmed_reads:
    correct_removed += len(read[3])
    error_removed += len(read[2])

print(f"Trimmomatic удалил {correct_removed} нуклеотидов, которые были считаны верно, и {error_removed} ошибочных нуклеотидов.")

# но Trimmomatic не работает корректно

""" 
Далее использую Lighter для восстановления ридов
для этого установил его и воспользовался lighter -r trimmed_reads.fastq -t 8 -k 21 -c 3 -o corrected_reads.fastq
-r или --reads: указывает путь к входному файлу FASTQ.
-t или --threads: указывает количество потоков для использования.
-k или --kmer-len: указывает длину k-меров, используемых в поиске ошибок.
-c или --correction-threshold: указывает пороговое значение для исправления ошибок. Этот параметр задает, сколько раз должен встречаться правильный нуклеотид в позиции, чтобы его можно было использовать для исправления ошибки.
-o или --output: указывает путь и имя файла, в который будут записаны скорректированные прочтения.
"""

""" 
Чтобы посчитать TP (true positive), FP (false positive), TN (true negative) и FN (false negative) 
я сравниваю скорректированные прочтения, полученные с помощью Lighter, с эталонными последовательностями.

Для этого использовал программу для выравнивания прочтений Bowtie2 чтобы сопоставить скорректированные прочтения 
с эталонной последовательностью. Затем, используя результаты выравнивания рассчитал метрики TP, FP, TN и FN.

Выравнивание сделал через команду bowtie2 -x reads.fastq -U corrected_reads.fastq -S aligned.sam

Для подсчета метрик использовал программу SAMtools, она считает в том числе все требуемые метрики
Сделал это через команду samtools flagstat aligned.sam

"""
