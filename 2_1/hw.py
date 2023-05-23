import subprocess
from Bio import SeqIO

import random
import numpy as np


num = 50000
# 1 задание

def simulate_reed_sequencing():
    alphabet = ['A', 'T', 'G', 'C']
    dna = ''.join(random.choices(alphabet, k=num))
    reads = []
    ground_truths = []
    for i in range(num):
        read_length = int(np.random.normal(250, 30))
        start = random.randint(0, len(dna) - read_length)
        read = dna[start:start+read_length]
        ground_truths.append(read)

        errors = []
        p_values = []
        for j in range(read_length):
            N = random.randint(0, 60)
            if N > 0:
                true_nucleotide = read[j]
                other_nucleotides = [n for n in ['A', 'T', 'G', 'C'] if n != true_nucleotide]
                read = read[:j] + random.choice(other_nucleotides) + read[j+1:]
                errors.append((j, true_nucleotide))
                p_values.append(0.1)
            else:
                p_values.append(0.01)

        p = len(errors) / read_length
        reads.append((read, errors, p_values))

    with open('reads.fastq', 'w') as f:
        for i, (read, errors, p_values) in enumerate(reads):
            f.write(f'@READ_{i}\n')
            f.write(f'{read}\n')
            f.write('+\n')
            qualities = []
            for j in range(len(read)):
                p = p_values[j]
                quality = -10 * np.log10(p)
                qualities.append(chr(int(quality) + 33))
            f.write(''.join(qualities) + '\n')

    return ground_truths
#2 задание
def get_data(file_input, file_output='output.fastq'):
    trimmomatic_command = [
        'java', '-jar', 'trimmomatic-0.39.jar', 'SE', '-phred33',
        file_input, file_output, 'LEADING:1', 'TRAILING:1',
        'SLIDINGWINDOW:5:2', 'MINLEN:250'
    ]
    subprocess.run(trimmomatic_command)

    file_output = file_output if file_output.endswith('.fastq') else file_output + '.fastq'
    data = list(SeqIO.parse(file_output, 'fastq'))
    return [str(record.seq) for record in data]


def evaluate_trimmomatic_performance(ground_truths):
    result = get_data('D:\\Study\\BioinformaticsCourse2023\\2_1\\reads.fastq')
    trimmed_reads = list(SeqIO.parse('output.fastq', 'fastq'))
    correct_removals = 0
    incorrect_removals = 0
    for i in range(len(trimmed_reads)):
        trimmed_read = str(trimmed_reads[i].seq)
        ground_truth_read = ground_truths[i]
        if len(trimmed_read) < len(ground_truth_read):
            for j in range(len(trimmed_read)):
                if trimmed_read[j] != ground_truth_read[j]:
                    incorrect_removals += 1
            correct_removals += len(ground_truth_read) - len(trimmed_read)
        else:
            for j in range(len(ground_truth_read)):
                if trimmed_read[j] != ground_truth_read[j]:
                    incorrect_removals += 1

    print(f'Correct removals: {correct_removals}')
    print(f'Incorrect removals: {incorrect_removals}')

ground_truths = simulate_reed_sequencing()
evaluate_trimmomatic_performance(ground_truths)

# 4 задание

def compare_reads(reads_file, output_file):
    TP = 0
    TN = 0
    FP = 0
    FN = 0

    with open(reads_file, 'r') as reads, open(output_file, 'r') as output:
        for read_line, output_line in zip(reads, output):
            if read_line.startswith('@'):
                continue
            for read_base, output_base in zip(read_line.strip(), output_line.strip()):
                if read_base == output_base:
                    if read_base == 'N':
                        FN += 1
                    else:
                        TN += 1
                else:
                    if read_base == 'N':
                        FP += 1
                    else:
                        TP += 1

    return TP, TN, FP, FN

tp, tn, fp, fn = compare_reads('reads.fastq', 'pollux-master/output/output.fastq.low')
print(f"TP: {tp}, TN: {tn}, FP: {fp}, FN: {fn}")