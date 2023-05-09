import random
import numpy as np
alphabet = ["A", "T", "G", "C"]

def write_fastq(seq):
    with open('my_seq.fastq', 'w') as file:
        lines =[]
        for i in range(len(seq)):
            lines.append('@ Sequense #' + str(i+1) + '\n' + str(seq[i])+'\n'+'+'+'\n'+ "".join('~'*15)+'\n')
        file.write("\n".join(str(line) for line in lines))

sint_line = []
for i in range(100):
    sint_line.append(random.choice(alphabet))

seq = ''.join(sint_line)

my_reads = []
for i in range(85):
    my_reads.append(seq[i:i+15])

write_fastq(my_reads)

