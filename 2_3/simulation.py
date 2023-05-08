import math
import random
import sys
from typing import Tuple

import numpy as np


ALPHABET = ['A', 'C', 'G', 'T']
SEQUENCE_LEN = 1_000
SUBSTR_LEN_E = 150
SUBSTR_LEN_S = 30
SEQ_CNT = 100
SEQ_MAX_ERROR_CNT = 66
READ_CNT = 20_000


# Speeding up series of random.choices error nucleotides
def letter_seq_simulation(correct_letter: str) -> Tuple[int, str]:
    errors_cnt = np.random.randint(0, SEQ_MAX_ERROR_CNT)
    error_alphabet = [letter for letter in ALPHABET if letter != correct_letter]
    error_distribution = [0] * 3
    error_distribution[0] = np.random.binomial(errors_cnt, 1 / 3)
    error_distribution[1] = np.random.binomial(errors_cnt - error_distribution[0], 1 / 2)
    error_distribution[2] = errors_cnt - (error_distribution[0] + error_distribution[1])
    seq_letters = [(SEQ_CNT - errors_cnt, correct_letter)] + list(zip(error_distribution, error_alphabet))
    return max(seq_letters)


# Meta file contains full correct sequence
def main(output_fastq_path, output_meta_path):
    sequence = "".join(random.choices(ALPHABET, k=SEQUENCE_LEN))
    with open(output_meta_path, 'w') as file:
        file.write(sequence)
    read_with_errors = 0
    with open(output_fastq_path, 'w') as file:
        for read_id in range(READ_CNT):
            substr_len = int(np.random.normal(SUBSTR_LEN_E, SUBSTR_LEN_S))
            start_pos = np.random.randint(0, SEQUENCE_LEN - substr_len - 1)
            seq_str = []
            quality_str = []
            for j in range(substr_len):
                max_letter_frequency, max_letter = letter_seq_simulation(sequence[start_pos + j])
                seq_str.append(max_letter)
                error_p = 1 - max_letter_frequency / SEQ_CNT
                quality = -10 * math.log10(max(1 / SEQ_CNT, error_p))  # Avoiding log(0)
                quality_str.append(chr(int(quality) + 33))
            read_str = "".join(seq_str)
            errors = sum(read_str[j] != sequence[start_pos + j] for j in range(substr_len))
            file.write(f"@FAKE.{read_id + 1} {read_id + 1} length={substr_len}\n")
            file.write(f"{read_str}\n")
            file.write(f"+ {start_pos} {errors}\n")
            file.write("".join(quality_str) + "\n")
            if errors > 0:
                read_with_errors += 1
    print(f"{read_with_errors} of {READ_CNT} reads contain errors")


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
