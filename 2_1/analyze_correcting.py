import sys
from utils import Read


def main(original_fastq_path, corrected_fastq_path, meta_path):
    with open(meta_path, 'r') as file:
        sequence = file.readline().strip()
    tp, tn, fp, fn = 0, 0, 0, 0
    with open(original_fastq_path, 'r') as original_file:
        with open(corrected_fastq_path, 'r') as corrected_file:
            while True:
                try:
                    original_read = Read(original_file)
                    corrected_read = Read(corrected_file)
                    start_pos = original_read.start_pos
                    substr_len = len(original_read.seq)
                    meta_seq = sequence[start_pos:start_pos + substr_len]
                    if original_read.seq == meta_seq and corrected_read.seq == meta_seq:
                        tn += 1
                    elif original_read.seq == meta_seq and corrected_read.seq != meta_seq:
                        fp += 1
                    elif original_read.seq != meta_seq and corrected_read.seq == meta_seq:
                        tp += 1
                    elif original_read.seq != meta_seq and corrected_read.seq != meta_seq:
                        fn += 1
                except EOFError as e:
                    break
    print(f"FP={fp} - correct before, incorrect after")
    print(f"TP={tp} - incorrect before, correct after")
    print(f"FN={fn} - incorrect before, incorrect after")
    print(f"TN={tn} - correct before, correct after")
    print(f"Incorrect sequences before: {tp + fn}")
    print(f"Incorrect sequences after: {fp + fn}")


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
