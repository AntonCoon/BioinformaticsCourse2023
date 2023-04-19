import sys
from utils import Read


def main(original_fastq_path, trimmed_fastq_path):
    survived_in_trimmed = set()
    with open(trimmed_fastq_path, 'r') as trimmed_file:
        while True:
            try:
                trimmed_read = Read(trimmed_file)
                survived_in_trimmed.add(trimmed_read.seq_id)
            except EOFError as e:
                break
    tp, tn, fp, fn = 0, 0, 0, 0
    with open(original_fastq_path, 'r') as original_file:
        while True:
            try:
                original_read = Read(original_file)
                if original_read.errors_cnt == 0 and original_read.seq_id in survived_in_trimmed:
                    tn += 1
                elif original_read.errors_cnt == 0 and original_read.seq_id not in survived_in_trimmed:
                    fp += 1
                elif original_read.errors_cnt > 0 and original_read.seq_id not in survived_in_trimmed:
                    tp += 1
                elif original_read.errors_cnt > 0 and original_read.seq_id in survived_in_trimmed:
                    fn += 1
            except EOFError as e:
                break
    print(f"FP={fp} - dropped incorrectly")
    print(f"TP={tp} - dropped correctly")
    print(f"FN={fn} - survived incorrectly")
    print(f"TN={tn} - survived correctly")
    print(f"Incorrect sequences before: {tp + fn}")
    print(f"Incorrect sequences after: {fp + fn}")


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
