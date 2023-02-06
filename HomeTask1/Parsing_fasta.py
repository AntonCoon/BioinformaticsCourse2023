# парсинг .fasta файлов

from Bio import SeqIO


def read_fasta(file_name):
    handle = open(file_name)
    seq_list = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    sequences_ids = []
    sequences_lines = []

    for i in range(len(seq_list)):
        sequences_ids.append(str(seq_list[i].id))
        sequences_lines.append(str(seq_list[i].seq))

    return sequences_ids, sequences_lines

file_name = input('Print path to file:', )