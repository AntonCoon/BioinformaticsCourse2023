import numpy as np
from Levenshtein import distance as lev_distance
from Bio import pairwise2
from Bio.Align import substitution_matrices


'''
Если вдруг не работает с pairwise2 (разработчики пакета Bio вроде как потерли эту фукцию с последней версии пакета), то вот переделка через другую библиотеку

from NWalign import align

def align_sequences(sequences):
    matrix = substitution_matrices.load("BLOSUM62")  # выбираем матрицу замен
    gap_open = -10  # штраф за гэп
    gap_extend = -0.5  # штраф за продолжение гэпа
    aligned_sequences = []
    for i in range(len(sequences)-1):
        seq1 = sequences[i]
        seq2 = sequences[i+1]
        # глобальное выравнивание с матрицей замен и штрафами за гэпы
        alignments = align(seq1, seq2, matrix, gap_open, gap_extend)
        # добавляем выровненную последовательность в список
        aligned_sequences.append(alignments[0])
        # добавляем выровненную последовательность в список
        aligned_sequences.append(alignments[1])
    return aligned_sequences
'''



def align_sequences(sequences):
    matrix = substitution_matrices.load("BLOSUM62")  # выбираем матрицу замен
    aligned_sequences = []
    for i in range(len(sequences)-1):
        seq1 = sequences[i]
        seq2 = sequences[i+1]
        # глобальное выравнивание, тут штраф за гэп -10 и за продолжение гэпа -0.5
        alignments = pairwise2.align.globalds(seq1, seq2, matrix, -10, -0.5)
        # добавляем выровненную последовательность в список
        aligned_sequences.append(alignments[0][0])
        # добавляем выровненную последовательность в список
        aligned_sequences.append(alignments[0][1])
    return aligned_sequences


def get_profile(aligned_sequence):
    profile = []
    for i in range(len(aligned_sequence[0])):
        count = {}
        for sequence in aligned_sequence:
            char = sequence[i]
            if char in count:
                count[char] += 1
            else:
                count[char] = 1
        profile.append(count)
    return profile


# на вход подаются уже выровненные строки
def generate_consensus(aligned_sequences):
    consensus = ""
    profile = get_profile(aligned_sequences)
    for count in profile:
        consensus += max(count, key=count.get)
    return consensus


def replace_with_consensus(strings, insertion, deletion):
    # Находим две самые близкие строки
    min_distance = -10**10
    min_pair = (0,0)
    for i in range(len(strings)):
        for j in range(i+1, len(strings)):
            distance = lev_distance(strings[i], strings[j], weights=(
                insertion, deletion, 1))
            if distance < min_distance:
                min_distance = distance
                min_pair = (i, j)
    # Получаем консенсусную строку, предварительно выровняв строки
    consensus = generate_consensus(align_sequences([strings[i], strings[j]]))
    # Заменяем две строки на их консенсусную строку
    new_strings = []
    for i, s in enumerate(strings):
        if i not in min_pair:
            new_strings.append(s)
        else:
            new_strings.append(consensus)

    return new_strings


def multiple_alignment(strings, insertion=-1, deletion=-1):
    while len(strings) > 1:
        strings = replace_with_consensus(
            strings, insertion, deletion)
    return align_sequences(strings)


strings = ["ACGTACGT", "CGTACGT", "GTACGT", "AACGTT"]
# strings = ['AAG', 'AAT', 'AGT']
aligned_strings = multiple_alignment(strings)
for s in aligned_strings:
    print(s)
