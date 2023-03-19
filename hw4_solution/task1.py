# Import the SeqIO module from Biopython
from Bio import SeqIO

# Read the DNA sequence from the input file
with open("islands.fasta", "r") as file:
    records1 = list(SeqIO.parse(file, "fasta"))
    sequence1 = str(records1[0].seq)

with open("nonislands.fasta", "r") as file:
    records2 = list(SeqIO.parse(file, "fasta"))
    sequence2 = str(records2[0].seq)

# Count the occurrences of each nucleotide within and outside the CpG islands
def count_single(sequence: list, total_count: dict, cg_count: int):
    for i in range(len(sequence)):
        nucleotide = sequence[i]
        total_count[nucleotide] += 1
        if nucleotide == "C" and i < len(sequence) - 1 and sequence[i+1] == "G":
            cg_count += 1

    return total_count, cg_count


total_count_inside = {"A": 0, "C": 0, "G": 0, "T": 0}
cg_count_inside = 0
total_count_outside = {"A": 0, "C": 0, "G": 0, "T": 0}
cg_count_outside = 0
total_count_inside, cg_count_inside = count_single(
    sequence1, total_count_inside, cg_count_inside)
total_count_outside, cg_count_outside = count_single(
    sequence2, total_count_outside, cg_count_outside)
print(total_count_inside, 'CG count inside:', cg_count_inside)
print(total_count_outside, 'CG count outside:', cg_count_outside)

# Count the occurrences of pairs of nucleotides within and outside the CpG islands
def count_pairs(sequence: list, total_count: dict):
    for i in range(0, len(sequence)-1, 2):
        pair = sequence[i:i+2]
        total_count[pair] += 1
    return total_count


# AGCT
pair_inside = {'AA': 0, 'AG': 0, 'GA': 0, 'AC': 0, 'CA': 0, 'AT': 0, 'TA': 0,
               'GG': 0, 'GC': 0, 'CG': 0, 'GT': 0, 'TG': 0, 'CC': 0, 'CT': 0,
               'TC': 0, 'TT': 0}
pair_outside = {'AA': 0, 'AG': 0, 'GA': 0, 'AC': 0, 'CA': 0, 'AT': 0, 'TA': 0,
                'GG': 0, 'GC': 0, 'CG': 0, 'GT': 0, 'TG': 0, 'CC': 0, 'CT': 0,
                'TC': 0, 'TT': 0}
pair_inside = count_pairs(sequence1, pair_inside)
pair_outside = count_pairs(sequence2, pair_outside)
print(pair_inside)
print(pair_outside)


# Calculate the frequencies of each nucleotide within and outside the CpG islands
total_inside = sum(total_count_inside.values())
total_outside = sum(total_count_outside.values())
freq_inside = {nucleotide: round(count/total_inside, 4) for nucleotide,
               count in total_count_inside.items()}
freq_outside = {nucleotide: round(count/total_outside, 4) for nucleotide,
                count in total_count_outside.items()}

print(freq_inside)
print(freq_outside)


# Calculate the frequencies of each nucleotide within and outside the CpG islands
total_inside_pair = sum(pair_inside.values())
total_outside_pair = sum(pair_outside.values())
pair_freq_inside = {pair: round(count/total_inside_pair, 4) for pair,
                    count in pair_inside.items()}
pair_freq_outside = {pair: round(count/total_outside_pair, 4) for pair,
                     count in pair_outside.items()}

print(pair_freq_inside)
print(pair_freq_outside)
