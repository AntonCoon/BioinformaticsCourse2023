def count_occurrences(motifs):
    count = {'A': [], 'C': [], 'G': [], 'T': []}
    k = len(motifs[0])
    for i in range(k):
        occurrences = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for motif in motifs:
            occurrences[motif[i]] += 1
        count['A'].append(occurrences['A'])
        count['C'].append(occurrences['C'])
        count['G'].append(occurrences['G'])
        count['T'].append(occurrences['T'])
    return count

def profile_most_probable_kmer(dna, k, profile):
    max_prob = -1
    most_probable = ""
    for i in range(len(dna) - k + 1):
        kmer = dna[i:i+k]
        prob = 1
        for j, nucleotide in enumerate(kmer):
            prob *= profile[nucleotide][j]
        if prob > max_prob:
            max_prob = prob
            most_probable = kmer
    return most_probable

def score(motifs):
    count = count_occurrences(motifs)
    consensus = ""
    k = len(motifs[0])
    for i in range(k):
        max_count = -1
        most_frequent = ""
        for nucleotide in ['A', 'C', 'G', 'T']:
            if count[nucleotide][i] > max_count:
                max_count = count[nucleotide][i]
                most_frequent = nucleotide
        consensus += most_frequent
    score = 0
    for motif in motifs:
        for i in range(k):
            if motif[i] != consensus[i]:
                score += 1
    return score

def greedy_motif_search(dna, k, t):
    best_motifs = [string[:k] for string in dna]
    n = len(dna[0])
    for i in range(n - k + 1):
        motifs = [dna[0][i:i+k]]
        for j in range(1, min(t, len(dna))): 
            profile = count_occurrences(motifs)
            motifs.append(profile_most_probable_kmer(dna[j], k, profile))
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
    return best_motifs


# Read input
k, t = 12, 25
dna = []
with open('rosalind_ba2d.txt', 'r') as file:
    dna = [line.strip() for line in file.readlines()]
    
dna_new = dna[1:]


result = greedy_motif_search(dna_new, k, t)


with open('result.txt', 'w') as file:
    for motif in result:
        file.write(motif + '\n')
