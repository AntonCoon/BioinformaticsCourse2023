def count_frequence(motifs):
    count = {'A': [], 'C': [], 'G': [], 'T': []}
    k = len(motifs[0])
    for i in range(k):
        frequences = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for motif in motifs:
            frequences[motif[i]] += 1
        count['A'].append(frequences['A'])
        count['C'].append(frequences['C'])
        count['G'].append(frequences['G'])
        count['T'].append(frequences['T'])
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
    count = count_frequence(motifs)
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
            profile = count_frequence(motifs)
            motifs.append(profile_most_probable_kmer(dna[j], k, profile))
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
    return best_motifs

dna = []
with open('rosalind_ba2d.txt', 'r') as file:
    dna = [line.strip() for line in file.readlines()]

k, t = map(int, dna[0].split())
dna_new = dna[1:]
result = greedy_motif_search(dna_new, k, t)
for motif in result:
    print(motif)