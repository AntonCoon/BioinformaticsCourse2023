import random
import csv
def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = [[0 for i in range(k)] for j in range(4)]
    for i in range(t):
        for j in range(k):
            if Motifs[i][j] == 'A':
                profile[0][j] += 1
            elif Motifs[i][j] == 'C':
                profile[1][j] += 1
            elif Motifs[i][j] == 'G':
                profile[2][j] += 1
            elif Motifs[i][j] == 'T':
                profile[3][j] += 1
    for i in range(4):
        for j in range(k):
            profile[i][j] = profile[i][j] / t
    return profile


def ProfileMostProbablePattern(Text, k, Profile):
    n = len(Text)
    maxprob = -1
    for i in range(n - k + 1):
        prob = 1
        Pattern = Text[i:i + k]
        for j in range(k):
            if Pattern[j] == 'A':
                prob *= Profile[0][j]
            elif Pattern[j] == 'C':
                prob *= Profile[1][j]
            elif Pattern[j] == 'G':
                prob *= Profile[2][j]
            elif Pattern[j] == 'T':
                prob *= Profile[3][j]
        if prob > maxprob:
            maxprob = prob
            MostProbablePattern = Pattern
    return MostProbablePattern


def Score(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    score = 0
    for j in range(k):
        count = [0, 0, 0, 0]
        for i in range(t):
            if Motifs[i][j] == 'A':
                count[0] += 1
            elif Motifs[i][j] == 'C':
                count[1] += 1
            elif Motifs[i][j] == 'G':
                count[2] += 1
            elif Motifs[i][j] == 'T':
                count[3] += 1
        score += t - max(count)
    return score
def RandomizedMotifSearch(Dna, k, t):
    M = []
    for i in range(t):
        n = len(Dna[i])
        r = random.randint(0,n-k)
        M.append(Dna[i][r:r+k])
    BestM = list(M)
    while True:
        profile = Profile(M)
        M = []
        for i in range(t):
            M.append(ProfileMostProbablePattern(Dna[i],k,profile))
        if Score(M) < Score(BestM):
            BestM = list(M)
        else:
            return BestM



with open('input_4.txt', 'r') as f:
    lines = f.read().strip().split('\n')
    k, t = map(int, lines[0].split())
    Dna = lines[1:]

best_motifs = RandomizedMotifSearch(Dna, k, t)
best_score = Score(best_motifs)

for i in range(1000):
    motifs = RandomizedMotifSearch(Dna, k, t)
    score = Score(motifs)
    if score < best_score:
        best_motifs = motifs
        best_score = score

for motif in best_motifs:
    print(motif)