def greedy_motif_search(Dna, k, t):
    # Функция для создания профиля из набора k-меров
    def create_profile(motifs):
        profile = {'A': [], 'C': [], 'G': [], 'T': []}
        motif_length = len(motifs[0])
        motif_count = len(motifs)

        for nucleotide in profile:
            profile[nucleotide] = [0] * motif_length

        for i in range(motif_count):
            for j in range(motif_length):
                nucleotide = motifs[i][j]
                profile[nucleotide][j] += 1

        for nucleotide in profile:
            for j in range(motif_length):
                profile[nucleotide][j] /= motif_count

        return profile

    # Функция для поиска наиболее вероятного k-мера в строке, используя профиль
    def find_most_probable_kmer(dna_string, k, profile):
        max_prob = -1
        most_probable_kmer = ""

        for i in range(len(dna_string) - k + 1):
            kmer = dna_string[i:i + k]
            prob = 1

            for j in range(k):
                nucleotide = kmer[j]
                prob *= profile[nucleotide][j]

            if prob > max_prob:
                max_prob = prob
                most_probable_kmer = kmer

        return most_probable_kmer

    # Функция для вычисления "скора" набора k-меров
    def compute_score(motifs):
        consensus = ""
        motif_length = len(motifs[0])
        motif_count = len(motifs)
        score = 0

        for j in range(motif_length):
            count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}

            for i in range(motif_count):
                nucleotide = motifs[i][j]
                count[nucleotide] += 1

            max_count = max(count.values())
            consensus_nucleotides = [nucleotide for nucleotide in count if count[nucleotide] == max_count]
            consensus += consensus_nucleotides[0]  # Используем первый нуклеотид в случае неоднозначности

            score += motif_count - max_count

        return score, consensus

    # Инициализация BestMotifs первыми k-мерами из каждой строки
    best_motifs = [Dna[i][:k] for i in range(t)]

    # Поиск лучших k-меров
    for i in range(len(Dna[0]) - k + 1):
        motifs = [Dna[0][i:i + k]]

        for j in range(1, t):
            profile = create_profile(motifs)
            most_probable_kmer = find_most_probable_kmer(Dna[j], k, profile)
            motifs.append(most_probable_kmer)

        if compute_score(motifs)[0] < compute_score(best_motifs)[0]:
            best_motifs = motifs

    return best_motifs


# Пример использования
k = 3
t = 5
Dna = [
    "GGCGTTCAGGCA",
    "AAGAATCAGTCA",
    "CAAGGAGTTCGC",
    "CACGTCAATCAC",
    "CAATAATATTCG"
]

best_motifs = greedy_motif_search(Dna, k, t)
for motif in best_motifs:
    print(motif)
