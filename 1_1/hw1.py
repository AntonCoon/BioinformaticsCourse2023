from typing import Tuple

from Bio import SeqIO


def hamming_distance(a: str, b: str) -> int:
    return sum(c1 != c2 for c1, c2 in zip(a, b))


def nearest_substring(a: str, b: str) -> Tuple[int, str, int]:
    if len(a) < len(b):
        a, b = b, a
    ans_dist, ans_index = min((hamming_distance(a[i:], b), i) for i in range(len(a) - len(b) + 1))
    return ans_index, a[ans_index:ans_index+len(b)], ans_dist


def levenshtein_distance(a: str, b: str) -> int:
    if len(a) < len(b):
        a, b = b, a
    dp = [
        [0 for _ in range(len(b) + 1)],
        [i for i in range(len(b) + 1)]
    ]
    for i in range(1, len(a) + 1):
        dp[0] = dp[1].copy()
        dp[1][0] = i
        for j in range(1, len(b) + 1):
            dp[1][j] = min(
                dp[0][j] + 1,
                dp[1][j - 1] + 1,
                dp[0][j - 1] + (a[i - 1] != b[j - 1])
            )
    return dp[1][len(b)]


if __name__ == '__main__':
    print("Test #1")
    records = list(SeqIO.parse('data/f8.fasta', 'fasta'))
    s1 = str(records[0].seq)
    s2 = str(records[1].seq)
    print("Nearest substring:", nearest_substring(s1, s2))
    print("Levenshtein distance:", levenshtein_distance(s1, s2))

    print("Test #2")
    records = list(SeqIO.parse('data/gattaca.fasta', 'fasta'))
    s1 = str(records[0].seq)
    s2 = str(records[1].seq)
    print("Hamming distance:", hamming_distance(s1, s2))
    print("Nearest substring:", nearest_substring(s1, s2))
    print("Levenshtein distance:", levenshtein_distance(s1, s2))
