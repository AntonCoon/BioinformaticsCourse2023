from Levenshtein import distance


def get_consensus_string(str1, str2):
    m, n = len(str1), len(str2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1):
        dp[i][0] = i * del_penalty
    for j in range(n + 1):
        dp[0][j] = j * ins_penalty
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if str1[i - 1] == str2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + match_price
            else:
                dp[i][j] = min(dp[i - 1][j] + del_penalty,
                               dp[i][j - 1] + ins_penalty,
                               dp[i - 1][j - 1] + mismatch_penalty)
    i, j = m, n
    aligned_str1 = []
    aligned_str2 = []
    while i > 0 or j > 0:

        if i > 0 and dp[i][j] == dp[i - 1][j] + del_penalty:
            aligned_str1.append(str1[i - 1])
            aligned_str2.append('-')
            i -= 1
        elif j > 0 and dp[i][j] == dp[i][j - 1] + ins_penalty:
            aligned_str1.append('-')
            aligned_str2.append(str2[j - 1])
            j -= 1
        else:
            aligned_str1.append(str1[i - 1])
            aligned_str2.append(str2[j - 1])
            i -= 1
            j -= 1
    consensus_string = []
    for c1, c2 in zip(reversed(aligned_str1), reversed(aligned_str2)):
        if c1 == c2:
            consensus_string.append(c1)
        elif c1 == '-':
            consensus_string.append(c2)
        elif c2 == '-':
            consensus_string.append(c1)
        else:

            consensus_string.append(c1)

    return ''.join(consensus_string), ''.join(reversed(aligned_str1)), ''.join(reversed(aligned_str2))


def greedy_multiple_alignment(strings):
    alignments = [[s] for s in strings]
    while len(strings) > 1:
        min_dist = float('inf')
        min_pair = None
        for i in range(len(strings)):
            for j in range(i + 1, len(strings)):
                dist = distance(strings[i], strings[j])
                if dist < min_dist:
                    min_dist = dist
                    min_pair = (i, j)

        consensus_string, aligned_str1, aligned_str2 = get_consensus_string(strings[min_pair[0]], strings[min_pair[1]])
        strings[min_pair[0]] = consensus_string
        del strings[min_pair[1]]
        alignments[min_pair[0]] = [aligned_str1] * len(alignments[min_pair[0]]) + [aligned_str2] * len(
            alignments[min_pair[1]])
        del alignments[min_pair[1]]
    return alignments[0]



strings = ['ACTG', 'ACGT', 'AGCT']
del_penalty = 1
ins_penalty = 1
mismatch_penalty = 1
match_price = 0
alignment = greedy_multiple_alignment(strings)
print(*alignment)