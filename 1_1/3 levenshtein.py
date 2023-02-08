def levenshtein(s1, s2):
    if len(s1) > len(s2):
        s2, s1 = s1, s2
    d = [[0 for _ in range(len(s1) + 1)] for _ in range(2)] #len(s1) <= len(s2)
    for j in range(len(s2) + 1):
        for i in range(len(s1) + 1):
            if min(i, j) == 0:
                d[j%2][i] = max(i, j)
            elif s1[i-1] == s2[j-1]:
                d[j%2][i] = d[j%2-1][i-1]
            else:
                d[j%2][i] = min(d[j%2-1][i], d[j%2][i-1], d[j%2-1][i-1]) + 1
    return d[len(s2)%2][-1]
