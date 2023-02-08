def hamming(s1, s2):
    l = 0
    for i, j in zip(s1, s2):
        if i != j:
            l += 1
    return l
