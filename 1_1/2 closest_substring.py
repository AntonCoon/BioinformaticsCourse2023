def closest_substring(s1, s2):
    if len(s1) > len(s2):
        s2, s1 = s1, s2
    st_pos = 0
    min_ham = float('inf')
    sub_s = ''
    for i in range(len(s2)):
        j = i + len(s1)
        if j > len(s2):
            break
        tmp = hamming(s2[i:j], s1)
        if tmp < min_ham:
            min_ham = tmp
            st_pos = i
            sub_s = s2[i:j]
    return st_pos, sub_s, min_ham
