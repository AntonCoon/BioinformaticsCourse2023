import Bio.Align as Align

matrix = Align.substitution_matrices.load("BLOSUM62")


def nidlman_vunsh(a: str, b: str, substitution_matrix: list, gap_cost: int):
    n, m = len(a), len(b)
    F = [[0] * (m+1) for _ in range(n+1)]

    for i in range(n+1):  # fill in the first column
        F[i][0] = gap_cost * i

    for j in range(m+1):  # fill in the first line
        F[0][j] = gap_cost * j

    for i in range(1, n+1):
        for j in range(1, m+1):
            match = F[i-1][j-1] + substitution_matrix[a[i-1]][b[j-1]]
            delete = F[i-1][j] + gap_cost
            insert = F[i][j-1] + gap_cost
            F[i][j] = max(match, delete, insert)

    new_a = ''
    new_b = ''
    i = n
    j = m
    while i > 0 and j > 0:
        score = F[i][j]
        score_diag = F[i-1][j-1]
        score_up = F[i][j-1]
        score_left = F[i-1][j]

        if score == score_diag + substitution_matrix[a[i-1]][b[j-1]]:
            new_a = a[i-1] + new_a
            new_b = b[j-1] + new_b
            i = i-1
            j = j-1

        elif score == score_left + gap_cost:
            new_a = a[i-1] + new_a
            new_b = '_' + new_b
            i = i - 1

        elif score == score_up + gap_cost:
            new_a = '_' + new_a
            new_b = b[j-1] + new_b
            j = j - 1

    while i > 0:
        new_a = a[i-1] + new_a
        new_b = '_' + new_b
        i = i - 1

    while j > 0:
        new_a = '_' + new_a
        new_b = b[j-1] + new_b
        j = j - 1
        
    return F[n][m], new_a, new_b

A = input()
B = input()
d = int(input())
print(nidlman_vunsh(A, B, matrix, d))
