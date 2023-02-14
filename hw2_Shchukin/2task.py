  import Bio.Align as Align

matrix = Align.substitution_matrices.load("BLOSUM62")


def aff_gaps(a: str, b: str, substitution_matrix: list, alpha: int, betta: int):
    n, m = len(a), len(b)
    F = [[0] * (m + 1) for _ in range(n + 1)]
    D = [[0] * (m + 1) for _ in range(n + 1)]  # deletions matrix
    I = [[0] * (m + 1) for _ in range(n + 1)]  # insertions matrix

    for i in range(n + 1):  # fill in the first column
        F[0][i] = alpha + betta * i
        D[0][i] = alpha + betta * i
        I[0][i] = alpha + betta * i

    for j in range(m + 1):  # fill in the first line
        F[j][0] = alpha + betta * j
        D[j][0] = alpha + betta * j
        I[j][0] = alpha + betta * j

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = F[i - 1][j - 1] + substitution_matrix[a[i - 1]][b[j - 1]]
            D[i][j] = max(D[i-1][j] + betta, F[i-1][j] + alpha + betta)
            I[i][j] = max(I[i][j-1] + betta, F[i][j-1] + alpha + betta)
            F[i][j] = max(match, D[i][j], I[i][j])

    new_a = ''
    new_b = ''
    i = n
    j = m
    while i > 0 and j > 0:
        score = F[i][j]
        score_diag = F[i - 1][j - 1]
        score_up = max(I[i][j-1] + betta, F[i][j-1] + alpha + betta)
        score_left = max(D[i-1][j] + betta, F[i - 1][j] + alpha + betta)

        if score == score_diag + substitution_matrix[a[i - 1]][b[j - 1]]:
            new_a = a[i - 1] + new_a
            new_b = b[j - 1] + new_b
            i = i - 1
            j = j - 1

        elif score == score_left:
            new_a = a[i - 1] + new_a
            new_b = '_' + new_b
            i = i - 1

        elif score == score_up:
            new_a = '_' + new_a
            new_b = b[j - 1] + new_b
            j = j - 1

    while i > 0:
        new_a = a[i - 1] + new_a
        new_b = '_' + new_b
        i = i - 1

    while j > 0:
        new_a = '_' + new_a
        new_b = b[j - 1] + new_b
        j = j - 1

    return F[n][m], new_a, new_b


A = input()
B = input()
al = int(input())
bet = int(input())
print(aff_gaps(A, B, matrix, al, bet))
