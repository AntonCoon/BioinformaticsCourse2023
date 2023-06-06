# ba5c

if __name__ == '__main__':
    s1 = input()
    s2 = input()
    dp = [[0 for _ in range(len(s2) + 1)] for _ in range(len(s1) + 1)]
    par = [[None for _ in range(len(s2) + 1)] for _ in range(len(s1) + 1)]
    for i in range(len(s1) + 1):
        for j in range(len(s2) + 1):
            if j == 0 or i == 0:
                dp[i][j] = 0
            else:
                sim = s1[i - 1] == s2[j - 1]
                diag = dp[i - 1][j - 1] + sim
                if diag >= dp[i - 1][j] and diag >= dp[i][j - 1] and sim == True:
                    par[i][j] = (i - 1, j - 1)
                    dp[i][j] = diag
                elif dp[i - 1][j] >= diag and dp[i - 1][j] >= dp[i][j - 1]:
                    par[i][j] = (i - 1, j)
                    dp[i][j] = dp[i - 1][j]
                else:
                    par[i][j] = (i, j - 1)
                    dp[i][j] = dp[i][j - 1]
    result = []
    cell = (len(s1), len(s2))
    while True:
        new_cell = par[cell[0]][cell[1]]
        if new_cell is None:
            break
        if cell[0] - new_cell[0] == 1 and cell[1] - new_cell[1] == 1:
            result.append(s1[cell[0] - 1])
        cell = new_cell
    print("".join(reversed(result)))
