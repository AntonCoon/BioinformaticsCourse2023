def neighbor_joining(D, n):
    if n == 2:
        T = [[0, 1, D[0][1]]]
        return T
    D_prime = neighbor_joining_matrix(D)
    min_val = float('inf')
    for i in range(len(D_prime)):
        for j in range(len(D_prime[0])):
            if i != j and D_prime[i][j] < min_val:
                min_val = D_prime[i][j]
                min_i = i
                min_j = j
    delta = (total_distance(D, min_i) - total_distance(D, min_j)) / (n - 2)
    limb_length_i = 0.5 * (D[min_i][min_j] + delta)
    limb_length_j = 0.5 * (D[min_i][min_j] - delta)
    m = len(D)
    D.append([0] * (m + 1))
    for k in range(m):
        D[k].append(0.5 * (D[k][min_i] + D[k][min_j] - D[min_i][min_j]))
        D[m][k] = D[k][m]
    D[m][m] = 0
    for i in sorted([min_i, min_j], reverse=True):
        del D[i]
        for row in D:
            del row[i]
    T = neighbor_joining(D, n - 1)
    T.append([m, min_i, limb_length_i])
    T.append([m, min_j, limb_length_j])
    return T

def neighbor_joining_matrix(D):
    n = len(D)
    total_distances = [total_distance(D, i) for i in range(n)]
    D_prime = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                D_prime[i][j] = (n - 2) * D[i][j] - total_distances[i] - total_distances[j]
    return D_prime

def total_distance(D, i):
    return sum(D[i])

D = [
 [0, 23, 27, 20],
 [23, 0, 30, 28],
 [27, 30, 0, 30],
 [20, 28, 30, 0]
]
n = len(D)
T = neighbor_joining(D,n)
#print(T)
for row in T:
    print(f'{row[0]}->{row[1]}:{row[2]:}')
