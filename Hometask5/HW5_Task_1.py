# Жадный алгоритм множественного выравнивания (7 баллов) Реализуйте алгоритм, который принимал бы на вход массив строк,
# штраф за удаления, вставки и несовпадения, а так-же цену совпадений. А возвращал бы множественное выравнивание.
# На первом шаге алгоритм должен выбрать две самые близкие по расстоянию Левенштейна строки и заменить их консенснусной
# строкой. При следующих шагах алгоритма выравниваться между собой могут так же и консенснусные строки.
# При этом стоит хранить для каждой строки не только ее саму но и профиль множественного выравнивания, чтобы в итоге
# правильно пересчитывать консенсус.
# Результат работы алгоритма - массив строк, соответствующий некоторому множественному выравниванию.

def hirsh_dist(a, b, ins, delete, match, mismatch):
    dp = [[0 for _ in range(len(b) + 1)],
        [i * ins for i in range(len(b) + 1)]]
    for i in range(1, len(a) + 1):
        dp[0] = dp[1].copy()
        dp[1][0] = i * delete
        for j in range(1, len(b) + 1):
            dp[1][j] = max(
                dp[0][j] + delete,
                dp[1][j - 1] + ins,
                dp[0][j - 1] + (mismatch if a[i - 1] != b[j - 1] else match)
            )
    return dp[1]

def hirschberg_algo(a, b, insert, delete, mismatch, match):
    if len(a) <= 1 or len(b) <= 1:
        return a if len(a) == 1 and len(b) == 1 else a + b

    a_split = len(a) // 2
    a_l, a_r = a[:a_split], a[a_split:]

    dp_l = hirsh_dist(a_l, b, insert, delete, mismatch, match)
    dp_r = hirsh_dist(a_r[::-1], b[::-1], insert, delete, mismatch, match)[::-1]
    common_dp = [d_l + d_r for d_l, d_r in zip(dp_l, dp_r)]
    b_split_index = max(range(len(common_dp)), key=common_dp.__getitem__)
    b_l, b_r = b[:b_split_index], b[b_split_index:]

    left = hirschberg_algo(a_l, b_l, insert, delete, mismatch, match)
    right = hirschberg_algo(a_r, b_r, insert, delete, mismatch, match)
    return left + right

def edit_distance(line_1, line_2):
    lenght_1, lenght_2 = len(line_1), len(line_2)

    if lenght_1 > lenght_2:
        line_1, line_2 = line_2, line_1
        lenght_1, lenght_2 = lenght_2, lenght_1

    rmmbr = [[0 for i in range(lenght_1 + 1)] for j in range(2)]

    for i in range(lenght_1 + 1):
        rmmbr[0][i] = i

    for i in range(1, lenght_2 + 1):
        for j in range(0, lenght_1 + 1):
            if j == 0:
                rmmbr[i % 2][j] = i
            elif (line_1[j - 1] == line_2[i - 1]):
                rmmbr[i % 2][j] = rmmbr[(i - 1) % 2][j - 1]
            else:
                rmmbr[i % 2][j] = (1 + min(rmmbr[(i - 1) % 2][j], min(rmmbr[i % 2][j - 1], rmmbr[(i - 1) % 2][j - 1])))

    return rmmbr[lenght_2 % 2][lenght_1]

def multi_alignment(seq: list, match, insert, delete, mismatch):
    consensus = seq.copy()
    n = len(consensus)
    for _ in range(n-1):
        pair = None
        min_dist = None
        for i in range(len(consensus)):
            for j in range(i):
                dist = edit_distance(consensus[i], consensus[j])
                if min_dist is None or dist < min_dist:
                    pair = i, j
                    min_dist = dist

        i, j = pair
        best_cons = hirschberg_algo(consensus[i], consensus[j], insert, delete, mismatch, match)
        print(consensus[i], consensus[j], min_dist, best_cons)
        consensus.pop(i)
        consensus.pop(j)
        consensus.append(best_cons)
    return consensus[0]


sequenses = ["ATA", "AT", "TAT"]
insert_cost, delete_cost, mismatch_cost = -1, -1, -1
match_cost = 1

alig_sequenses = multi_alignment(sequenses, match_cost, insert_cost, delete_cost, mismatch_cost)
print(alig_sequenses)
