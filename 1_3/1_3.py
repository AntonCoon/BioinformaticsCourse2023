# Вызовы Хиршберга (от 3 до 8 баллов)
# Реализуйте алгоритм, который принимал бы на вход две строки, штраф за удаления, вставки и несовпадения, а так-же цену совпадений.
# А возвращал бы дерево рекурсивных вызовов алгоритма Хоршберга. В вершинах дерева должны храниться пары строк, на которых будут
# происходить вызовы алгоритма Хиршберга.
# За 3 балла можно реализовать алгоритм используя квадрат памяти. При использовании линейного размера памяти задача стоит 8 баллов.
import numpy as np
import graphviz
from binarytree import build, build2


def is_match(s1: str, s2: str):
    if s1 == s2:
        return match
    else:
        return mismatch


def make_table(a: str, b: str):
    m = len(a)
    n = len(b)
    dp = [[i * ins for i in range(n + 1)], [0 for _ in range(n + 1)]]
    for i in range(1, m + 1):
        dp[1][0] = i * delete
        for j in range(1, n + 1):
            match_score = dp[0][j - 1] + is_match(a[i - 1], b[j - 1])
            delete_score = dp[0][j] + delete
            insert_score = dp[1][j - 1] + ins
            dp[1][j] = max(match_score, delete_score, insert_score)
        dp[0] = dp[1].copy()

    return dp[-1]



def hirshberg(a: str, b: str):
    a_len = len(a)
    b_len = len(b)

    mid = a_len // 2
    edit_distance = make_table(a[:mid], b)
    reversed_edit_distance = make_table(a[:mid - 1:-1], b[::-1])
    d = [edit_distance[i] + reversed_edit_distance[-1 - i] for i in range(len(edit_distance))]
    i = mid
    j = d.index(max(d))

    s1 = str(a[:i]) + " " + str(b[:j])
    result.append(s1)

    s2 = str(a[i:]) + " " + str(b[j:])
    result.append(s2)

    if len(b[:j]) > 1 and len(a[:i]) > 1:
        hirshberg(a[:i], b[:j])
    if len(b[j:]) > 1 and len(a[i:]) > 1:
        hirshberg(a[i:], b[j:])

    # if len(b[:j]) <= 1 and len(a[:i]) <= 1:
    #     result.append(None)
    # if len(b[j:]) <= 1 and len(a[i:]) <= 1:
    #     result.append(None)

result = []
result1 = []
result2 = []

a = "AGTACGCA"
b = "TATGC"

result.append("AGTACGCA" + " " + "TATGC")

delete = -2
ins = -2
match = 2
mismatch = -1

hirshberg(a, b)
print(result)
print(result.sort(key=lambda e: (e is None, len(e)), reverse=True))

binary_tree = build2(result)
print(binary_tree)
