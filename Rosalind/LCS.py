# https://rosalind.info/problems/ba5c/

from typing import Dict, Tuple

def reverse(s: str):
    ans = list(s)
    u, v = 0, len(ans) - 1
    while u < v:
        ans[u], ans[v] = ans[v], ans[u]
        u += 1
        v -= 1
    return "".join(ans)

def max_str(a: str, b: str):
    return a if len(a) > len(b) else b

def LCS_main(a: str, b: str):
    if not a or not b:
        return ""
    dp_i = (a, b)

    if dp_i in vs:
        return dp[dp_i]
    else:
        vs[dp_i] = 1

    if a[-1] == b[-1]:
        ans = a[-1] + LCS_main(a[:-1], b[:-1])
        dp[dp_i] = ans
        return ans

    ans = max_str(LCS_main(a[:-1], b), LCS_main(a, b[:-1]))
    dp[dp_i] = ans
    return ans

def find_LCS(a: str, b: str):
    return reverse(LCS_main(a, b))

seq1 = input()
seq2 = input()

dp: Dict[Tuple[str, str], str] = {}
vs: Dict[Tuple[str, str], int] = {}

print(find_LCS(seq1, seq2))