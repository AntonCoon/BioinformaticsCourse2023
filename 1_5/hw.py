import numpy as np
from functools import lru_cache
from typing import List, Tuple


def nw(a: str, b: str, del_pen: float, ins_pen: float, mismatch_pen: float, match: float) -> np.ndarray:
    a = "-" + a 
    b = "-" + b
    f = np.zeros((2, len(b)))
    for i in range(len(a)):
        for j in range(len(b)):
            if i == 0:
                f[1, j] = ins_pen*j
            elif j == 0:
                f[1, 0] = del_pen*i
            else:
                delete = f[0, j] + del_pen
                insert = f[1, j-1] + ins_pen
                f[1, j] = max(delete, insert, f[0, j-1] + (mismatch_pen if a[i]!=b[j] else match))
        f[0, :] = f[1, :]
    return f[1, :]

def hirschberg(a: str, b: str, del_pen: float, ins_pen: float, mismatch_pen: float, match: float) -> str:
    if len(a) <= 1 or len(b) <= 1:
        return a if len(a) == 1 and len(b) == 1 else a + b
    mid = len(a)//2
    scores = nw(a[:mid], b, del_pen, ins_pen, mismatch_pen, match) + \
        nw(a[mid:][::-1], b[::-1], del_pen, ins_pen, mismatch_pen, match)[::-1]
    cut_ind = np.argmax(scores)
    lc = hirschberg(a[:mid], b[:cut_ind], del_pen, ins_pen, mismatch_pen, match)
    rc = hirschberg(a[mid:], b[cut_ind:], del_pen, ins_pen, mismatch_pen, match)
    return lc + rc

@lru_cache(maxsize=None)
def dist_and_alignment(s1: str, s2: str, del_pen: float, ins_pen: float, mismatch_pen: float, match: float) -> Tuple[float, str]:
    score = nw(s1, s2, del_pen, ins_pen, mismatch_pen, match)[-1]
    consensus = hirschberg(s1, s2, del_pen, ins_pen, mismatch_pen, match)
    return score, consensus

def find_max_score(all_seqs: List[str], curr_in_use: List[bool], del_pen: float, ins_pen: float, mismatch_pen: float, match: float) -> Tuple[int, int, str]:
    min_score = -float('inf')
    first_seq_id, second_seq_id = 0, 0
    best_consensus = ''
    for i in range(len(all_seqs)):
        for j in range(i+1, len(all_seqs)):
            if curr_in_use[i] and curr_in_use[j]:
                score, consensus = dist_and_alignment(all_seqs[i], all_seqs[j], del_pen, ins_pen, mismatch_pen, match)
                if score > min_score:
                    min_score = score
                    first_seq_id, second_seq_id = i, j
                    best_consensus = consensus
    return first_seq_id, second_seq_id, best_consensus

def multiple_seq_alignment(gens: List[str], del_pen: float, ins_pen: float, mismatch_pen: float, match: float) -> str:
    all_seqs = gens.copy()
    curr_in_use = [True for _ in range(len(gens))]
    for iter in range(len(gens)-1):
        i, j, consensus = find_max_score(all_seqs, curr_in_use, del_pen, ins_pen, mismatch_pen, match)
        all_seqs.append(consensus)
        curr_in_use[i] = curr_in_use[j] = False
        curr_in_use.append(True)
    return all_seqs[-1]


def main():
    del_pen = -2
    ins_pen = -2
    mismatch_pen = -2
    match = 2
    gens = ["GATTACA", "GTTACA", "GATTAA", "GATTATA"]
    print(multiple_seq_alignment(gens, del_pen, ins_pen, mismatch_pen, match))


if __name__ == "__main__":
    main()