import graphviz
from typing import NamedTuple, Optional, Tuple, List


# String with modified substring function
# Substring time/memory complexity is O(1) instead of O(n)
class StrHolder:
    # view for a[l ... r] = a[l : r + 1]
    def __init__(self, a: str, l: int = 0, r: int = None, reversed: bool = False):
        if r is None:
            r = len(a) - 1
        assert 0 <= l <= r + 1 <= len(a), "out of bounds"
        self._a = a
        self._l = l
        self._r = r
        self._reversed = reversed

    def __len__(self):
        return self._r - self._l + 1

    def _get_real_index(self, ind):
        if not self._reversed:
            return self._l + ind
        else:
            return self._r - ind

    def _get_real_substr_idx(self, l, r):
        new_l, new_r = self._get_real_index(l), self._get_real_index(r)
        if not self._reversed:
            return new_l, new_r
        else:
            return new_r, new_l

    def __getitem__(self, ind):
        if type(ind) is slice:
            start = ind.start if ind.start is not None else 0
            end = ind.stop - 1 if ind.stop is not None else self.__len__() - 1
            step = ind.step if ind.step is not None else 1
            assert step in [-1, 1], "step must be either -1 or 1"
            reversed = step == -1

            new_l, new_r = self._get_real_substr_idx(start, end)
            return StrHolder(self._a, new_l, new_r, self._reversed != reversed)
        else:
            assert ind >= 0, "negative index is not allowed"
            return self._a[self._get_real_index(ind)]

    # O(n) function, but it's used only for printing
    def __str__(self):
        step = 1 if not self._reversed else -1
        return self._a[self._l:self._r + 1:step]


class HirschbergConfig(NamedTuple):
    del_cost: int
    ins_cost: int
    match_cost: int
    mismatch_cost: int


class TracedHirschbergReport:
    def __init__(self, a: StrHolder, b: StrHolder, children: Optional[Tuple]):
        self.a = a
        self.b = b
        self.children = children

    def get_dot(self):
        dot = graphviz.Digraph('hirschberg', comment='Hirschberg algorithm')

        def _fill_dot(v, index=1):
            node_idx = f'V_{index}'
            dot.node(node_idx, f'({v.a},{v.b})')
            if v.children is not None:
                left_child, right_child = v.children
                left_idx, right_idx = 2 * index, 2 * index + 1
                _fill_dot(left_child, left_idx)
                _fill_dot(right_child, right_idx)
                dot.edge(node_idx, f'V_{left_idx}')
                dot.edge(node_idx, f'V_{right_idx}')
        _fill_dot(self)

        return dot

    def view(self):
        self.get_dot().view()

    def save(self, filename):
        self.get_dot().render(filename)


def calculate_distance_matrix(a: StrHolder, b: StrHolder, config: HirschbergConfig) -> List[int]:
    dp = [
        [0 for _ in range(len(b) + 1)],
        [i * config.ins_cost for i in range(len(b) + 1)]
    ]
    for i in range(1, len(a) + 1):
        dp[0] = dp[1].copy()
        dp[1][0] = i * config.del_cost
        for j in range(1, len(b) + 1):
            dp[1][j] = max(
                dp[0][j] + config.del_cost,
                dp[1][j - 1] + config.ins_cost,
                dp[0][j - 1] + (config.mismatch_cost if a[i - 1] != b[j - 1] else config.match_cost)
            )
    return dp[1]


def traced_hirschberg(a: StrHolder, b: StrHolder, config: HirschbergConfig) -> TracedHirschbergReport:
    if len(a) <= 1 or len(b) <= 1:
        return TracedHirschbergReport(a, b, None)

    a_split_index = len(a) // 2
    a_l, a_r = a[:a_split_index], a[a_split_index:]

    dp_l = calculate_distance_matrix(a_l, b, config)
    dp_r = calculate_distance_matrix(a_r[::-1], b[::-1], config)[::-1]
    common_dp = [d_l + d_r for d_l, d_r in zip(dp_l, dp_r)]
    b_split_index = max(range(len(common_dp)), key=common_dp.__getitem__)
    b_l, b_r = b[:b_split_index], b[b_split_index:]

    left_children = traced_hirschberg(a_l, b_l, config)
    right_children = traced_hirschberg(a_r, b_r, config)
    return TracedHirschbergReport(a, b, (left_children, right_children))


def hirschberg_call_tree(a: str, b: str, del_cost: int, ins_cost: int, match_cost: int, mismatch_cost: int):
    config = HirschbergConfig(del_cost, ins_cost, match_cost, mismatch_cost)
    report = traced_hirschberg(StrHolder(a), StrHolder(b), config)
    report.save('hirschberg')


if __name__ == '__main__':
    hirschberg_call_tree(
        "AGTACGCA",
        "TATGC",
        del_cost=-2,
        ins_cost=-2,
        match_cost=2,
        mismatch_cost=-1
    )
