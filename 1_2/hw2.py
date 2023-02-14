from collections import defaultdict
from typing import Dict, Tuple
import Bio.Align


class DPCell:
    def __init__(self):
        self._initialized: bool = False
        self._value: int = 0
        self._parent: Tuple = ()

    def update(self, value: int, parent: Tuple):
        if not self.is_init() or self._value < value:
            self._initialized = True
            self._value = value
            self._parent = parent

    def is_init(self):
        return self._initialized

    def value(self):
        assert self.is_init()
        return self._value

    def parent(self):
        assert self.is_init()
        return self._parent


class DPField:
    def __init__(self, dp_shape):
        def create_mdim_list(shape):
            if len(shape) == 0:
                return DPCell()
            else:
                return [create_mdim_list(shape[1:]) for _ in range(shape[0])]
        self.field = create_mdim_list(dp_shape)

    def __getitem__(self, coords):
        return self.get_cell(coords)

    def get_cell(self, coords: Tuple):
        def extract_cell(current_field, rec_coords):
            if len(rec_coords) == 0:
                return current_field
            else:
                return extract_cell(current_field[rec_coords[0]], rec_coords[1:])
        return extract_cell(self.field, coords)

    def try_update(self, source_coords: Tuple, diff: int, target_coords: Tuple):
        source_cell = self.get_cell(source_coords)
        target_coords = self.get_cell(target_coords)
        if source_cell.is_init():
            target_coords.update(source_cell.value() + diff, source_coords)


def needleman_wunsch(a: str, b: str, replacement_matrix: Dict[str, Dict[str, int]],
                     gap_price: int) -> Tuple[str, str, int]:
    n, m = len(a), len(b)
    dp = DPField((n + 1, m + 1))
    dp[0, 0].update(0, None)
    for i in range(n + 1):
        for j in range(m + 1):
            if i > 0:
                dp.try_update((i - 1, j), gap_price, (i, j))
            if j > 0:
                dp.try_update((i, j - 1), gap_price, (i, j))
            if i > 0 and j > 0:
                dp.try_update((i - 1, j - 1), replacement_matrix[a[i - 1]][b[j - 1]], (i, j))

    finish_coords = (n, m)
    weight = dp[finish_coords].value()

    a_star_rev = ""
    b_star_rev = ""
    current_coords = finish_coords
    while current_coords != (0, 0):
        i, j = current_coords
        ni, nj = dp[i, j].parent()
        a_star_rev += "_" if i == ni else a[i - 1]
        b_star_rev += "_" if j == nj else b[j - 1]
        current_coords = (ni, nj)

    a_star = a_star_rev[::-1]
    b_star = b_star_rev[::-1]

    return a_star, b_star, weight


def affine_gaps(a: str, b: str, replacement_matrix: Dict[str, Dict[str, int]],
                gap_start_price: int, gap_continue_price: int) -> Tuple[str, str, int]:
    n, m = len(a), len(b)
    # Last dimension - bit mask of used gaps
    # 0 - no gaps at the ends, 1 - gap at the end of `a`, 2 - gap at the end of `b`
    dp = DPField((n + 1, m + 1, 3))
    dp[0, 0, 0].update(0, None)
    for i in range(n + 1):
        for j in range(m + 1):
            if i > 0:
                dp.try_update((i - 1, j, 0), gap_start_price, (i, j, 2))  # Start gap in `b`
                dp.try_update((i - 1, j, 1), gap_start_price, (i, j, 2))  # Start gap in `b`
                dp.try_update((i - 1, j, 2), gap_continue_price, (i, j, 2))  # Continue gap in `b`

            if j > 0:
                dp.try_update((i, j - 1, 0), gap_start_price, (i, j, 1))  # Start gap in `a`
                dp.try_update((i, j - 1, 1), gap_continue_price, (i, j, 1))  # Continue gap in `a`
                dp.try_update((i, j - 1, 2), gap_start_price, (i, j, 1))  # Start gap in `a`

            if i > 0 and j > 0:
                dp.try_update((i - 1, j - 1, 0), replacement_matrix[a[i - 1]][b[j - 1]], (i, j, 0))
                dp.try_update((i - 1, j - 1, 1), replacement_matrix[a[i - 1]][b[j - 1]], (i, j, 0))
                dp.try_update((i - 1, j - 1, 2), replacement_matrix[a[i - 1]][b[j - 1]], (i, j, 0))

    finish_coords = (n, m, 0)
    for k, possible_finish_cell in enumerate(dp[n, m]):
        if possible_finish_cell.value() > dp[finish_coords].value():
            finish_coords = (n, m, k)
    weight = dp[finish_coords].value()

    a_star_rev = ""
    b_star_rev = ""
    current_coords = finish_coords
    while current_coords != (0, 0, 0):
        i, j, k = current_coords
        ni, nj, nk = dp[i, j, k].parent()
        a_star_rev += "_" if i == ni else a[i - 1]
        b_star_rev += "_" if j == nj else b[j - 1]
        current_coords = (ni, nj, nk)

    a_star = a_star_rev[::-1]
    b_star = b_star_rev[::-1]

    return a_star, b_star, weight


if __name__ == '__main__':
    raw_matrix = Bio.Align.substitution_matrices.load('BLOSUM62')
    matrix = defaultdict(dict)
    for (k1, k2), v in raw_matrix.items():
        matrix[k1][k2] = v

    print(needleman_wunsch("AGTA", "ATA", matrix, -4))
    print(affine_gaps("AGTA", "ATA", matrix, -4, 4))
