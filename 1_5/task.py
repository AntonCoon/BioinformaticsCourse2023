import logging
from dataclasses import dataclass
from typing import List, Tuple

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger('multi-aligner')


@dataclass
class DPConfig:
    del_cost: int
    ins_cost: int
    match_cost: int
    mismatch_cost: int


def calculate_distance_matrix(a: str, b: str, config: DPConfig) -> List[int]:
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


def recursive_hirschberg(a: str, b: str, config: DPConfig) -> str:
    if len(a) <= 1 or len(b) <= 1:
        # If strings' length is equal to 1 then we choose letter from first one
        # Otherwise we take letter from nonempty string
        return a if len(a) == 1 and len(b) == 1 else a + b

    a_split_index = len(a) // 2
    a_l, a_r = a[:a_split_index], a[a_split_index:]

    dp_l = calculate_distance_matrix(a_l, b, config)
    dp_r = calculate_distance_matrix(a_r[::-1], b[::-1], config)[::-1]
    common_dp = [d_l + d_r for d_l, d_r in zip(dp_l, dp_r)]
    b_split_index = max(range(len(common_dp)), key=common_dp.__getitem__)
    b_l, b_r = b[:b_split_index], b[b_split_index:]

    left_children = recursive_hirschberg(a_l, b_l, config)
    right_children = recursive_hirschberg(a_r, b_r, config)
    return left_children + right_children


def hirschberg(a: str, b: str, config: DPConfig) -> Tuple[int, str]:
    optimal_distance = calculate_distance_matrix(a, b, config)[-1]
    consensus = recursive_hirschberg(a, b, config)
    return optimal_distance, consensus


class HirschbergCacheProxy:
    # Used for caching queries to hirschberg function
    def __init__(self, config: DPConfig):
        self.hirschberg_cache = {}
        self.config = config

    def get(self, a: str, b: str):
        key = a, b
        if key not in self.hirschberg_cache:
            logger.debug(f"New hirschberg calculation for {a} and {b}")
            self.hirschberg_cache[key] = hirschberg(a, b, self.config)
        return self.hirschberg_cache[key]


def greedy_multi_alignment(genomes: List[str],
                           del_cost: int, ins_cost: int, match_cost: int, mismatch_cost: int) -> str:
    genomes = genomes.copy()

    n = len(genomes)
    hirschberg_proxy = HirschbergCacheProxy(DPConfig(del_cost, ins_cost, match_cost, mismatch_cost))
    for _ in range(n - 1):
        best_pair, best_distance = None, None
        for i in range(len(genomes)):
            for j in range(i):
                distance, consensus = hirschberg_proxy.get(genomes[i], genomes[j])
                logger.debug(f"Distance between {genomes[i]} and {genomes[j]}: {distance}. Consensus: {consensus}")
                if best_distance is None or distance > best_distance:
                    best_pair = i, j
                    best_distance = distance
        i, j = best_pair
        _, best_consensus = hirschberg_proxy.get(genomes[i], genomes[j])
        logger.info(f"Nearest pair: {genomes[i]}, {genomes[j]}. Distance: {best_distance}. Consensus: {best_consensus}")
        genomes.pop(i)
        genomes.pop(j)
        genomes.append(best_consensus)
    assert len(genomes) == 1, "Only one genome must be in array at the end"
    return genomes[0]


def main():
    del_cost, ins_cost, match_cost, mismatch_cost = -1, -1, 1, -1
    genomes = ["GATTCA", "GTCTGA", "GATATT", "GTCAGC"]
    alignment = greedy_multi_alignment(genomes, del_cost, ins_cost, match_cost, mismatch_cost)
    print("Input genomes:", *genomes)
    print("Best alignment:", alignment)


if __name__ == '__main__':
    main()
