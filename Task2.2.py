from typing import List, Tuple


def shortest_common_supersequence(strings: List[str], suffix_size: int) -> List[str]:
    overlap_graph = build_overlap_graph(strings, suffix_size)
    return condense_graph(overlap_graph, suffix_size)


def build_overlap_graph(strings: List[str], min_suffix_len: int) -> dict[str, List[Tuple[str, int]]]:
    graph = {string_i: [(string_j, len(get_longest_common_suffix(string_i, string_j, min_suffix_len)))
                        for j, string_j in enumerate(strings) if i != j and (suffix := get_longest_common_suffix(string_i, string_j, min_suffix_len))]
             for i, string_i in enumerate(strings)}
    return graph


def get_longest_common_suffix(string1: str, string2: str, min_suffix_len: int) -> str:
    return next((string2[:i] for i in range(len(string2), 0, -1) if string1.endswith(string2[:i]) and i >= min_suffix_len), "")


def get_max_weighted_neighbor(graph, node):
    return max(graph[node], key=lambda x: x[1])[0]


def merge_nodes(node1: str, node2: str) -> str:
    for i in range(min(len(node1), len(node2)), 0, -1):
        if node1[-i:] == node2[:i]:
            return node1 + node2[i:]
    return node1 + node2


def condense_graph(graph: dict[str, List[Tuple[str, int]]], threshold: int) -> List[str]:
    while len(graph) > 1:
        node, neighbor = next(((node, get_max_weighted_neighbor(graph, node)) for node in graph), (None, None))
        if neighbor is None:
            node, neighbor = next(((node, n_node) for node in graph for n_node in graph if node != n_node), (None, None))
        condensed_node = merge_nodes(node, neighbor)
        all_nodes = [condensed_node] + [n_node for n_node in graph if n_node not in (node, neighbor)]
        graph = build_overlap_graph(all_nodes, threshold)
        if len(graph) == 1:
            return list(graph)
    return list(graph)


strings = ['AAA', 'TTT', 'AAT', 'ATT', 'TTA']
suffix_size = 1

result = shortest_common_supersequence(strings, suffix_size)
print(f"Input strings: {', '.join(strings)}")
print(f"Suffix size: {suffix_size}")
print(f"Shortest common supersequence: {''.join(result)}")

#Output
#Input strings: AAA, TTT, AAT, ATT, TTA
#Suffix size: 1
#Shortest common supersequence: AAATTTA
