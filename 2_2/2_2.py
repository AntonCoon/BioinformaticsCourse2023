import math

def greedy_scs(strings):
    graph = {}
    for i in range(len(strings)):
        graph[strings[i]] = []
        for j in range(len(strings)):
            if i != j:
                overlap_len = overlap(strings[i], strings[j])
                if overlap_len > 0:
                    graph[strings[i]].append((strings[j], overlap_len))

    while len(graph) > 1:
        max_overlap = -math.inf
        keys_to_remove = []
        for a in graph:
            for b, overlap_len in graph[a]:
                if overlap_len > max_overlap:
                    max_overlap = overlap_len
                    max_a, max_b = a, b
        keys_to_remove.append(max_a)
        keys_to_remove.append(max_b)
        for key in keys_to_remove:
            graph.pop(key, None)

        if max_a.endswith(max_b[:max_overlap]):
            new_string = max_a + max_b[max_overlap:]
        else:
            new_string = max_b + max_a[max_overlap:]
        graph[new_string] = []
        for a in list(graph.keys()):
            if a != new_string:
                overlap_len = overlap(new_string, a)
                if overlap_len > 0:
                    graph[new_string].append((a, overlap_len))
                    graph[a].append((new_string, overlap_len))

    return list(graph.keys())[0]


def overlap(a, b):
    max_overlap = 0
    for i in range(1, min(len(a), len(b))):
        if a[-i:] == b[:i]:
            max_overlap = i
    return max_overlap


strings = ["AAA", "AAT", "ATT", "TTA", "TTT"]
print(greedy_scs(strings))
