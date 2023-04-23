# well, we build graph
def build_overlap_graph(S: list, t: int) -> dict:
    graph = {}
    for i in range(len(S)):
        graph[S[i]] = []
        for j in range(len(S)):
            if i != j:
                suffix = get_longest_suffix(S[i], S[j], t)
                if suffix:
                    graph[S[i]].append((S[j], len(suffix)))
    return graph


# let`s unite graph nodes
def compress_graph(graph: dict, t: int) -> list:
    while len(graph) > 1:
        nodes_list = list(graph)
        new_node = 0
        neighbor = 0
        new_graph = None
        for node in graph:
            neighbor = get_max_weight_edge(graph, node)
            if neighbor != None:
                new_node = merge_vertices(node, neighbor) 
                del nodes_list[0:2]
                nodes_list.insert(0, new_node)
                new_graph = build_overlap_graph(nodes_list, t)
                break
            else:
                for next in graph:
                    if node == next:
                        continue
                    else:
                        new_node = merge_vertices(node, next)
                        del nodes_list[0:2]
                        nodes_list.insert(0, new_node)
                        new_graph = build_overlap_graph(nodes_list, t)
                        if len(new_graph) == 1:
                            return list(new_graph)
                        break
        graph = new_graph
    return list(graph)
            

# name of this function says all for itself
def get_longest_suffix(s1: str, s2: str, threshold: int) -> str:
    for i in range(len(s2), 0, -1):
        suffix = s2[:i]
        if s1.endswith(suffix) and len(suffix) >= threshold:
            return suffix
    return ""


# make 1 node from 2 separate, a.e. 'AAA' and 'AAT' -> 'AAAT'
def merge_vertices(vertex1: str, vertex2: str) -> str:
    common_suffix = ""
    for i in range(1, min(len(vertex1), len(vertex2))+1):
        if vertex1[-i:] == vertex2[:i]:
            common_suffix = vertex1[-i:]
    return vertex1 + vertex2[len(common_suffix):]
    

# find the biggest neighbor
def get_max_weight_edge(graph: dict, vertex: str) -> str:
    max_weight = 0
    max_weight_neighbor = None
    for neighbor, weight in graph[vertex]:
        if weight > max_weight:
            max_weight = weight
            max_weight_neighbor = neighbor
    return max_weight_neighbor


subseq = ['AAA', 'AAT', 'ATT', 'TTA', 'TTT']
threshold = 1

G = build_overlap_graph(subseq, threshold)
print(compress_graph(G, threshold))

