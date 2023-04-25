import heapq
from typing import List


# DSU with O(1) for operation
class DSU:
    def __init__(self):
        self.parent: List[int] = []
        self.sizes: List[int] = []

    def add_class(self):
        new_idx = len(self.parent)
        self.parent.append(new_idx)
        self.sizes.append(1)

    def is_merged(self, x: int, y: int) -> bool:
        return self.get_class(x) == self.get_class(y)

    def merge_classes(self, x: int, y: int) -> int:
        class_x = self.get_class(x)
        class_y = self.get_class(y)
        if class_x == class_y:
            return class_x
        if self.sizes[class_x] < self.sizes[class_y]:
            class_x, class_y = class_y, class_x
        self.parent[class_y] = class_x
        self.sizes[class_x] += self.sizes[class_y]
        return class_x

    def get_class(self, x: int):
        if self.parent[x] != x:
            self.parent[x] = self.get_class(self.parent[x])
        return self.parent[x]


class GraphEdge:
    def __init__(self, source: int, dest: int, weight: int):
        self.source = source
        self.dest = dest
        self.weight = weight

    # For using heap of edges
    def __lt__(self, other):
        return self.weight > other.weight


class GraphNode:
    def __init__(self, id: int, text: str):
        self.id = id
        self.text = text
        self.can_be_left = True
        self.can_be_right = True

    def cannot_be_left(self):
        self.can_be_left = False  # Vertex have right neighbour already

    def cannot_be_right(self):
        self.can_be_right = False  # Vertex have left neighbour already


class Graph:
    def __init__(self):
        self.nodes: List[GraphNode] = []
        # DSU keeps info about merged vertexes
        # Assembly of union is in the parent node
        self.dsu = DSU()
        # Heap of edges with the most large overlap at top
        self.edges_heap: List[GraphEdge] = []

    def create_node(self, s: str):
        idx = len(self.nodes)
        node = GraphNode(idx, s)
        self.nodes.append(node)
        self.dsu.add_class()
        print(f"Vertex {idx} with text <{s}> is created")

    def add_edge(self, source: int, dest: int, weight: int):
        edge = GraphEdge(source, dest, weight)
        heapq.heappush(self.edges_heap, edge)
        print(f"Edge {source} -> {dest} with weight {weight} is created")

    def get_active_nodes(self) -> List[GraphNode]:
        return [node for idx, node in enumerate(self.nodes) if self.dsu.get_class(idx) == idx]

    # O(1)
    def is_merged(self, x: int, y: int) -> bool:
        return self.dsu.is_merged(x, y)

    def merge_by_largest_edge(self) -> bool:
        largest_edge = None
        while len(self.edges_heap) > 0:
            cur_edge = heapq.heappop(self.edges_heap)
            left_node = self.nodes[cur_edge.source]
            right_node = self.nodes[cur_edge.dest]
            if not self.is_merged(left_node.id, right_node.id) and left_node.can_be_left and right_node.can_be_right:
                # Vertexes must be in different unions
                # Left node must have no right neighbour, right one must have no left neighbour
                largest_edge = cur_edge
                break
        if largest_edge is None:
            return False
        print(f"- Merge by edge: source={largest_edge.source}, dest={largest_edge.dest}, weight={largest_edge.weight}")
        self._merge_vertexes(largest_edge)
        return True

    # O(string concatenation)
    def _merge_vertexes(self, edge: GraphEdge):
        left_parent = self.nodes[self.dsu.get_class(edge.source)]
        right_parent = self.nodes[self.dsu.get_class(edge.dest)]
        common_text_length = edge.weight

        merged_node_idx = self.dsu.merge_classes(left_parent.id, right_parent.id)
        self.nodes[merged_node_idx].text = left_parent.text + right_parent.text[common_text_length:]
        print(f"- Text <{self.nodes[merged_node_idx].text}> saved to vertex {merged_node_idx}")

        self.nodes[edge.source].cannot_be_left()
        self.nodes[edge.dest].cannot_be_right()


def z_function(s) -> List[int]:
    s_len = len(s)
    res = [0] * s_len
    l = 0
    r = -1
    for i in range(1, s_len):
        current_sz = min(r - i + 1, res[i - l]) if i <= r else 0
        while i + current_sz < s_len and s[current_sz] == s[i + current_sz]:
            current_sz += 1
        res[i] = current_sz
        if current_sz > r - i + 1:
            l, r = i, i + current_sz - 1
    return res


def build_graph(reads: List[str], threshold: int) -> Graph:
    graph = Graph()
    for read in reads:
        graph.create_node(read)
    for i in range(len(reads)):
        for j in range(len(reads)):
            if i == j:
                continue
            z_str = reads[i] + "#" + reads[j]
            z_val = z_function(z_str)
            common_length = 0
            for k in range(1, min(len(reads[i]), len(reads[j]))):
                if z_val[-k] == k:
                    common_length = k
            if common_length >= threshold:
                graph.add_edge(j, i, common_length)
    return graph


def compute_assembly(reads: List[str], threshold: int) -> str:
    # O(n^2 * l)
    graph = build_graph(reads, threshold)
    # Summary O(E * (log E + string concatenation))
    while graph.merge_by_largest_edge():
        pass
    active_nodes = [node.text for node in graph.get_active_nodes()]
    print(f"Nodes {active_nodes} will be concatenated")
    concatenated_text = "".join(active_nodes)  # Usually single vertex is left
    return concatenated_text


def test(input: List[str], threshold: int):
    print(f"Test input: {input}")
    print(f"Test output: {compute_assembly(input, threshold)}")
    print()


if __name__ == '__main__':
    test(["AAA", "AAT", "ATT", "TTA", "TTT"], 1)
    test(["AAA", "CTA", "CCC", "TTA"], 1)
