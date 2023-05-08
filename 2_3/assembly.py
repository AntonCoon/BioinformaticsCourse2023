import os
import subprocess
from collections import defaultdict
from typing import List, Dict

import graphviz
import matplotlib.pyplot as plt
import numpy as np


class PairedStr:
    def __init__(self, first: str, second: str):
        self.first = first
        self.second = second
        self.h = hash(f"{self.first},{self.second}")
        assert len(self.first) == len(self.second)

    def __len__(self):
        return len(self.first)

    def __add__(self, other):
        return PairedStr(self.first + other.first, self.second + other.second)

    def substr(self, i, j):
        return PairedStr(self.first[i:j], self.second[i:j])

    def __str__(self):
        return f"{self.first},{self.second}"

    def __hash__(self):
        return self.h

    def __eq__(self, other):
        return self.first == other.first and self.second == other.second


class PairedStrPointer:
    def __init__(self, paired_str: PairedStr, start: int, length: int):
        self.paired_str = paired_str
        self.start = start
        self.length = length

    def shift(self, x: int):
        self.start += x

    def is_valid_start(self) -> bool:
        return self.start + self.length <= len(self.paired_str)

    def prefix(self):
        return PairedStrPointer(self.paired_str, self.start, self.length - 1)

    def suffix(self):
        return PairedStrPointer(self.paired_str, self.start + 1, self.length - 1)

    def last_symbol(self):
        return PairedStrPointer(self.paired_str, self.start + self.length - 1, 1)

    def unpack(self):
        return self.paired_str.substr(self.start, self.start + self.length)


class GraphNode:
    def __init__(self, id: int, paired_str: PairedStr):
        self.id = id
        self.edges: Dict[PairedStr, GraphEdge] = {}
        self.paired_str = paired_str
        self.in_edge_count = 0

    def add_edge(self, edge):
        if edge.paired_str in self.edges:
            self.edges[edge.paired_str].weight += edge.weight
        else:
            self.edges[edge.paired_str] = edge
            edge.dest.in_edge_count += 1

    def remove_edge(self, label):
        edge = self.edges[label]
        edge.dest.in_edge_count -= 1
        self.edges.pop(label)

    # def __eq__(self, other)  # Comparing by addresses is good


class GraphEdge:
    def __init__(self, source: GraphNode, dest: GraphNode, paired_str: PairedStr, weight: float):
        self.source = source
        self.dest = dest
        self.paired_str = paired_str
        self.weight = weight

    # For using heap of edges
    def __lt__(self, other):
        return self.weight > other.weight


class Graph:
    def __init__(self):
        self.vertexes: Dict[PairedStr, GraphNode] = {}

    def get_or_create_vertex(self, paired_str: PairedStr) -> GraphNode:
        if paired_str not in self.vertexes:
            idx = len(self.vertexes)
            result = GraphNode(idx, paired_str)
            self.vertexes[paired_str] = result
        else:
            result = self.vertexes[paired_str]
        return result

    def add_edge(self, text: PairedStrPointer, weight: int):
        prefix_node = self.get_or_create_vertex(text.prefix().unpack())
        suffix_node = self.get_or_create_vertex(text.suffix().unpack())
        edge = GraphEdge(prefix_node, suffix_node, text.last_symbol().unpack(), weight)
        prefix_node.add_edge(edge)

    def remove_node(self, node: GraphNode):
        self.vertexes.pop(node.paired_str)

    def clear_isolated_vertexes(self):
        self.vertexes = {k: v for k, v in self.vertexes.items() if v.in_edge_count > 0 or len(v.edges) > 0}

    def count_vertexes(self) -> int:
        return len(self.vertexes)

    def count_edges(self) -> int:
        return sum(node.in_edge_count for node in self.vertexes.values())

    def count_sources(self) -> int:
        return sum(1 for node in self.vertexes.values() if node.in_edge_count == 0)

    def count_sinks(self) -> int:
        return sum(1 for node in self.vertexes.values() if len(node.edges) == 0)

    def count_components(self) -> int:
        used: Dict[GraphNode, bool] = defaultdict(bool)
        rev_edges: Dict[GraphNode, List[GraphEdge]] = defaultdict(list)
        for node in self.vertexes.values():
            for edge in node.edges.values():
                rev_edges[edge.dest].append(edge)

        def dfs(v: GraphNode):
            if used[v]:
                return
            used[v] = True
            for edge in v.edges.values():
                dfs(edge.dest)
            for rev_edge in rev_edges[v]:
                dfs(rev_edge.source)

        res = 0
        for node in self.vertexes.values():
            if not used[node]:
                res += 1
                dfs(node)
        return res

    def info(self) -> str:
        return f"Graph metrics:\n" \
               f"|V| = {self.count_vertexes()};\n" \
               f"|E| = {self.count_edges()}\n" \
               f"Components count = {self.count_components()}\n" \
               f"Sources count = {self.count_sources()}\n" \
               f"Sinks count = {self.count_sinks()}"

    def plot_degrees(self, name):
        def plot_histogram(values, axis, xlabel, ylabel, xlim=None):
            counts, bins = np.histogram(values, bins=25,
                                        range=(0, xlim) if xlim is not None else None)
            axis.stairs(counts, bins)
            axis.set_xlabel(xlabel)
            axis.set_ylabel(ylabel)

        v_degrees = [len(node.edges) for node in self.vertexes.values()]
        e_lengths = [len(edge.paired_str) for node in self.vertexes.values() for edge in node.edges.values()]
        e_weights = [edge.weight for node in self.vertexes.values() for edge in node.edges.values()]
        figure, axis = plt.subplots(2, 2)
        plot_histogram(v_degrees, axis[0, 0], "Степени вершин", "Количество вершин")
        plot_histogram(e_lengths, axis[1, 0], "Длина ребра", "Количество ребер")
        plot_histogram(e_weights, axis[0, 1], "Вес ребер", "Количество ребер", 10.0)
        figure.tight_layout()
        plt.savefig(name)
        plt.show()

    def get_dot(self, name=None) -> graphviz.Digraph:
        name = name or 'assembly'
        dot = graphviz.Digraph(name, comment='Genome assembly graph')
        for node in self.vertexes.values():
            cur_str = str(node.paired_str)
            dot.node(f"V_{cur_str}", f"{cur_str}")
        for node in self.vertexes.values():
            for edge in node.edges.values():
                dot.edge(f"V_{str(node.paired_str)}",
                         f"V_{str(edge.dest.paired_str)}",
                         f"{round(edge.weight, 2)} {str(edge.paired_str)}")
        return dot

    def view(self, name=None):
        dot = self.get_dot(name)
        print(f'Visualization {name}')
        dot.view()

    def save(self, filename):
        dot = self.get_dot()
        print(f'Dumping to file')
        dot.render(filename)


def read_universal(file_1: str, file_2: str, pair_mode: bool, step: int) -> List[PairedStr]:
    with open(file_1, 'r') as file:
        lines_1 = file.readlines()
    with open(file_2, 'r') as file:
        lines_2 = file.readlines()
    result = []
    if pair_mode:
        assert len(lines_1) == len(lines_2)
        for i in range(0, len(lines_1), step):
            if len(lines_1[i]) == 0:
                continue
            result.append(PairedStr(lines_1[i + 1], lines_2[i + 1]))
    else:
        for i in range(0, len(lines_1), step):
            result.append(PairedStr(lines_1[i + 1], lines_1[i + 1]))
        for i in range(0, len(lines_2), step):
            result.append(PairedStr(lines_2[i + 1], lines_2[i + 1]))
    return result


def read_fasta(file_1: str, file_2: str, pair_mode: bool) -> List[PairedStr]:
    return read_universal(file_1, file_2, pair_mode, 2)


def read_fastq(file_1: str, file_2: str, pair_mode: bool) -> List[PairedStr]:
    return read_universal(file_1, file_2, pair_mode, 4)


def build_graph(paired_reads: List[PairedStr], k: int) -> Graph:
    graph = Graph()
    for i, paired_read in enumerate(paired_reads):
        pointer = PairedStrPointer(paired_read, 0, k + 1)
        while pointer.is_valid_start():
            # if not pointer.contains_trash():
            graph.add_edge(pointer, 1)
            pointer.shift(1)
        if i % 20_000 == 0:
            print(f"Iteration {i} is finished")
    return graph


def compress_graph(graph: Graph):
    removed_count = 0
    for start_node in graph.vertexes.values():
        edges = list(start_node.edges.values())
        for start_edge in edges:
            middle_node = start_edge.dest
            while start_node is not middle_node and len(middle_node.edges) == 1 and middle_node.in_edge_count == 1:
                end_edge = list(middle_node.edges.values())[0]
                end_node = end_edge.dest
                a, b = len(start_edge.paired_str), len(end_edge.paired_str)
                new_weight = (a * start_edge.weight + b * end_edge.weight) / (a + b)
                new_edge = GraphEdge(start_node, end_node, start_edge.paired_str + end_edge.paired_str, new_weight)
                start_node.remove_edge(start_edge.paired_str)
                middle_node.remove_edge(end_edge.paired_str)
                start_node.add_edge(new_edge)
                removed_count += 1
                middle_node = end_node
                start_edge = new_edge
    graph.clear_isolated_vertexes()
    print(f"Compression: Removed {removed_count} edges")


def remove_tips(graph: Graph):
    tips = []
    for start_node in graph.vertexes.values():
        edges = list(start_node.edges.values())
        for edge in edges:
            end_node = edge.dest
            if start_node is end_node:
                continue  # Edge to itself
            if len(end_node.edges) == 0:
                tips.append(edge)
    tips.sort(key=lambda x: x.weight)
    will_be_deleted = int(len(tips) * 0.3)
    for edge in tips[:will_be_deleted]:
        edge.source.remove_edge(edge.paired_str)
    graph.clear_isolated_vertexes()
    print(f"Tips: Removed {will_be_deleted} edges")


def remove_single_edges(graph: Graph):
    removed_count = 0
    for start_node in graph.vertexes.values():
        edges = list(start_node.edges.values())
        for edge in edges:
            if edge.weight < 1 + 1e-8:
                start_node.remove_edge(edge.paired_str)
                removed_count += 1
    graph.clear_isolated_vertexes()
    print(f"Single edges: Removed {removed_count} edges")


def remove_bubbles(graph: Graph, k: int):
    removed_count = 0
    for start_node in graph.vertexes.values():
        edges_grouped_by_dest: Dict[GraphNode, List[GraphEdge]] = defaultdict(list)
        for edge in start_node.edges.values():
            edges_grouped_by_dest[edge.dest].append(edge)
        for dest, edges in edges_grouped_by_dest.items():
            edges.sort(key=lambda x: x.weight)
            for edge in edges[:-1]:  # At least one edge must be left
                if edge.weight < 1 + 1e-5 and len(edge.paired_str) <= k:
                    start_node.remove_edge(edge.paired_str)
                    removed_count += 1
    graph.clear_isolated_vertexes()
    print(f"Bubbles: Removed {removed_count} edges")


def main(file_1: str, file_2: str, k: int, pair_mode: bool, file_preparation: bool):
    print(f"K={k}, pair_mode={pair_mode}, file_preparation={file_preparation}")
    if file_preparation:
        _, file_1_name = os.path.split(file_1)
        _, file_2_name = os.path.split(file_2)
        f_1_unpaired = f"{file_1_name}_unpaired.fastq"
        f_1_paired = f"{file_1_name}_paired.fastq"
        f_2_unpaired = f"{file_2_name}_unpaired.fastq"
        f_2_paired = f"{file_2_name}_paired.fastq"
        subprocess.run(
            f"java -jar trimmomatic-0.39.jar PE -phred33 {file_1} {file_2} "
            f"{f_1_paired} {f_1_unpaired} {f_2_paired} {f_2_unpaired} "
            f"SLIDINGWINDOW:{k}:37 MINLEN:100".split(" ")
        )  # Drop low quality sequences
        file_1 = f_1_paired
        file_2 = f_2_paired
        subprocess.run(
            f"./pollux -p -i {file_1} {file_2} -o pollux_out -k {k} -n false -d false -h false".split(" ")
        )  # Repair sequences
        file_1 = f"pollux_out/{file_1}.corrected"
        file_2 = f"pollux_out/{file_2}.corrected"
    if file_1.endswith('.fasta') and file_2.endswith('.fasta'):
        reads = read_fasta(file_1, file_2, pair_mode)
    elif (file_1.endswith('.fastq') or file_1.endswith('.fastq.corrected')) and \
            (file_2.endswith('.fastq') or file_2.endswith('.fastq.corrected')):
        reads = read_fastq(file_1, file_2, pair_mode)
    else:
        raise RuntimeError('Unrecognized file format')
    print(f"Reads were scanned: {len(reads)}")
    print("Graph is building")
    graph = build_graph(reads, k)
    print(f"Graph was built")
    print(graph.info())
    print()
    graph.plot_degrees('begin-stats.png')

    print("Graph is compressing")
    compress_graph(graph)
    print(f"Graph is compressed")
    print(graph.info())
    print()

    print("Tips are removing")
    remove_tips(graph)
    compress_graph(graph)
    print(f"Tails are removed")
    print(graph.info())
    print()

    print("Bubbles are removing")
    remove_bubbles(graph, k)
    compress_graph(graph)
    print(f"Bubbles are removed")
    print(graph.info())
    print()

    graph.plot_degrees('middle-stats.png')

    print("Single edges are removing")
    remove_single_edges(graph)
    compress_graph(graph)
    print(f"Single edges are removed")
    print(graph.info())
    print()
    graph.plot_degrees('final-stats.png')
    graph.save('result-assembly')


if __name__ == '__main__':
    file_1 = "EAS20_8/s_6_1.fastq"
    file_2 = "EAS20_8/s_6_2.fastq"

    k = int(input("Enter k: "))  # Vertexes <= 4^(2k), Edges <= 4^(2k + 2)
    preparation = True
    pair_mode = True
    main(file_1, file_2, k, pair_mode, preparation)

