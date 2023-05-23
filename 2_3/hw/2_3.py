from collections import defaultdict
import random
import numpy as np
from Bio import SeqIO
import subprocess


# 1 задание


def read_fastq(file):
    with open(file, 'r') as f:
        while True:
            id_line = f.readline().strip()
            if not id_line:
                break
            seq_line = f.readline().strip()
            plus_line = f.readline().strip()
            quality_line = f.readline().strip()
            yield id_line, seq_line, plus_line, quality_line


def build_de_bruijn_graph(fastq_file, k):
    graph = defaultdict(set)
    coverage = defaultdict(int)
    edges = {}
    for id_line, seq_line, plus_line, quality_line in read_fastq(fastq_file):
        read = seq_line
        for i in range(len(read) - k + 1):
            kmer = read[i:i + k]
            prefix = kmer[:-1]
            suffix = kmer[1:]
            graph[prefix].add(suffix)
            if (prefix, suffix) not in edges:
                edges[(prefix, suffix)] = (kmer[0], kmer[-1])
            coverage[(prefix, suffix)] += 1
    return graph, coverage, edges


# 2 задание

def compress_de_bruijn_graph(graph, coverage, edges):
    compressed_graph = defaultdict(set)
    compressed_coverage = {}
    compressed_edges = {}
    graph_items = graph.items()
    for node, neighbors in graph_items:
        if len(neighbors) == 1:
            one_neighbors = [k for k, v in graph_items if node in v]
            if len(one_neighbors) == 1:
                prev_node = one_neighbors[0]
                next_node = list(neighbors)[0]
                compressed_graph[prev_node].add(next_node)
                a = edges[(prev_node, node)]
                b = edges[(node, next_node)][-1]

                if not isinstance(a, tuple):
                    a = tuple(a)
                if not isinstance(b, tuple):
                    b = tuple(b)

                new_edge = a + b
                # new_edge = edges[(prev_node, node)] + edges[(node, next_node)][-1]
                compressed_edges[(prev_node, next_node)] = new_edge
                new_coverage = (coverage[(prev_node, node)] * len(edges[(prev_node, node)]) + coverage[
                    (node, next_node)] * len(edges[(node, next_node)])) / len(new_edge)
                compressed_coverage[(prev_node, next_node)] = new_coverage
            else:
                compressed_graph[node] = neighbors
        else:
            compressed_graph[node] = neighbors
        for neighbor in neighbors:
            compressed_coverage[(node, neighbor)] = coverage[(node, neighbor)]
            compressed_edges[(node, neighbor)] = edges[(node, neighbor)]
    return dict(compressed_graph), compressed_coverage, compressed_edges


# 3 задание

def remove_low_coverage_tails(graph, coverage, edges):
    tail_coverages = []
    for node, neighbors in graph.items():
        if len(neighbors) == 1 and sum(1 for n in graph.values() if node in n) == 0:
            print(type(coverage), type(edges))
            print(type(neighbors))
            neighbors_list = list(neighbors)
            print(type(neighbors_list))
            tail_coverage = coverage[(node, neighbors_list[0])] * len(edges[(node, neighbors_list[0])])
            tail_coverages.append(tail_coverage)
    if tail_coverages:
        threshold = np.percentile(tail_coverages, 30)
        for node, neighbors in list(graph.items()):
            if len(neighbors) == 1 and sum(1 for n in graph.values() if node in n) == 0:
                print(type(coverage), type(edges))
                print(type(neighbors))
                neighbors_list = list(neighbors)
                print(type(neighbors_list))
                tail_coverage = coverage[(node, neighbors_list[0])] * len(edges[(node, neighbors_list[0])])
                if tail_coverage < threshold:
                    del graph[node]
                    del coverage[(node, neighbors_list[0])]
                    del edges[(node, neighbors_list[0])]
    return graph, coverage, edges


# 4 задача
def remove_bubbles(graph, coverage, edges, k):
    new_graph = graph.copy()
    new_coverage = coverage.copy()
    new_edges = edges.copy()

    for edge in edges:
        if coverage[edge] <= 2 * k:
            if edge[0] in new_graph and edge[1] in new_graph[edge[0]]:
                new_graph[edge[0]].remove(edge[1])
                del new_coverage[edge]
                del new_edges[edge]

    return new_graph, new_coverage, new_edges


def process_genom(read, k):
    with open('../intermediate_calculations.txt', 'w') as f:
        print("0")
        graph, coverage, edges = build_de_bruijn_graph(read, k)
        f.write("graph: {}\n".format(graph))
        f.write("coverage: {}\n".format(coverage))
        print("1")
        f.write("edges: {}\n".format(edges))
        compressed_graph, compressed_coverage, compressed_edges = compress_de_bruijn_graph(graph, coverage, edges)
        f.write("compressed_graph: {}\n".format(compressed_graph))
        f.write("compressed_coverage: {}\n".format(compressed_coverage))
        f.write("compressed_edges: {}\n".format(compressed_edges))
        print("2")
        compressed_graph_without_tails, compressed_coverage_without_tails, compressed_edges_without_tails = remove_low_coverage_tails(
            compressed_graph, compressed_coverage, compressed_edges)
        f.write("compressed_graph_without_tails: {}\n".format(compressed_graph_without_tails))
        f.write("compressed_coverage_without_tails: {}\n".format(compressed_coverage_without_tails))
        f.write("compressed_edges_without_tails: {}\n".format(compressed_edges_without_tails))
        graph_without_bubbles, coverage_without_bubbles, edges_without_bubbles = remove_bubbles(
            compressed_graph_without_tails, compressed_coverage_without_tails,
            compressed_edges_without_tails, k)
        print("3")
        f.write("graph_without_bubbles: {}\n".format(graph_without_bubbles))
        f.write("coverage_without_bubbles: {}\n".format(coverage_without_bubbles))
        f.write("edges_without_bubbles: {}\n".format(edges_without_bubbles))
    return graph_without_bubbles


# genome = "".join(random.choices(['A', 'C', 'G', 'T'], k=100000))
# reads = [genome[i:i + 150] for i in range(0, len(genome), 150)]
# print(process_genom(reads, 10))

print(process_genom("output.fastq", 10))
