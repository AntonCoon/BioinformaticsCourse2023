from collections import defaultdict
import random
import numpy as np
from Bio import SeqIO
import subprocess

# 1 задание
def build_de_bruijn_graph(reads, k):
    graph = defaultdict(list)
    coverage = defaultdict(int)
    edges = {}
    for read in reads:
        kmers = (read[i:i + k] for i in range(len(read) - k + 1))
        for kmer in kmers:
            prefix = kmer[:-1]
            suffix = kmer[1:]
            graph[prefix].append(suffix)
            coverage[(prefix, suffix)] += 1
            edges[(prefix, suffix)] = kmer
    return graph, coverage, edges


# 2 задание
def compress_de_bruijn_graph(graph, coverage, edges):
    compressed_graph = {}
    compressed_coverage = {}
    compressed_edges = {}
    graph_items = graph.items()
    for node, neighbors in graph_items:
        if len(neighbors) == 1:
            one_neighbors = [k for k, v in graph_items if node in v]
            if len(one_neighbors) ==1:
                prev_node = one_neighbors[0]
                next_node = neighbors[0]
                compressed_graph.setdefault(prev_node, []).append(next_node)
                new_edge = edges[(prev_node, node)] + edges[(node, next_node)][-1]
                compressed_edges[(prev_node, next_node)] = new_edge
                new_coverage = (coverage[(prev_node, node)] * len(edges[(prev_node, node)]) + coverage[(node, next_node)] * len(edges[(node, next_node)])) / len(new_edge)
                compressed_coverage[(prev_node, next_node)] = new_coverage
            else:
                compressed_graph[node] = neighbors
        else:
            compressed_graph[node] = neighbors
        for neighbor in neighbors:
            compressed_coverage[(node, neighbor)] = coverage[(node, neighbor)]
            compressed_edges[(node, neighbor)] = edges[(node, neighbor)]
    return compressed_graph, compressed_coverage, compressed_edges


# 3 задание

def remove_low_coverage_tails(graph, coverage, edges):
    tail_coverages = []
    for node, neighbors in graph.items():
        if len(neighbors) == 1 and sum(1 for n in graph.values() if node in n) == 0:
            tail_coverage = coverage[(node, neighbors[0])] * len(edges[(node, neighbors[0])])
            tail_coverages.append(tail_coverage)
    if tail_coverages:
        threshold = np.percentile(tail_coverages, 30)
        for node, neighbors in list(graph.items()):
            if len(neighbors) == 1 and sum(1 for n in graph.values() if node in n) == 0:
                tail_coverage = coverage[(node, neighbors[0])] * len(edges[(node, neighbors[0])])
                if tail_coverage < threshold:
                    del graph[node]
                    del coverage[(node, neighbors[0])]
                    del edges[(node, neighbors[0])]
    return graph, coverage, edges


# 4 задача
def remove_bubbles(graph, coverage, edges, k):
    new_graph = graph.copy()
    for edge in edges:
        if coverage[edge] <= 2 * k:
            if edge[0] in new_graph and edge[1] in new_graph[edge[0]]:
                new_graph[edge[0]].remove(edge[1])
    return new_graph


def process_genom(read, k):
    with open('intermediate_calculations.txt', 'w') as f:
        graph, coverage, edges = build_de_bruijn_graph(read, k)
        f.write("graph: {}\n".format(graph))
        f.write("coverage: {}\n".format(coverage))
        f.write("edges: {}\n".format(edges))
        compressed_graph, compressed_coverage, compressed_edges = compress_de_bruijn_graph(graph, coverage, edges)
        f.write("compressed_graph: {}\n".format(compressed_graph))
        f.write("compressed_coverage: {}\n".format(compressed_coverage))
        f.write("compressed_edges: {}\n".format(compressed_edges))
        compressed_graph_without_tails, compressed_coverage_without_tails, compressed_edges_without_tails = remove_low_coverage_tails(
            compressed_graph, compressed_coverage, compressed_edges)
        f.write("compressed_graph_without_tails: {}\n".format(compressed_graph_without_tails))
        f.write("compressed_coverage_without_tails: {}\n".format(compressed_coverage_without_tails))
        f.write("compressed_edges_without_tails: {}\n".format(compressed_edges_without_tails))
        graph_without_bubbles = remove_bubbles(compressed_graph_without_tails, compressed_coverage_without_tails,
                                               compressed_edges_without_tails, k)
        f.write("graph_without_bubbles: {}\n".format(graph_without_bubbles))
    return graph_without_bubbles



# genome = "".join(random.choices(['A', 'C', 'G', 'T'], k=100000))
# reads = [genome[i:i + 150] for i in range(0, len(genome), 150)]
# print(process_genom(reads, 10))
sequences = []
path = 'D:\\Study\\BioinformaticsCourse2023\\2_1\\Trimmomatic-0.39\\trimmomatic-0.39.jar'
trim_command = f'java -jar {path} SE -phred33 ERR008613.fastq output.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 ' \
               f'TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'

result = subprocess.run(trim_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
with open("output.fastq", "r") as f, open("sequences.txt", "w") as out:
    for record in SeqIO.parse(f, "fastq"):
        sequences.append(record.seq)
        out.write("{}\n".format(record.seq))
        print("I am alive")
print("I am here!!!!!!!!!!!!!!!!!!")
print(process_genom(sequences, 20))
