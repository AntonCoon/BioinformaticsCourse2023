import numpy as np
from Bio import SeqIO
import collections
import subprocess


# Читаем файл с ридами 
reads = []
with open("D:/PYTHON/ALGORITHMS_IN_BIOINFORMATICS/10HW/reads.fastq", "r") as handle:
    for record in SeqIO.parse(handle, "fastq"):
        reads.append(record.seq)
        
# Обработка ридов Trimmomatic
# subprocess.run(["java", "-jar", "trimmomatic.jar", "SE", "D:/PYTHON/ALGORITHMS_IN_BIOINFORMATICS/10HW/reads.fastq", "trimmed_reads.fastq", "ILLUMINACLIP:TruSeq3-SE.fa:2:30:10", "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"])
        
# ПАРАМЕТР ДЛИНЫ К-МЕРОВ
k = 3

# Функция разбиения рида на k-меры
def read_to_kmers(read: str, k: int) -> list:
    kmers = []
    for i in range(len(read) - k + 1):
        kmers.append(read[i:i+k])
    return kmers

# Преобразуем риды в набор k-меров
kmers = []
for read in reads:
    kmers.extend(read_to_kmers(str(read), k))
    
def get_kmer_and_prev_kmer_sets(reads: list, k: int) -> set:
    kmer_set = set()
    prev_kmer_set = set()
    for read in reads:
        kmers = read_to_kmers(read, k)
        kmer_set.update(kmers)
        prev_kmer_set.update(kmer[:-1] for kmer in kmers)
    return kmer_set, prev_kmer_set

kmer_set, kmer_minus_one_set = get_kmer_and_prev_kmer_sets(reads, k)
coverage = {kmer: 1 for kmer in kmer_set} # Начальное покрытие всех ребер равно 1
coverage = dict(coverage)
graph_dict = {kmer_minus_one: [] for kmer_minus_one in kmer_minus_one_set}

# Создаем граф Де Брюина
for kmer in kmer_set:
    kmer_minus_one = kmer[:-1]
    if kmer_minus_one in graph_dict:
        graph_dict[kmer_minus_one].append(kmer)

# сжимаем граф
def compress_graph(graph_dict: dict, kmer_set: set) -> dict:
    compressed_graph = {}
    for kmer_minus_one in graph_dict.keys():
        neighbors = graph_dict[kmer_minus_one]
        if len(neighbors) == 1: # Если у вершины только один входящий и один исходящий сосед
            neighbor = neighbors[0]
            merged_kmer = kmer_minus_one + neighbor[-1] # Объединяем две вершины
            if len(kmer_minus_one) + len(neighbor[-1:]) == len(merged_kmer): # Если длины совпадают, то мы можем объединить ребра
                compressed_graph[kmer_minus_one] = [merged_kmer] # Обновляем граф
                kmer_set.remove(kmer_minus_one) # Удаляем старые вершины из множества k-меров
                kmer_set.remove(neighbor)
                # Обновляем покрытие ребер, пересчитывая взвешенное среднее
                weight_sum = len(kmer_minus_one) + len(neighbor[-1:])
                compressed_coverage = (len(kmer_minus_one) / weight_sum) * coverage[kmer_minus_one] + (len(neighbor[-1:]) / weight_sum) * coverage[neighbor]
                coverage[merged_kmer] = compressed_coverage
        else:
            compressed_graph[kmer_minus_one] = neighbors # Добавляем вершину без изменений
    return compressed_graph, kmer_set

# удаляем хвосты
def remove_low_coverage_tails(graph, coverage, threshold=0.3):
    ggraph = graph
    tail_coverage = {}  # собираем информацию о покрытиях концов хвостов

    for node in graph.keys():
        for neighbor in graph[node]:
            # получаем хвосты
            node_tail = node[1:]
            neighbor_tail = neighbor[1:]

            # проверяем покрытие и запоминаем
            for elem in coverage:
                node_tail_coverage = node_tail in elem
                neighbor_tail_coverage = neighbor_tail in elem
                tail_coverage[node_tail] = tail_coverage.get(node_tail, []) + [node_tail_coverage]
                tail_coverage[neighbor_tail] = tail_coverage.get(neighbor_tail, []) + [neighbor_tail_coverage]

    # определяем порог покрытия для хвостов
    cutoff = sorted(val for key, val in tail_coverage.items() if len(val) > 1)[int(threshold * len(tail_coverage))]

    # удаляем низкопокрытые хвосты
    for node in graph:
        new_neighbors = []
        for neighbor in graph[node]:
            # получаем хвосты
            node_tail = node[1:]
            neighbor_tail = neighbor[1:]

            # проверяем покрытие
            node_tail_coverage = node_tail in coverage
            neighbor_tail_coverage = neighbor_tail in coverage

            # сохраняем только хвосты с высоким покрытием
            if node_tail_coverage and neighbor_tail_coverage:
                tail_coverage.setdefault(node_tail, []).append(node_tail_coverage)
                tail_coverage.setdefault(neighbor_tail, []).append(neighbor_tail_coverage)
                new_neighbors.append(neighbor)

        graph[node] = new_neighbors

    return ggraph


# удаляем пузыри
def remove_bubbles(graph_dict, k, coverage):
    kmers_to_remove = set()
    for kmer in graph_dict.keys():
        if len(graph_dict[kmer]) == 1: # Если у вершины только один входящий и один исходящий сосед
            neighbor = graph_dict[kmer][0]
            path = [kmer, neighbor]
            while len(graph_dict[neighbor]) == 2: # Продолжаем поиск до тех пор, пока у соседа нет ветвлений
                next_neighbor = [n for n in graph_dict[neighbor] if n != kmer][0]
                path.append(next_neighbor)
                kmer, neighbor = neighbor, next_neighbor
            if len(path) <= 2 * k: # Если найден путь меньше или равный 2k, то удаляем вершины
                for p in path:
                    kmers_to_remove.add(p)
                # Обновляем граф
                start_node = path[0][:-1]
                end_node = path[-1][1:]
                compressed_kmer = start_node + end_node[-1]
                if len(start_node) + len(end_node[-1:]) == len(compressed_kmer):
                    graph_dict[start_node].remove(path[0])
                    graph_dict[end_node].remove(path[-1])
                    if compressed_kmer in graph_dict:
                        graph_dict[compressed_kmer].append(graph_dict[neighbor][0])
                    else:
                        graph_dict[compressed_kmer] = [graph_dict[neighbor][0]]
                    # Обновляем покрытие ребра
                    weight_sum = len(start_node) + len(end_node[-1:])
                    compressed_coverage = (len(start_node) / weight_sum) * coverage[start_node] + (len(end_node[-1:]) / weight_sum) * coverage[end_node]
                    coverage[compressed_kmer] = compressed_coverage
    # Удаляем к-меры из множества к-меров и из графа
    for kmer in kmers_to_remove:
        if kmer in graph_dict:
            del graph_dict[kmer]
        if kmer in coverage:
            del coverage[kmer]
    return graph_dict, coverage


compressed_graph, coverage = compress_graph(graph_dict, kmer_set)
compressed_graph = remove_low_coverage_tails(compressed_graph, coverage, 0.3)
compressed_graph = remove_bubbles(compressed_graph, k, coverage)

