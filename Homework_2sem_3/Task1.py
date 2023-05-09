# 1. Построение графа Де Брюина (5 баллов)
# По заданному набору ридов в формате FASTQ и параметру k, который соответствует длине k-меров, построить
# граф Де Брюина, некоторый путь в котором соответствовал бы возможной подстроке в исходном геноме.
# Не забывайте про запоминание покрытия каждого k-мера, а так же про сами подстроки, которые соответствуют
# каждому ребру. В остальном граф полностью соответствует тому, что был описан в лекции.



def read_reads(file):
    f = open(file, 'r')
    lines = f.readlines()
    f.close()
    reads = []

    for line in lines:
        if line[0] == 'A' or line[0] == 'C' or line[0] == 'T' or line[0] == 'G':
            reads = reads + [line.rstrip()]

    return reads

def construct_graph(reads):
    nodes = []
    edges = []
    for r in reads:
        r1 = r[:-1]
        r2 = r[1:]
        nodes.append(r1)
        nodes.append(r2)
        edges.append([r1,r2])
    return (nodes, edges)

def make_node_edge_map(edges):
    node_edge_map = {}
    for e in edges:
        n = e[0]
        if n in node_edge_map:
            node_edge_map[n].append(e[1])
        else:
            node_edge_map[n] = [e[1]]
        if e[1] not in node_edge_map:
            node_edge_map[e[1]] = []
    return node_edge_map

reads = read_reads('my_seq.fastq')

# идеальное покрытие (из файла)
g = construct_graph(reads)
#print(g)

m = make_node_edge_map(g[1])
print(m)

# повторы (примеры ридов из лекции)
reads1 = ["AAAAA", "AAAAA", "AAAAC", 'AAACA', 'ACATG', 'CATGG', 'ATGGG', 'TGGGA', 'GGGAT','GATGT', 'ATGTT']
g1 = construct_graph(reads1)
#print(g1)

m1 = make_node_edge_map(g1[1])
print(m1)
