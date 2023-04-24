# 1. Жадный алгоритм поиска SCS (5 баллов)
# Напоминание: задача поиска SCS (shortest common supersequence) заключается в том, чтобы по заданному набору подстрок
# найти самую короткую строку, которая содержит заданные строки в качестве подстрок. Поиск точного решения - NP трудная
# задача, но существует жадный алгоритм, который находит приближенное решение. В данной домашней работе вам предстоит
# его реализовать.
# Сам алгоритм был описан в лекции 2.2 и заключается в том, чтобы построить граф перекрытий и итеративно сливать в нем
# те вершины, которые сильнее всего перекрываются. Построение графа можно реализовать наивно, за O(n^2*l^2))
# где n - число строк, а l это оценка сверху на их длину.

def built_overlap_graph(sss, t):

    def count_longest_suff(si, sj, t):
        for i in range(len(sj), 0, -1):
            suffix = sj[:i]
            if si.endswith(suffix) and len(suffix) >= t:
                return suffix
        return ""

    g = {}
    for i in range(len(sss)):
        g[sss[i]] = []
        for j in range(len(sss)):
            if i != j:
                suff = count_longest_suff(sss[i], sss[j], t)
                if suff:
                    g[sss[i]].append((sss[j], len(suff)))
    return g

def condensate_graph(g, t):

    def edge_w(g, node):
        max_w, max_n = 0, None
        for n, w in g[node]:
            if w > max_w:
                max_w, max_n = w, n
        return max_n

    def condensate_nodes(node1, node2):
        common_suff = ""
        length = len(node1) if len(node1) < len(node2) else len(node2)
        for i in range (1, length + 1):
            if node1[-i:] == node2[:i]:
                common_suff = node1[-i:]
        return node1 + node2[len(common_suff):]


    while len(g) > 1:
        all_nodes = list(g)
        conden_node, neighbor, conden_g = 0, 0, None
        for node in g:
            neighbor = edge_w(g, node)
            if neighbor != None:
                conden_node = condensate_nodes(node, neighbor)
                all_nodes.remove(node)
                all_nodes.remove(neighbor)
                all_nodes.insert(0, conden_node)
                conden_g = built_overlap_graph(all_nodes, t)
                break
            else:
                for n_node in g:
                    if node == n_node:
                        continue
                    else:
                        conden_node = condensate_nodes(node, n_node)
                        all_nodes.remove(node)
                        all_nodes.remove(n_node)
                        all_nodes.insert(0, conden_node)
                        conden_g = built_overlap_graph(all_nodes, t)
                        if len(conden_g) == 1:
                            return list(conden_g)
                        break
        g = conden_g
    return list(g)

def shortest_common_supersequence(sss, t):
    over_graph = built_overlap_graph(sss, t)
    return condensate_graph(over_graph, t)

ss1 = ['AAA', 'TTT', 'AAT', 'ATT', 'TTA']
ss1_2 = ['AAA', 'AAT', 'ATT', 'TTA', 'TTT']
ss2 = ['ACA', 'CGG', 'AAC', 'GTT', 'TTA', 'TGT', 'AA']
ss2_2 = ['ACA', 'AAC', 'TTA', 'TGT', 'GTT', 'AA', 'CGG']

suff_size = 1

print(*ss1)
print(*shortest_common_supersequence(ss1, suff_size))

print(*ss1_2)
print(*shortest_common_supersequence(ss1_2, suff_size))
print()
print(*ss2)
print(*shortest_common_supersequence(ss2, suff_size))
print(*ss2_2)
print(*shortest_common_supersequence(ss2_2, suff_size))

# Можно увидеть, что т.к. это жадный алгоритм, результат работы очень сильно зависит от порядка обхода элементов.
# Пример - набор из подпоследовательностей ss1 и ss1_2. При одних и тех же элементах результирующая строка разная
# Хотя для ss2 и ss2_2 строка при перестановке одинаковая. Возможно, чем больше различность подстрок,
# тем алгоритм работает лучше (выдает более выгодную строку), т.к. вариантов совмещения меньше
