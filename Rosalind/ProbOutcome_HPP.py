# https://rosalind.info/problems/ba10b/

def readFromFile(file_name):
    f = open(file_name, 'r')
    data = f.read().split('\n--------\n')
    x = data[0]
    alphabet = list(data[1].split())
    path = data[2]
    states = list(data[3].split())
    init_states = data[4].split()
    emission = {}
    corection = 1
    for i in range(len(states)):
        emission[states[i]] = {}
        for j in range(len(alphabet)):
            emission[states[i]][alphabet[j]] = float(init_states[len(alphabet)+i+j+corection])
        corection+=len(alphabet)
    f.close()
    return x, path, emission

def calculate_probability(x, path, emission):
    prob = 1
    for i in range(len(x)):
        prob *= emission[path[i]][x[i]]
    return prob

file_name = 'rosalind_ba10b.txt'
x, path, emission = readFromFile(file_name)
print(calculate_probability(x, path, emission))