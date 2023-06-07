import numpy as np
def viterby(x, pi, probs, probs_in, probs_out):
    T = [[0 for i in range(len(x))] for j in range(len(pi))]
    P = [[0 for i in range(len(x))] for j in range(len(pi))]
    for pi_i in range(len(pi)):
        T[pi_i][0] = probs[pi_i] * probs_out[pi[pi_i]][alphabet[x[0]]]

    for x_i in range(1, len(x)):
        for pi_i in range(len(pi)):
            k = np.argmax([(T[k][x_i-1] * probs_in[pi[k]][pi[pi_i]] * probs_out[pi[pi_i]][alphabet[x[x_i]]]) for k in range(len(pi))])
            T[pi_i][x_i] = T[k][x_i-1] * probs_in[pi[k]][pi[pi_i]] * probs_out[pi[pi_i]][alphabet[x[x_i]]]
            P[pi_i][x_i] = k
    best_path = []
    k = np.argmax([T[k][len(x)-1] for k in range(len(pi))])
    for x_i in range(len(x) -1, -1, -1):
        best_path.insert(0, pi[k])
        k = P[k][x_i]
    return best_path




x = 'zyxxyzyzxyxyxzzxzzxxxyyxzzxyyzxyzxzyzyxxzxxzyyxzzyxzxxyzxzzzzxyzyxzyzxxxyyyyzzyzxyxxzzzzxzxzyyzxyzyx'
alphabet = ['x', 'y', 'z']
states = ['A', 'B']
x_num = [alphabet.index(i) for i in x]
transition = {'A': {'A': 0.654, 'B': 0.346}, 'B': {'A': 0.687, 'B': 0.313}}
emission = {'A': {'x': 0.182, 'y': 0.092, 'z': 0.726}, 'B': {'x': 0.448, 'y': 0.449, 'z': 0.103}}

initial_probs = [1/len(states) for _ in states]

path = viterby(x_num, states, initial_probs, transition, emission)

print(''.join(path))

