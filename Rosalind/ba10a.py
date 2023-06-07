def hidden_path_prob(seq, transition):
    states = list(transition.keys())
    prob = 1 / len(states)
    for i in range(len(seq) - 1):
        prob *= transition[seq[i]][seq[i + 1]]
    return prob

seq = 'AAABBBBBABABABABBABBAAAAAAABAAABABAAABABBBAAAABABA'
transition = {
    'A': {'A': 0.612, 'B': 0.388},
    'B': {'A': 0.258, 'B': 0.742}
}

prob = hidden_path_prob(seq, transition)
print(prob)
