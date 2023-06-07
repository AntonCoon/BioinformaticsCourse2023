def outcome_given_path(x, seq, emission_matrix):
    probability = 1
    for i in range(len(x)):
        emission_prob = emission_matrix[seq[i]][x[i]]
        probability *= emission_prob
    return probability


x = 'yxyyxzzxyyxyxyxyxxyxyxyyyxzxzxzyxzxzxzyyyzxzyzyyzz'
seq = 'AABBAABAAABAAAAABABBABBBABABAAAAAABAABABBABAAABAAA'
emission_matrix = {
    'A': {'x': 0.256, 'y': 0.372, 'z': 0.372},
    'B': {'x': 0.244, 'y': 0.293, 'z': 0.463}
}

prob = outcome_given_path(x, seq, emission_matrix)
print(prob)
