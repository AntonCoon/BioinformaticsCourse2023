# https://rosalind.info/problems/ba10a/

def calculate_probability(line, trans):
    prob = 0.5 #initial probability
    for i in range(len(line) - 1):
        prob *= trans[line[i]][line[i+1]]
    return prob

path_pi = "AABBBAABABAAAABBBBAABBABABBBAABBAAAABABAABBABABBAB"
states = ['A', "B"]
transition = {'A':{'A': 0.194, 'B':0.806}, 'B':{'A': 0.273, 'B':0.727}}

print(calculate_probability(path_pi, transition))


