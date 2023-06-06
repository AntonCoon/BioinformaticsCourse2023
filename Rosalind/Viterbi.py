# https://rosalind.info/problems/ba10c/

import sys
from math import *
import numpy as np

def readFromFile():
    f = open('rosalind_ba10c.txt', 'r')
    data = f.read().split()
    x = data[0]
    ind = [i for i in range(len(data)) if '--------' == data[i]]
    alphabet = data[ind[0] + 1:ind[1]]
    states = data[ind[1] + 1:ind[2]]
    stateDict = {i: states[i] for i in range(len(states))}
    transitionLog = {
        i: {k: log(float(data[ind[2] + len(states) + 2 + i * (len(states) + 1) + k])) for k in range(len(states))} for i
        in range(len(states))}
    emissionLog = {i: {alphabet[k]: log(float(data[ind[3] + len(alphabet) + 2 + i * (len(alphabet) + 1) + k])) for k in
                       range(len(alphabet))} for i in range(len(states))}
    f.close()
    return x, transitionLog, emissionLog, stateDict

def viterbi(x, transionLog, emissionLog, stateDict):
    n = len(x)
    l = len(transionLog)
    s = [[0 for _ in range(l)] for __ in range(n)]
    backTrack = [[0 for _ in range(l)] for __ in range(n)]
    for k in range(l):
        s[0][k] = log(1 / l) + emissionLog[k][x[0]]
    for i in range(1, n):
        for k in range(l):
            currS = [s[i - 1][kpre] + transionLog[kpre][k] + emissionLog[k][x[i]] for kpre in range(l)]
            ind = np.argmax(currS)
            backTrack[i][k] = ind
            s[i][k] = currS[ind]

    currState = np.argmax(s[n - 1])
    stateList = [currState]
    for i in range(n - 1, 0, -1):
        currState = backTrack[i][currState]
        stateList.insert(0, currState)
    path = ''.join([stateDict[state] for state in stateList])
    return path

x, transionLog, emissionLog, stateDict = readFromFile()
path = viterbi(x, transionLog, emissionLog, stateDict)
print(path)