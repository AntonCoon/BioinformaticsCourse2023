import numpy as np
import binarytree


def tree_build(a: str, b: str, dl: int, ins: int, mism: int, match: int):
    
    def hirshberg(tree: binarytree, a: str, b: str):
        
        if len(a) <= 1 or len(b) <= 1:
            return
        
        mid = len(a)//2
        idx = np.argmax(Needleman_Wunsch_Algorithm(a[:mid], b, dl, ins, mism, match) + Needleman_Wunsch_Algorithm(a[mid:][::-1], b[::-1], dl, ins, mism, match)[::-1])
        tree.left = binarytree.Node(str((a[:mid], b[:idx])))
        tree.right = binarytree.Node(str((a[mid:], b[idx:])))
        hirshberg(tree.left, a[:mid], b[:idx])
        hirshberg(tree.right, a[mid:], b[idx:])
        
    tree = binarytree.Node(str((a, b)))
    hirshberg(tree, a, b)
    return tree



def Needleman_Wunsch_Algorithm(a: str, b: str, dl: int, ins: int, mism: int, match: int):
    a = "-" + a
    b = "-" + b
    M = np.zeros((2, len(b)))
    
    for i in range(1, len(a)):
        M[1, 0] = dl*i

    for j in range(1, len(b)):
        M[1, j] = ins*j

    for i in range(len(a)):
        for j in range(len(b)):
            D = M[0, j] + dl
            I = M[1, j-1] + ins
            M[1, j] = max(D, I, M[0, j-1] + (mism if a[i] != b[j] else match))
        M[0, :] = M[1, :]
        
    return M[1, :]


if __name__ == "__main__":
    a, b = input(), input()
    delete_cost, insert_cost, mismatch_cost,  match_reward= map(int, input().split())

    print(tree_build(a, b, delete_cost, insert_cost, mismatch_cost, match_reward))
