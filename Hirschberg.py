import argparse
import time

gapPenalty = {'A': 1, 'T': 2, 'G': 1, 'C': 3}
simMatrix = {'A': {'A': 0, 'T': 1, 'G': 5, 'C': 1},
              'T': {'A': 1, 'T': 0, 'G': 9, 'C': 1},
              'G': {'A': 5, 'T': 9, 'G': 0, 'C': 1},
              'C': {'A': 1, 'T': 1, 'G': 1, 'C': 0}}

def nw(A, B, simMatrix, gapPenalty):
    # The Needleman-Wunsch algorithm
    # Stage 1: Create a zero matrix and fills it via algorithm
    n, m = len(A), len(B)
    mat = []
    for i in range(n+1):
        mat.append([0]*(m+1))
    for j in range(1, m+1):
        mat[0][j] = gapPenalty[B[j-1]]
        mat[0][j] += mat[0][j-1]
    for i in range(1, n+1):
        mat[i][0] = gapPenalty[A[i-1]]
        mat[i][0] += mat[i-1][0]
    for i in range(1, n+1):
        for j in range(1, m+1):
            mat[i][j] = min(mat[i-1][j-1] + simMatrix[A[i-1]][B[j-1]], mat[i][j-1] + gapPenalty[B[j-1]], mat[i-1][j] + gapPenalty[A[i-1]])

    # Stage 2: Computes the final alignment, by backtracking through matrix
    alignmentA = ""
    alignmentB = ""
    i, j = n, m
    while i and j:
        score, scoreDiag, scoreUp, scoreLeft = mat[i][j], mat[i-1][j-1], mat[i-1][j], mat[i][j-1]
        if score == scoreDiag + simMatrix[A[i-1]][B[j-1]]:
            alignmentA = A[i-1] + alignmentA
            alignmentB = B[j-1] + alignmentB
            i -= 1
            j -= 1
        elif score == scoreUp + gapPenalty[A[i-1]]:
            alignmentA = A[i-1] + alignmentA
            alignmentB = '-' + alignmentB
            i -= 1
        elif score == scoreLeft + gapPenalty[B[j-1]]:
            alignmentA = '-' + alignmentA
            alignmentB = B[j-1] + alignmentB
            j -= 1
    while i:
        alignmentA = A[i-1] + alignmentA
        alignmentB = '-' + alignmentB
        i -= 1
    while j:
        alignmentA = '-' + alignmentA
        alignmentB = B[j-1] + alignmentB
        j -= 1
    return [alignmentA, alignmentB, mat[n][m]]

def forwards(x, y, simMatrix, gapPenalty):
    # This is the forwards subroutine.
    n, m = len(x), len(y)
    mat = []
    for i in range(n+1):
        mat.append([0]*(m+1))
    for j in range(1, m+1):
        mat[0][j] = gapPenalty[y[j-1]]
        mat[0][j] += mat[0][j-1]
    for i in range(1, n+1):
        mat[i][0] = mat[i-1][0] + gapPenalty[x[i-1]]
        for j in range(1, m+1):
            mat[i][j] = min(mat[i-1][j-1] + simMatrix[x[i-1]][y[j-1]],
                            mat[i-1][j] + gapPenalty[x[i-1]],
                            mat[i][j-1] + gapPenalty[y[j-1]])
        # Now clear row from memory.
        mat[i-1] = []
    return mat[n]

def backwards(x, y, simMatrix, gapPenalty):
    # This is the backwards subroutine.
    n, m = len(x), len(y)
    mat = []
    for i in range(n+1):
        mat.append([0]*(m+1))
    for j in range(m+1):
        mat[0][j] = gapPenalty[y[j-1]]
        mat[0][j] += mat[0][j-1]
    for i in range(1, n+1):
        mat[i][0] = mat[i-1][0] + gapPenalty[x[i-1]]
        for j in range(1, m+1):
            mat[i][j] = min(mat[i-1][j-1] + simMatrix[x[n-i]][y[m-j]],
                            mat[i-1][j] + gapPenalty[x[i-1]],
                            mat[i][j-1] + gapPenalty[y[j-1]])
        # Now clear row from memory.
        mat[i-1] = []
    return mat[n]

def hirschberg(x, y, simMatrix, gapPenalty):
    # This is the main Hirschberg routine.
    n, m = len(x), len(y)
    if n<2 or m<2:
        # In this case we just use the N-W algorithm.
        return nw(x, y, simMatrix, gapPenalty)
    else:
        # Make partitions, call subroutines.
        F, B = forwards(x[:n//2], y, simMatrix, gapPenalty), backwards(x[n//2:], y, simMatrix, gapPenalty)
        partition = [F[j] + B[m-j] for j in range(m+1)]
        cut = partition.index(min(partition))
        # Clear all memory now, so that we don't store data during recursive calls.
        F, B, partition = [], [], []
        # Now make recursive calls.
        callLeft = hirschberg(x[:n//2], y[:cut], simMatrix, gapPenalty)
        callRight = hirschberg(x[n//2:], y[cut:], simMatrix, gapPenalty)
        # Now return result in format: [1st alignment, 2nd alignment, similarity]
        return [callLeft[r] + callRight[r] for r in range(3)]
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("file1", help="file 1 path")
    parser.add_argument("file2", help="file 2 path")
    args = parser.parse_args()
    file1_path = args.file1
    file2_path = args.file2
    
    # toy case
    print("*********** A Toy Example *************")
    str1 = 'ATGCATGC'
    str2 = 'TGCAGC'
    
    st = time.time()
    AA, BB, ss = hirschberg(str1, str2, simMatrix, gapPenalty)
    et = time.time()
    print("time cost: {}s".format(et - st))
    print("score: ", ss)
    print(AA)
    print(len(str1)*'|')
    print(BB)
    print()
    
    # long string case
    print("************ Begin Real Inputs *************")
    with open(file1_path, 'r') as fin:
        a_all = ''
        while(fin.readline()):
            a_all += fin.readline().strip()
    with open(file2_path, 'r') as fin:
        b_all = ''
        while(fin.readline()):
            b_all += fin.readline().strip()

    st = time.time()
    AA, BB, ss = hirschberg(a_all, b_all, simMatrix, gapPenalty)
    et = time.time()
    print("time cost: {}s".format(et - st))
    print("score: ", ss)
    print(A)
    print(len(A)*'|')
    print(B)

