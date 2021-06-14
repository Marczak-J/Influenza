def MerDict(seq, k=3):
    d = {}
    N = len(seq)
    for i in range(N-k+1):
        mer = seq[i:(i+k)]
        if mer not in d.keys():
            d[mer] = 1
        else:
            d[mer] += 1         
    return d
    
def CompareSequences(seq1, seq2, k = 3):
    m1 = MerDict(seq1, k)
    m2 = MerDict(seq2, k)
    mers = []
    for elem in m1.keys():
        mers.append(elem)
    for elem in m2.keys():
        if elem not in mers:
            mers.append(elem)    
    seq1M = []
    seq2M = []
    for elem in mers:
        if elem in m1.keys():
            seq1M.append(m1[elem])
        else:
            seq1M.append(0)
            
        if elem in m2.keys():
            seq2M.append(m2[elem])
        else:
            seq2M.append(0)
    d = 0
    for i in range(len(mers)):
        d += (seq1M[i]-seq2M[i])**2
    d = d**0.5
        
        
    return d
    
    
def MerMatrix(alignment, similarityMatrix = DNA_2, gap = 3):
    import numpy as np
    
    n = len(alignment)
    A = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            M = max(CompareSequences(alignment[i], alignment[i], gap), CompareSequences(alignment[j], alignment[j], gap))
            A[i,j] = M - CompareSequences(alignment[i], alignment[j], gap) #the smaler value the more similar sequences are
    return A