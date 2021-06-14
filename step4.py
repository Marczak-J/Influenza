DNA_2 = {'G': { 'G': 1, 'C':-3, 'A':-3, 'T':-3, 'N':0 },
'C': { 'G':-3, 'C': 1, 'A':-3, 'T':-3, 'N':0 },
'A': { 'G':-3, 'C':-3, 'A': 1, 'T':-3, 'N':0 },
'T': { 'G':-3, 'C':-3, 'A':-3, 'T': 1, 'N':0 },
'N': { 'G': 0, 'C': 0, 'A': 0, 'T': 0, 'N':0 }}


def WordSeq(seq):
    d=[]
    l=len(seq)
    i=0
    k=1
    n=0
    while i<l:
        while seq[i:i+k] in d and i+k<l:
            k+=1
    
        if seq[i:i+k] not in d:     
            d.append(seq[i:i+k])  
        i+=k   
        k=1
    return d
    
    
def LZcomplexity(seq1, seq2):
    N1 = len(WordSeq(seq1))
    N2 = len(WordSeq(seq2))
    N = len(WordSeq(seq1+seq2))
    C = (N - min(N1,N2))/max(N1,N2)
    return C
    
def LZMatrix(alignment, similarityMatrix = DNA_2, gap = 3):
    import numpy as np
    
    n = len(alignment)
    A = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            M = max(LZcomplexity(alignment[i], alignment[i]), LZcomplexity(alignment[j], alignment[j]))
            A[i,j] = M - LZcomplexity(alignment[i], alignment[j]) 
    return A