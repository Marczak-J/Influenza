#!/usr/bin/env python
# coding: utf-8


DNAdict = {
    'G': { 'G': 1, 'C':-3, 'A':-3, 'T':-3, 'N':0 },
    'C': { 'G':-3, 'C': 1, 'A':-3, 'T':-3, 'N':0 },
    'A': { 'G':-3, 'C':-3, 'A': 1, 'T':-3, 'N':0 },
    'T': { 'G':-3, 'C':-3, 'A':-3, 'T': 1, 'N':0 },
    'N': { 'G': 0, 'C': 0, 'A': 0, 'T': 0, 'N':0 }
}


def SequenceAlign(seqA, seqB, similarityMatrix=DNAdict, insert=8, extend=4):
    
    import numpy as np
    
    numI = len(seqA) + 1
    numJ = len(seqB) + 1
    
    SMatrix = np.zeros((numI, numJ))
    RMatrix = np.zeros((numI, numJ))
    
    for i in range(1, numI):
        RMatrix[i, 0] = 1
        
    for j in range(1, numJ):
        RMatrix[0, j] = 2
    
    for i in range(1, numI):
        for j in range(1, numJ):
            
            penalty1 = insert
            penalty2 = insert
            
            if RMatrix[i-1, j] == 1:
                penalty1 = extend
                
            elif RMatrix[i, j-1] == 2:
                penalty2 = extend
                
            similarity = similarityMatrix[seqA[i-1]][seqB[j-1]]
            
            paths = [SMatrix[i-1, j-1] + similarity,
                     SMatrix[i-1, j] - penalty1,
                     SMatrix[i, j-1] - penalty2]
        
            best = max(paths)         #maximum value of path list
            route = paths.index(best) #index where maximum value
        
            SMatrix[i, j] = best  
            RMatrix[i, j] = route
                    
        alignA = []
        alignB = []
        
        i = numI-1
        j = numJ-1
            
        score = SMatrix[i, j]
        
        while i > 0 or j > 0:
            route = RMatrix[i, j]
            
            if route == 0: 
                alignA.append( seqA[i-1] )
                alignB.append( seqB[j-1] )
                i -= 1
                j -= 1
                
            elif route == 1:
                alignA.append( seqA[i-1] )
                alignB.append( '-' )
                i -= 1
                
            elif route == 2: 
                alignA.append( '-' )
                alignB.append( seqB[j-1] )
                j -= 1
                
    alignA.reverse()
    alignB.reverse()
    
    alignA = ''.join(alignA)
    alignB = ''.join(alignB)
    
    return score, alignA, alignB


def Profile(alignment):
    
    n = len(alignment[0])
    nS = len(alignment)
    profile = []
    
    for i in range(n):
        
        how_many = {}
        
        for seq in alignment:
            residue = seq[i]
            
            if residue == '-':
                continue
                
            how_many[residue] = how_many.get(residue, 0) + 1 
            
        for residue in how_many:
            how_many[residue] /= nS
            
        profile.append(how_many)
        
    return profile


def ProfileAlign(profileA, profileB, simiarityMatrix = DNAdict, insert=8, extend=4):
    
    import numpy as np
    
    numI = len(profileA) + 1
    numJ = len(profileB) + 1
    
    SMatrix = np.zeros((numI, numJ))
    RMatrix = np.zeros((numI, numJ))
    
    for i in range(1, numI):
        RMatrix[i,0] = 1
        
    for j in range(1, numJ):
        RMatrix[0,j] = 2

    for i in range(1, numI):
        for j in range(1, numJ):
            
            penalty1 = insert
            penalty2 = insert
            
            if RMatrix[i-1, j] == 1:
                penalty1 = extend
                
            elif RMatrix[i, j-1] == 2:
                penalty2 = extend
                
            frac_A = profileA[i-1]
            frac_B = profileB[j-1]
            
            similarity = 0
            totalWeight = 0
            
            for residue_A in frac_A:
                for residue_B in frac_B:
                    
                    weight = frac_A[residue_A] * frac_B[residue_B]
                    totalWeight += weight
                    similarity += weight * simiarityMatrix[residue_A][residue_B]
                    
            penalty1 *= totalWeight
            penalty2 *= totalWeight
            
            paths = [SMatrix[i-1, j-1] + similarity,
                     SMatrix[i-1, j] - penalty1,
                     SMatrix[i, j-1] - penalty2]
            
            best = max(paths) #maximum value of paths list
            route = paths.index(best) #index where maximum value
            
            SMatrix[i, j] = best
            RMatrix[i, j] = route
            
        pA = []
        pB = []
        
        i = numI-1
        j = numJ-1
        
        score = SMatrix[i, j]
        
        while i > 0 or j > 0:
            route = RMatrix[i, j]
            if route == 0: 
                pA.append(profileA[i-1])
                pB.append(profileB[j-1])
                i -= 1
                j -= 1
            elif route == 1: 
                pA.append(profileA[i-1])
                pB.append(None)
                i -= 1
            elif route == 2: 
                pA.append(None)
                pB.append(profileB[j-1])
                j -= 1
                
    pA.reverse()
    pB.reverse()
    
    return score, pA, pB


def ProfileMultipleAlignment(seqs, similarityMatrix = DNAdict):
    """
    This function returns Multiple Sequence Alignment (MSA)
    for a given list of sequences using profiles.
    """
    n = len(seqs)
    
    score, alignA, alignB = SequenceAlign(seqs[0], seqs[1], similarityMatrix) #alignment for two first sequences
    
    MSA = [alignA, alignB]
    
    for i in range(2,n):
        
        profA = Profile(MSA)
        toAdd = [seqs[i],] #the next sequence (3d, 4th and so on)
        profB = Profile(toAdd)
        score, alignA, alignB = ProfileAlign(profA, profB, similarityMatrix) #alignment between two profiles
        
        gaps = []
        
        for j, frac in enumerate(alignA):
            if frac is None:
                gaps.append(j)
                
        for j, seq in enumerate(MSA):
            for gap in gaps:
                seq = seq[:gap] + '-' + seq[gap:]
                MSA[j] = seq
                
        gaps = []
        
        for j, frac in enumerate(alignB):
            if frac is None:
                gaps.append(j)
                
        for j, seq in enumerate(toAdd):
            for gap in gaps:
                seq = seq[:gap] + '-' + seq[gap:]
            toAdd[j] = seq
    
        MSA.extend(toAdd)
                
    return MSA


def Score(align1, align2, similarityMatrix = DNAdict, gap = 3):
    n = len(align1)
    S = 0
    for i in range(n):
        if align1[i] == '-' or align2[i] == '-':
            S -= gap
        else:
            S += similarityMatrix[align1[i]][align2[i]]
    return S
            

def SimMatrix(alignment, similarityMatrix = DNAdict, gap = 3):
    import numpy as np
    
    n = len(alignment)
    A = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            M = max(Score(alignment[i], alignment[i], similarityMatrix, gap), Score(alignment[j], alignment[j], similarityMatrix, gap))
            A[i,j] = M - Score(alignment[i], alignment[j], similarityMatrix, gap) #the smaler value the more similar sequences are
    return A
