#! /usr/bin/env python

# Source: https://rosettacode.org/wiki/Longest_increasing_subsequence#Python

def longest_increasing_subsequence(X):
    """Returns the Longest Increasing Subsequence in the Given List/Array"""
    N = len(X)
    P = [0] * N
    M = [0] * (N+1)
    L = 0
    for i in range(N):
       lo = 1
       hi = L
       while lo <= hi:
           mid = (lo+hi)//2
           if (X[M[mid]] < X[i]):
               lo = mid+1
           else:
               hi = mid-1
 
       newL = lo
       P[i] = M[newL-1]
       M[newL] = i
 
       if (newL > L):
           L = newL
 
    S = []
    k = M[L]
    for i in range(L-1, -1, -1):
        S.append([X[k], k])
        k = P[k]
    return S[::-1]
