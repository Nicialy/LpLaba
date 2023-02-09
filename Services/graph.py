from sympy import *
from fractions import Fraction

class GaussAlgorithm:
    
    def __init__(self, A, rhs):
        self.A = A.copy()
        self.rhs = rhs.copy()

    def forward_step(self, row_idx):
        A, rhs = self.A, self.rhs
        i = row_idx

        if A[i, i] == 0:        
            for j in range(i, A.rows):
                if A[j, i] != 0:
                    A[i, :], A[j, :] = A[j, :], A[i, :]
                    rhs[i, :], rhs[j, :] = rhs[j, :], rhs[i, :]
            
                    break


        rhs[i, :] /= A[i, i]
        A[i, :] /= A[i, i]



        for j in range(i+1, A.rows):
            if A[j, i] != 0:
                rhs[j, :] -= A[j, i] * rhs[i, :]
                A[j, :] -= A[j, i] * A[i, :]

        return A, rhs

    def backward_step(self, row_idx):
        A, rhs = self.A, self.rhs
        i = row_idx
        for j in range(i-1, -1, -1):
            rhs[j, :] -= A[j, i] * rhs[i, :]
            A[j, :] -= A[j, i] * A[i, :]


    def doit(self):
        A, rhs = self.A, self.rhs
        for row in range(A.rows):
            self.forward_step(row)
        for row in range(A.rows - 1, -1, -1):
            self.backward_step(row)
        return rhs

    