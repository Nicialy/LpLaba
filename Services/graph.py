from sympy import *
from fractions import Fraction
import copy

class GaussAlgorithm:
    
    def __init__(self, A, rhs, f0):
        self.A = A.copy()
        self.rhs = rhs.copy()
        self.f0 = f0

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
    

    def check_matrix(self):   
        self.indexs = []
        A, rhs = self.A, self.rhs
        countOne = []
        for col in range(A.cols):  
            rowflag = 0
            i = 0
            for row in range(A.rows):
              if A[row,col] > 0 or A[row,col] < 0 :
                 i+=1
                 rowflag = row
              if  A[row,col] == 1:
                  countOne.append((rowflag,col))
            if i == 1:
                self.indexs.append((rowflag,col))
        
        if A.cols - 2 == len(self.indexs) or len(countOne) == A.cols - 2:
            return True
        
        return False
    
    def podstanovka(self):
        self.expressions = []
        self.constraint = []
        self.create_symbols()
        for expression in self.indexs:
            x = self.symbolsList[expression[1]]
            y = self.symbolsList[expression[1]]
            row = self.A.row(expression[0])
            for col,koeff in enumerate(row) :
                if col != expression[1] and koeff != 0:
                    x += -1*koeff * self.symbolsList[col]
                    y +=  1*koeff * self.symbolsList[col]
            x -=  self.symbolsList[expression[1]]
            y -=  self.symbolsList[expression[1]]
            self.constraint.append(y <= S(self.rhs.row(expression[0])[0]))
            x += S(self.rhs.row(expression[0])[0])
            
            self.symbolsList[expression[1]] = x
        self.create_new_f()
        self.expressions = copy.deepcopy(self.symbolsList)
        for index,expression in enumerate(self.symbolsList):
            self.symbolsList[index] = expression 
           
            
    def create_new_f(self):
        f = symbols("f")
        for index ,koeff in enumerate(self.f0):
            f +=  koeff * ( self.symbolsList[index])
        self.f0 = f -symbols("f")  

    def create_symbols(self):
        self.symbolsList = []
        for col in range(self.A.cols):
            x = symbols(f"x_{col + 1}")  
            self.symbolsList.append(x)

    def recursion(self, n=0):
        if n == len(self.symbolsList) - 1:
            return self.symbolsList[n]
        return self.symbolsList[n] & self.recursion(n+1)



# проверка
if __name__ == "__main__":
    dsa = GaussAlgorithm(Matrix([[1,2,5,-1],[1,-1,-1,2]]),Matrix([4,1]),Matrix([-2,-1,-3,-1]))
    dsa.doit()
    if dsa.check_matrix():
        dsa.podstanovka()
        #dsa.symbolsList.append()
        from cvxopt.modeling import variable, op   
        print(dsa.constraint)
        my_str = [str(expr).replace('*', '').replace('-','- ').replace(' x',' 1x').replace('-  ','- ') for expr in dsa.constraint]
        my_str = [strr if strr[0] != 'x' else '1' + strr[:] for strr in my_str ]
        print(my_str)
        if str(dsa.f0)[0] == 'x':
           f0 = '1' + str(dsa.f0)[:]
        else: f0 = str(dsa.f0)
        from simplex import Simplex
        from tkinter import *
        objective = ('min',f0)
        from tkinter import scrolledtext
        text = scrolledtext.ScrolledText( state='disable')
        Lp = Simplex(count_vars=2, constraints=my_str, objective_fun=objective,text = text, char = END)
        dsa2 = Lp.solution
        dsa2.pop("val")
        print(dsa2)
        print(dsa.f0.subs([(dsa.expressions[-2] , int(dsa2["x_1"])),(dsa.expressions[-1] , int(dsa2["x_2"]))]))

