import re
from sympy import symbols
from scipy.optimize import linprog

class  Equation:

    def  __init__(self, task, ogran):
        self.task = self.parse_str_task(task)
        self.ogran = self.parse_str_ogra(ogran)
        self.symvols = symbols(f'x:{len(self.task)}')

    def parse_str_task(self,task):
        splitx1 = re.split('x\d',task.get())
        # TODO убрать
        self.signs = re.split('\d',''.join(splitx1))
        splitx = list(filter(None,re.split('[+,-]',''.join(splitx1))))
        return splitx

    def parse_str_ogra(self,ogran):
        ogrns = []
        for ogra in ogran:
            splitx1 = re.split('x\d',ogra.get())
            signs = re.split('\d',''.join(splitx1))
            splitx = re.split('[+,-,=,<,>]',''.join(splitx1))
            ogrns.append((splitx,signs))
        return ogrns