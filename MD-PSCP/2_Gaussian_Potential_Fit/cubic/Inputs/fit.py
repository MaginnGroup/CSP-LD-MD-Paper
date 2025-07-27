
#----------------------------------------------------------------------------------

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pymoo.optimize import minimize
from pymoo.factory import get_termination

from pymoo.core.problem import ElementwiseProblem
from pymoo.algorithms.soo.nonconvex.pso import PSO

#----------------------------------------------------------------------------------

df = pd.read_csv('dados_dr0_P.txt', sep=' \t ', skiprows=0, names=['dr0', 'P'])

dr0 = df['dr0'].values
A_observado = df['P'].values

n = dr0.shape[0]

T = 350.0 #Kelvin
kb = 1.9872e-3 #kcal/mol/K
beta = 1/(kb*T)

#----------------------------------------------------------------------------------

def erro(x):
    a = x[0]
    f = 0
    for i in range(n):      
        A_analitica = (beta*a/np.pi)**(3/2)*4*np.pi*(dr0[i])**2*np.exp(-beta*a*(dr0[i])**2)
        f += (A_analitica - A_observado[i])**2   
    return f
    
#----------------------------------------------------------------------------------

class MyProblem(ElementwiseProblem):
            
    def __init__(self):
        super().__init__(n_var=1,
                          n_obj=1,
                          n_constr=0,
                          xl=np.array([0]),
                          xu=np.array([20]))

    def _evaluate(self, x, out, *args, **kwargs):
        f1 = erro(x)
        out["F"] = [f1]
    
problem = MyProblem()

algorithm = PSO(pop_size=100)


res = minimize(problem,
                algorithm,
                ("n_gen", 200),
                seed=1,
                verbose=False)
    
x_PSO = res.X
f_PSO = res.F

print('a = %e' % (x_PSO[0]))

#----------------------------------------------------------------------------------

a = x_PSO[0]

A_analitica = np.zeros(n)
for i in range(n):   
    A_analitica[i] = (beta*a/np.pi)**(3/2)*4*np.pi*(dr0[i])**2*np.exp(-beta*a*(dr0[i])**2) 

fig = plt.figure(figsize=(5,5))
plt.plot(dr0, A_analitica, 'k-', label='')
plt.plot(dr0, A_observado, 'bo', label='')
plt.show()

#----------------------------------------------------------------------------------



