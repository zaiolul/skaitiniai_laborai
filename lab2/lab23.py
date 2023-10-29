import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from PyFunkcijos import *

def draw(base, new):
    limit = 10
    plt.plot(base[:,0],base[:,1], 'or', label='Esamos parduotuvės', linestyle = 'None')
    plt.plot(new[:,0],new[:,1], 'ob',   label='Naujos parduotuvės', linestyle = 'None')
   
    plt.xlim([-limit-2, limit+2])
    plt.ylim([-limit-2, limit+2])
    plt.plot([-limit, -limit, limit, limit, -limit], [-limit, limit, limit, -limit, -limit],'--k') 
    plt.grid()
    plt.legend()
    plt.show()

def distance(pos1, pos2):
    x1 = pos1[0]
    x2 = pos2[0]
    y1 = pos1[1]
    y2 = pos2[1]

    return np.exp(-0.3 *((x1 - x2)**2 + (y1 - y2)**2))

def loc_price(pos):
    x1 = pos[0]
    y1 = pos[1]

    return ((x1**4 + y1**4) / 1000) + ((np.sin(x1) + np.cos(y1)) / 5) + 0.4

def target(base, new):
    n_base = len(base)
    n_new = len(new)
    sum = 0
    for i in range(n_new):
        sum += loc_price(new[i])
        for j in  range(n_base):
            sum += distance(new[i], base[j])
            # rez = rez + (np.linalg.norm([X[i] - X[j], Y[i] - Y[j]]) - dist) ** 2 # sitoj vietoj tikslo funkcija

    for i in range(n_new):
        for j in range(i + 1, n_new):
            sum += distance(new[i], new[j])
    return sum

def numerical_gradient(base, new, h):
    f0 = target(base, new)
    G = np.zeros(shape=(len(new), 2))
    for i in range(len(new)):
        for j in range (2): # x and y coordinates
            # temp = np.zeros(shape=(len(new), 2))
            temp = np.copy(new)
            temp[i,j] += h
            f1 = target(base, temp)
            G[i, j] = (f1 - f0) / h;
    
    norm = np.linalg.norm(G)
    G = G / norm

    return G 


np.random.seed(11)

n_base = 5

n_new = 3

base_shops = np.random.uniform(-10, 10, (n_base, 2)) 
new_shops = np.random.uniform(-10, 10, (n_new, 2)) 

draw(base_shops, new_shops)

itmax = 1000
step = 0.2

h = 0.001

fx = target(base_shops, new_shops)

f_vals = []
its = []
for i in range(itmax):
    G = numerical_gradient(base_shops, new_shops, h)
    new_shops -=  step * G
    fx1 = target(base_shops, new_shops)
    f_vals.append(fx1)
    its.append(i)
    if fx1 > fx:
        new_shops +=  step * G
        step = step * 0.9
        print('fx1 = ', fx1, 'step = ', step, 'i =', i)
       
    else:
        fx = fx1

    if step < 1e-8:
        print("optimizavimas baigtas, fx =", fx, ",iteraciju skaicius:", i )
        break

draw(base_shops, new_shops)
print(new_shops)

plt.plot(its, f_vals)
plt.xlabel("iteracijos")
plt.ylabel("tikslo funkcijos reikšmė")
plt.show()