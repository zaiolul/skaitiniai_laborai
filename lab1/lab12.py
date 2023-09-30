#18 −111 sin(𝑥) + 1 + 𝑥2 ; −2 ≤ 𝑥 ≤ 8

import sympy as sp
import matplotlib.pyplot as plt
import numpy as np
import math
from lab11 import scan, get_roots

def hx(x):
    return -111 * sp.sin(x) + 1 + x**2

def taylor(f, x, x0, n):
    exp = f.subs(x, x0)
    df = f
    for i in range(1, n):
        df = df.diff(x);
        exp = exp + df.evalf(subs={x: x0}) * (x - x0) ** i / math.factorial(i)
    return exp


def plot_func(f, interval, vals = 500, line_type = '-k', limits = True, show = False):
    f_lin = np.linspace(interval[0], interval[1], vals)
    f_points = []
    for i in range(vals):
        f_points.append(f.evalf(subs = {'x': f_lin[i]}))
    if limits:
        plt.ylim([-200, 200])
    plt.plot(f_lin, f_points, line_type)
    if show:
        plt.xlabel("x")
        plt.ylabel("y")
        plt.show()

def find_series(x, x0, f_roots, interval):
    n = 1
    print("Finding Taylor series count...")
    all_roots = []
    while True:
        exp = taylor(hx(x), x, x0, n)
        exp_pairs = scan(exp, interval[0], interval[1], sym = True, step = 0.1)
        exp_roots = get_roots(exp, exp_pairs, 100, 1e-6, method = "secant", method_sym = True)
        all_roots.append(exp_roots)
        n = n + 1
        if len(exp_roots) != len(f_roots):
            continue
        
        found = True
        for i in range(len(f_roots)):
            if abs(f_roots[i] - exp_roots[i]) > 1e-4:
                found = False
        
        if found:
            return exp, all_roots, n
        

def main():
    x, exp = sp.symbols("x exp")
    x0 = 3 # -2 <= x <= 8

    hx_range = [-2, 8]

    # 2.1, atidaro skirtingus langus grafiku

    for i in range(3, 6):
        exp = taylor(hx(x), x, x0, i)
        plt.plot(hx_range, [0, 0], '--k')
        plot_func(hx(x), hx_range, line_type = '-b')
        plot_func(exp, hx_range, line_type = '-r')
        plt.legend(['Ox', 'hx', 'TE'])
        plt.xlabel("x")
        plt.ylabel("y")
        plt.show()
 
    # 2.2
    hx_pairs = scan(hx, hx_range[0], hx_range[1], step = 0.1)
    hx_roots = get_roots(hx, hx_pairs, 200, 1e-6, method = 'secant')
    print('hx roots (secant):', hx_roots)
    exp, exp_roots, n = find_series(x, x0, hx_roots, hx_range)
    print('Teiloro eilutes numeris: {0}\nTeiloro eilutes saknys: {1}\nTeiloro eilute: {2}'.format(n, exp_roots, exp))

    new_range = [-15, 15]
    plt.plot(new_range, [0, 0], '--k')
    plot_func(hx(x), new_range, line_type = '-b')
    plot_func(exp, new_range, line_type = '-r', limits = True)
    plt.legend(['Ox', 'hx', 'TE'])
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()
    
    #2.3 spausdina

    #2.4 a

    root_n = []
    for i in range(n - 1):
        root_n.append(len(exp_roots[i]))
    
    plt.plot(np.linspace(1,n - 1, n-1), root_n, '-r')
    plt.xticks(range(0, n))
    plt.xlabel("TE eilė")
    plt.ylabel("Šaknų skaičius")
    plt.show()

    # 2.4 b

    diffs = []
    for i in range(len(hx_roots)):
        series = [0]
        for j in range(1, len(exp_roots)):
            closest = abs(exp_roots[j][0]) # default value to compare against
            for k in range(len(exp_roots[j])):
                abs_diff = abs(hx_roots[i] - exp_roots[j][k])
                if abs_diff < closest:
                    closest = abs_diff
            series.append(closest)
        diffs.append(series)
        
        plt.plot(np.linspace(1, n - 1, n - 1), diffs[i], '-r')
        plt.xticks(range(0, n))
        plt.xlabel("TE eilė")
        plt.ylabel("{0} h(x) šaknies ir artimiausios TE šaknies skirtumas".format(i))
        plt.show()

        
 
            

    

if __name__ == "__main__":
    main()