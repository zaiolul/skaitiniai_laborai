import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
import math
import sympy as sp

#grubus ivertis
def fx(x):
    return -0.95 * x**4 + 10.19 * x**3 - 37.83 * x**2 + 55.58 * x - 24.49 

def gx(x):
    return sp.sin(x)**2 * sp.log(x) - (x / 4)

def get_roots(f, pairs, n, eps, method = "chords", method_sym = False):
    ans = []
    if method == 'chords':
        for pair in pairs:
            ans.append(chords(f, pair[0], pair[1], n, eps, method_sym))
    elif method == 'secant':
        for pair in pairs:
            ans.append(secant(f, pair[0], pair[1], n, eps, method_sym))
    return ans

def scan(f, x_start, x_end, sym = False, step = 0.001):
    bracket_pairs = []
    x_current = x_start

    if not sym:
        while x_current <= x_end:
            if np.sign(f(x_current)) !=np.sign(f(x_current + step)):
                bracket_pairs.append([x_current, x_current + step])
            x_current += step
    else:
        while x_current <= x_end:
            x = sp.symbols('x')
            if np.sign(f.evalf(subs={x: x_current})) !=np.sign(f.evalf(subs={x: x_current + step})):
                bracket_pairs.append([x_current, x_current + step])
            x_current += step
    return bracket_pairs

#STYGU METODAS
def chords(f, x_start, x_end, n, eps, sym = False):
    count = 0
    f_start, f_end = 0,0
    x = sp.symbols('x')
    for i in range(n):
        if not sym:
            f_start = f(x_start)
            f_end = f(x_end)
        else:
            f_start = f.evalf(subs={x: x_start})
            f_end = f.evalf(subs={x: x_end})

        k = np.abs(f_start / f_end)
        x_mid = (x_start + k * x_end) / (1 + k)

        if np.sign(x_mid) != np.sign(x_end):
            x_end = x_mid
        else:
            x_start = x_mid
        count = count + 1
        if abs(x_end - x_start) < eps:
            break
        
    print(count, "chords iteraciju")
    return x_mid

#KIRSTINIU METODAS
def secant(f, x_start, x_end, n, eps, sym = False):
    count = 0
    f_start, f_end = 0, 0
    x = sp.symbols("x")
    for i in range(n):
        if not sym:
            f_start = f(x_start)
            f_end = f(x_end)
        else:
            f_start = f.evalf(subs={x: x_start})
            f_end = f.evalf(subs={x: x_end})
    
        xi = x_end - f_end * (x_end - x_start) / (f_end - f_start)
        x_start = x_end
        x_end = xi
        count = count + 1
        if abs(x_end - x_start) < eps:
            
            break
    print(count, "secant iteraciju")    
    return xi

def print_func_value(f, roots, syms = False):
    if not syms:
        for root in roots:
            print(f(root))
    else:
        x = sp.symbols('x')
        for root in roots:
            print(f(x).evalf(subs={x:root}))

def main():
    #vyriausio laipsnio narys neigiamas, pertvarkoma
    # 0.95 * x**4 - 10.19 * x**3 + 37.83 * x**2 - 55.58 * x + 24.49 
    R = (55.58 / 0.95) + 1 

    #tikslus ivertis teigiamoms saknims
    k = 4 - 3  # vyriausias neigiamas 3 eiles (10.19), vyriausias ketvirtas
    B = 55.58
    R_pos = (55.58 / 0.95) + k 

    #tikslus ivertis neigiamoms saknims

    # f(-x) = -0.95 * x**4 - 10.19 * x**3 - 37.83 * x**2 - 55.58 * x - 24.49
    # 0.95 * x**4 + 10.19 * x**3 + 37.83 * x**2 + 55.58 * x + 24.49
    #nera neigiamu koeficientu, todel R = 0
    R_neg = 0

    x_start = -1 * min(R, R_neg)
    x_end = min(R, R_pos)

    print("Galutinis ivertis: ", x_start, " <= x <= ", x_end)

    # PIESIMAS

    #f(x)
    fx_range = [-60, 60]
    fx_lin = np.linspace(fx_range[0], fx_range[1], 500)

    plt.ylim(-5,5)
    plt.plot(fx_lin, fx(fx_lin), '-b')
    plt.plot(fx_range, [0, 0], 'k--')
    plt.plot(-R, 0, 'sc')
    plt.plot(R, 0, 'sc')
    plt.plot(R_neg, 0, 'or')
    plt.plot(R_pos, 0, 'or')
    plt.title("f(x)")
    plt.legend(['f(x)', 'Ox', 'Grubus įv', 'Tikslesnis įv.'])
    plt.show()

    # g(x)
    gx_range = [1, 10]
    gx_lin = np.linspace(gx_range[0], gx_range[1], 500)
    x = sp.symbols("x")
    gx_points = []
    for i in range(500):
        gx_points.append(gx(x).subs(x, gx_lin[i]))

    plt.plot(gx_lin, gx_points, '-g')
    plt.plot(gx_range, [0, 0], 'k--')
    plt.legend(['g(x)', 'Ox'])
    plt.title("g(x)")
    plt.show()


    # INTERVALU SKENAVIMAS

    fx_pairs = scan(fx, x_start, x_end, step = 0.01)
    gx_pairs = scan(gx, 1, 10, step = 0.01)
    print("fx_pairs", fx_pairs)
    print("gx_pairs", gx_pairs)
    eps = 1e-6
    n_iter = 500

    fx_chords, fx_secant = get_roots(fx, fx_pairs, n_iter, eps, method = 'chords'), get_roots(fx, fx_pairs, n_iter, eps, method = 'secant')
    gx_chords, gx_secant = get_roots(gx, gx_pairs, n_iter, eps, method = 'chords'), get_roots(gx, gx_pairs, n_iter, eps, method = 'secant')
   
    print("fx chords", fx_chords)
    print_func_value(fx, fx_chords)
    print("fx secant",fx_secant)
    print_func_value(fx, fx_secant)
    print("fx numpy roots", np.roots([-0.95, 10.19, -37.83, 55.28, -24.49]))
    
    # print("gx chords", gx_chords)
    # print_func_value(gx, gx_chords, syms = True)
    # print("gx secant", gx_secant)
    # print_func_value(gx, gx_secant, syms = True)

if __name__ == "__main__":
    main()