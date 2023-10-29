import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from PyFunkcijos import *
from scipy.optimize import fsolve


def lf(x):
    m = np.array(
        [
            ((x[0] - 3) ** 3 + x[1] - 8),
            (((x[0] ** 2 + x[1] ** 2) / 2) - 6 * (np.cos(x[0]) + np.cos(x[1])) - 10),
        ]
    )
    m.shape = (2, 1)
    m = np.matrix(m)
    return m


def plot_surfaces(X, Y, Z):
    fig1 = plt.figure(1, figsize=plt.figaspect(0.5))
    ax1 = fig1.add_subplot(1, 2, 1, projection="3d")
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_ylabel("z")
    ax2 = fig1.add_subplot(1, 2, 2)
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    ax2.set_ylabel("z")

    ax1.plot_surface(
        X, Y, Z[:, :, 0], color="blue", alpha=0.4, linewidth=0.1, antialiased=True
    )
    ax1.contour(X, Y, Z[:, :, 0], [0], colors="b")
    ax1.plot_surface(
        X, Y, Z[:, :, 1], color="purple", alpha=0.4, linewidth=0.1, antialiased=True
    )
    ax1.contour(X, Y, Z[:, :, 1], [0], colors="g")
    ax2.contour(X, Y, Z[:, :, 0], [0], colors="b")
    ax2.contour(X, Y, Z[:, :, 1], [0], colors="g")
    plt.show()


def broiden(X, Y, Z, rn):
    n = 2
    x = np.matrix(np.zeros(shape=(n, 1)))
    maxiter = 50
    eps = 1e-6
    plt.ylim(-12, 12)
    plt.xlim(-12, 12)
    plt.xlabel("x")
    plt.ylabel("y")

    XX = np.linspace(rn[0], rn[1], 2)
    YY = XX
    XX, YY = np.meshgrid(XX, YY)
    ZZ = XX * 0

    plt.draw()

    sols = np.empty(shape=(rn[1] * 2, 2))
    sol_count = 0
    colors = ["r", "c", "k"]

    for i in range(rn[0], rn[1] + 1):
        for j in range(rn[0], rn[1] + 1):
            x[0] = i
            x[1] = j
            found = False

            dx = 0.1
            A = np.matrix(np.zeros(shape=(n, n)))
            x1 = np.zeros(shape=(n, 1))
            for k in range(0, n):
                x1 = np.matrix(x)
                x1[k] += dx
                A[:, k] = (lf(x1) - lf(x)) / dx

            ff = lf(x)

            for k in range(1, maxiter):
                deltax = -np.linalg.solve(A, ff)
                x1 = np.matrix(x + deltax)
                ff1 = lf(x1)
                A += (
                    (ff1 - ff - A * deltax)
                    * deltax.transpose()
                    / (deltax.transpose() * deltax)
                )
                tiksl = tikslumas(x, x1, ff, ff1, eps)
                ff = ff1
                x = x1
                if tiksl < 1e-12:
                    break

            if abs(lf(x)[0, 0]) > eps and abs(lf(x)[1, 0]) > eps:
                plt.plot(i, j, colors[-1] + ".", markersize=10)
                continue

            for k in range(sol_count):
                if abs(sols[k][0] - x[0, 0]) < eps and abs(sols[k][1] - x[1, 0]) < eps:
                    found = True
                    plt.plot(i, j, colors[k] + ".", markersize=10)
                    break

            if not found:
                print("Sprendinys: [%.5f %.5f] (%d %d)" % (x[0,0], x[1,0], i, j))
                sols[sol_count][0] = x[0]
                sols[sol_count][1] = x[1]

                plt.plot(i, j, colors[sol_count] + ".", markersize=10)
                plt.plot(
                    x[0, 0],
                    x[1, 0],
                    colors[sol_count] + "*",
                    markersize=10,
                    markeredgecolor=(0, 0, 0, 1),
                )
                sol_count = sol_count + 1

            plt.draw()
    plt.contour(X, Y, Z[:, :, 0], [0], colors="b")
    plt.contour(X, Y, Z[:, :, 1], [0], colors="g")

    # x[0] = 2
    # x[1] = 10
    # plt.plot(x[0], x[1] -0.4,"^")

    # dx=0.1
    # A = np.matrix(np.zeros(shape=(n,n)))
    # x1=np.zeros(shape=(n,1));
    # for k in range (0,n):
    #     x1 = np.matrix(x);
    #     x1[k]+=dx;
    #     A[:,k]=(lf(x1)-lf(x))/dx

    # ff=lf(x)

    # for k in range (1, maxiter):
    #     plt.plot(x[0],x[1], "k^-", markersize= 10 - k )
    #     deltax = -np.linalg.solve(A, ff);
    #     x1 = np.matrix(x + deltax);
    #     ff1 = lf(x1)
    #     A += (ff1 - ff - A * deltax) * deltax.transpose() / (deltax.transpose() * deltax);
    #     tiksl = tikslumas(x,x1,ff,ff1,eps);
    #     ff = ff1;
    #     x = x1;
    #     if tiksl < 1e-12:
    #         break;

    # print(lf(x))
    # print(x)

    # plt.plot(x[0], x[1]-0.4,"^")
    plt.show()


print("2a ir b")

rn = [-10, 10]
xx = np.linspace(rn[0], rn[1], 100)
yy = np.linspace(rn[0], rn[1], 100)
Z = np.zeros(shape=(len(xx), len(yy), 2))
X, Y = np.meshgrid(xx, yy)
for i in range(0, len(xx)):
    for j in range(0, len(yy)):
        Z[i, j, :] = lf([X[i][j], Y[i][j]]).transpose()

plot_surfaces(X,Y,Z)

print("\n2c")

broiden(X, Y, Z, rn)

print("\n2d")

def eq(x):
    res = lf(x)
    return [res[0,0], res[1,0]] 

def check():
    print("Tikrinimas su fsolve:")
    print(fsolve(eq, (-10, -9)))
    print(fsolve(eq, (-10, -7)))

check()

