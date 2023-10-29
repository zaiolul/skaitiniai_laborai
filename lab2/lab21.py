import numpy as np
import sympy as sp
import matplotlib.pyplot as plt


def print_matrix(mat):
    str = "[\n"
    for i in range(0, np.shape(mat)[0]):
        str = str + "["
        for j in range(np.shape(mat)[1]):
            str = str + "%20.5f " % mat[i, j]
        str = str + "]\n"
    str = str + "]"
    print(str)


def atspindzio(mat, free, eps=1e-12):
    n = (np.shape(mat))[0]
    nb = (np.shape(free))[1]

    for i in range(n - 1):
        z = mat[i:n, i]
        if not all(z):
            continue

        ref_z = np.zeros(np.shape(z))
        ref_z[0] = np.sign(z[0]) * np.linalg.norm(z)
        omega = (z - ref_z) / np.linalg.norm(z - ref_z)
        Q = np.identity(n - i) - 2 * omega * omega.transpose()
        mat[i:n] = Q.dot(mat[i:n])

    x = sp.zeros(n, nb)

    for i in range(n - 1, -1, -1):
        calc = mat[i, n : n + nb] - mat[i, i + 1 : n] * x[i + 1 : n, :]
        if abs(mat[i, i]) < eps and calc.evalf(subs={"p0:%d".format(n): 1})[0] < eps:
            p = sp.symbols("p0:{0}".format(n))
            x[i] = p[i]
        elif abs(mat[i, i]) < eps and calc.evalf(subs={"p0:%d".format(n): 1})[0] > eps:
            return -1
        else:
            x[i, :] = calc / mat[i, i]
    return x


def gauso_zeidelio(mat, free):
    alpha = np.array([1.5, 1.5, 1.5, 1.5])
    n = np.shape(mat)[0]
    mat_tld = np.diag(1.0 / np.diag(mat)).dot(mat) - np.diag(alpha)
    free_tld = np.diag(1.0 / np.diag(mat)).dot(free)

    n_iter = 1000
    eps = 1e-12

    res = np.zeros(shape=(n, 1))
    res1 = np.zeros(shape=(n, 1))

    x = []
    y = []
    prec = 0

    for i in range(n_iter):
        res1 = ((free_tld - mat_tld.dot(res)).transpose() / alpha).transpose()
        prec = np.linalg.norm(res1 - res) / (np.linalg.norm(res) + np.linalg.norm(res1))

        y.append(prec)
        x.append(i)
        if prec < eps:
            break5
        res[:] = res1[:]

    print(prec)
    print_matrix(res)
    plt.scatter(x, y, marker="*")
    plt.yscale("log")
    plt.show()
    return res

def check(mat, x):
    n = np.shape(mat)[0]
    p = sp.symbols("p0:{0}".format(n))
    if isinstance(x, sp.MutableDenseMatrix):
        res = mat[0:n, :] * x.evalf(subs={sym: 1 for sym in p})
        return np.array(res).astype(float).round(12)
    else:
        res = mat[0:n, :] * x
        return res.round(12)


def lu(mat, free, P):
    n = np.shape(mat)[0]
    nb = np.shape(free)[1]
    for i in range(0, n - 1):  # range pradeda 0 ir baigia n-2 (!)
        largest = abs(mat[i:n, i]).argmax()
        mat[[i, i + largest], :] = mat[[i + largest, i], :]  # sukeiciamos eilutes
        P[[i, i + largest]] = P[[i + largest, i]]  # sukeiciami eiluciu numeriai
        for j in range(i + 1, n):  # range pradeda i+1 ir baigia n-1
            r = mat[j, i] / mat[i, i]
            mat[j, i : n + nb] = mat[j, i : n + nb] - mat[i, i : n + nb] * r
            mat[j, i] = 0

    free = free[P, :]
    # 1-as atgalinis etapas, sprendziama Ly=b, y->b
    for i in range(1, n):
        free[i] = free[i] - mat[i, 0:i] * free[0:i]
    # 2-as atgalinis etapas , sprendziama Ux=b, x->b
    for i in range(n - 1, -1, -1):
        free[i] = (free[i] - mat[i, i + 1 : n] * free[i + 1 : n]) / mat[i, i]
    return free


# DATA
coef_1 = np.matrix([[3, 7, 1, 3], [-3, 8, 2, 1], [4, 4, -7, 1], [1, -6, 6, 9]]).astype(
    float
)
free_1 = (np.matrix([37, 0, 38, 11])).transpose().astype(float)
mat_1 = np.hstack((coef_1, free_1))

coef_14 = np.matrix(
    [[2, 4, 6, -2], [1, 3, 1, -3], [1, 1, 5, 1], [2, 3, -3, -2]]
).astype(float)
free_14 = (np.matrix([4, -7, 11, -4])).transpose().astype(float)
mat_14 = np.hstack((coef_14, free_14))

coef_20 = np.matrix(
    [[2, 4, 6, -2], [1, 3, 1, -3], [1, 1, 5, 1], [2, 3, -3, -2]]
).astype(float)
free_20 = (np.matrix([2, 1, 7, 2])).transpose().astype(float)
mat_20 = np.hstack((coef_20, free_20))


all_data = [
    (mat_1, coef_1, free_1),
    (mat_14, coef_14, free_14),
    (mat_20, coef_20, free_20),
]

print("1a")

print("Atspindžio metodas:\n")

for i in range(len(all_data)):
    print("\n\n\n%d matrica:" % (i + 1))
    tup = all_data[i]
    copy = np.copy(tup[0], order="K")
    res = atspindzio(tup[0], tup[2])
    if res == -1:
        print("Matrica neturi sprendinių")
        print("Išorinio įrankio rezultatas:")
        try:
            print(np.linalg.solve(tup[1], tup[2]))
        except:
            print("Singuliari matrica, neina spręsti su solve")
        continue
    print("Lygčių sistemos sprendiniai:")
    print(res)

    print("Išorinio įrankio rezultatas:")
    try:
        print(np.linalg.solve(tup[1], tup[2]))
    except:
        print("Singuliari matrica, neina spręsti su solve")
    print("Įstatyti sprendiniai:")
    print(check(tup[1], res))


print("\n Gauso-Zeidelo metodas:\n")
res = gauso_zeidelio(coef_1, free_1)
print(res)

try:
    print(np.linalg.solve(tup[1], tup[2]))
except:
    print("Singuliari matrica, neina spręsti su solve")

print("Įstatyti sprendiniai:")
check(coef_1, res)

print("\n1b")

coef_lu = np.matrix([[2, 3, 1, 1], [0, 2, 1, 3], [7, -4, 1, 1], [2, -12, 1, 1]]).astype(
    float
)

free_lu = np.matrix([[7, 38, -3.5], [6, 50, -4], [5, 1, -8.5], [9, -53, -1.25]]).astype(
    float
)

P = np.arange(0, np.shape(coef_lu)[0])
res = lu(coef_lu, free_lu, P)
print(res)

print("Įstatyti į sistemą: ")
for vals in res:
    print(check(coef_lu, res))
