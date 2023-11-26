import numpy as np
import matplotlib.pyplot as plt

def func(x):
    return (np.log(x) / (np.sin(2*x) + 1.5)) + (x / 5);

def cebysev_points(n):
    points = []
    for i in range(n):
        points.append(np.cos((np.pi *(2*i+1)) / (2*n)))
    return points

def expand(points, start, end):
    expanded = []
    for point in points:
        expanded.append((start + end) / 2 + (start - end) / 2 * point)
    return expanded

def plot_interpolated(interval,x ,y):
    
    points = np.linspace(interval[0], interval[1], 200);

    plt.plot(points, func(points), 'r')

    for i in range(len(x)):
        plt.plot(x[i], y[i], 'go')
    A = np.zeros((len(x), len(x)), dtype=float);

    for i in range(len(x)):
        A[:,i] = np.power(x, i); #fill columns
    
    print(A)

    cofs = np.linalg.solve(A, y)
    print(cofs)

    points_y = np.zeros(points.size, dtype=float)
    
    for i in range(len(x)):
        points_y += np.power(points, i) * cofs[i];
    
    plt.plot(points, points_y, 'b')
    plt.show()

interval = [2, 10];

num_points = 30;
# step = (float)(interval[1] - interval[0]) / num_points


x_even = np.array(np.linspace(interval[0], interval[1], num_points)) #evely spaced out points
y_even = np.array(list(map(func, x_even)));

x_cebysev = expand(cebysev_points(num_points), interval[0], interval[1]) #evely spaced out points
y_cebysev = np.array(list(map(func, x_cebysev)));

plot_interpolated(interval, x_even, y_even)
plot_interpolated(interval, x_cebysev, y_cebysev)

