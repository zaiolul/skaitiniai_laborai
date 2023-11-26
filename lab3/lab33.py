import numpy as np
import matplotlib.pyplot as plt

#Austria
data = [ "80294.93242","78470.02491","78694.10972","83001.77734","84500.38127","89432.43453","90322.80678","90356.77091","88152.69679","85204.26989","84953.81617","77767.40594","83550.95772","81738.65126","78292.56689","78909.03381","75143.18762","76429.62374","76780.91177","78699.13623","75582.15688" ]
scale = 1;

Y = list((float(point) / scale for point in data))
Y = np.array(Y)

X = np.arange(1998, 2019, 1)

def approx(deg):
    plot_points = np.linspace(X[0], X[-1], 200);

    G = np.zeros((len(X), deg+1), dtype=float)

    for i in range(deg+1):
        G[:, i] = np.power(X, i); #fill columns

    cofs = np.linalg.solve(np.dot(np.transpose(G), G), np.dot(np.transpose(G), Y))

    points_y = np.zeros(plot_points.size, dtype=float)
        
    for i in range(deg+1):
        points_y += np.power(plot_points, i) * cofs[i];
  
    plt.plot(plot_points, points_y, 'b')
    for i in range(len(X)):
        plt.plot(X[i], Y[i], 'ro')

    plt.show()



for i in [1, 2, 3, 5]:
    approx(i)
    