import numpy as np
import matplotlib.pyplot as plt

def dx_lang_2nd(x, x_start, y_start, x_cur, y_cur, x_end, y_end):
    return (((x - x_cur) + (x - x_end)) / ((x_start - x_cur) * (x_start - x_end))) * y_start + \
          (((x - x_start) + (x - x_end)) / ((x_cur - x_start) * (x_cur - x_end))) * y_cur + \
          (((x - x_start) + (x - x_cur)) / ((x_end - x_start) * (x_end - x_cur))) * y_end

def akima(X, Y):
    result = []
    n = len(X);
    for i in range(n):
        if i == 0:
            result.append(dx_lang_2nd(X[0], X[0], Y[0], X[1], Y[1], X[2], Y[2]))
        if i == n - 1:
            result.append(dx_lang_2nd(X[n - 1], X[n - 3], Y[n - 3], X[n - 2], Y[n - 2], X[n - 1], Y[n - 1]))
        else: 
            result.append(dx_lang_2nd(X[i], X[i - 1], Y[i - 1], X[i], Y[i], X[i + 1], Y[i + 1]))
    return result

#Austria
data = [ "80294.93242","78470.02491","78694.10972","83001.77734","84500.38127","89432.43453","90322.80678","90356.77091","88152.69679","85204.26989","84953.81617","77767.40594","83550.95772","81738.65126","78292.56689","78909.03381","75143.18762","76429.62374","76780.91177","78699.13623","75582.15688" ]
scale = 1;
Y = list((float(point) / scale for point in data))


# X = np.arange(1998, 2018, 1);
X = np.arange(1998, 2019, 1);
dY = akima(X, Y)

# for i in range(len(X)):
#     print(lagrange(X, X[i], i));

draw_points = 100;
plt.plot(X[0], Y[0], 'ro')
for i in range(len(X) - 1):
    plot_x = np.linspace(X[i], X[i + 1], draw_points);
    
    plot_y = []
    plt.plot(X[i+1], Y[i+1], 'ro')
    d = X[i+1] - X[i]
    for k in range(len(plot_x)):
        s = plot_x[k] - X[i];
        U1 = 1 - 3 * (s**2 / d**2) + 2 * (s**3 / d**3)
        V1 = s - 2 * (s**2 / d) + (s**3 / d**2)
        U2 = 3 * (s ** 2 / d **2) - 2 * (s**3 / d ** 3)
        V2 = -1 * (s ** 2 / d) + (s ** 3 / d ** 2)

        f = U1*Y[i] + V1*dY[i]
        f += U2*Y[i + 1] + V2*dY[i + 1]
        plot_y.append(f);
    # print("Plotting ", X[i], " - ", X[i+1])
    plt.plot(plot_x, plot_y, 'b')
    plt.draw()


plt.show()
        

