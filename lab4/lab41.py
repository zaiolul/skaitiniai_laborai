#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate 

def draw(x, y, xlabel="x", ylabel="y", title="sample text", color="r", style=""):
    plt.plot(x, y, color + style)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()

def _func(t, velocity):
    if t <= 25:
        return (-9.8 + 0.1 * velocity ** 2 / (60 + 15)) 
    else:
        return (-9.8 + 7 * velocity ** 2 / (60 + 15)) 
    
def func(velocity, t, dt):
    return _func(t, velocity) * dt

def eulers(height, velocity, max_time, dt, plot_results=True):
    time_intervals = int(max_time / dt)
    times = np.linspace(0, max_time, time_intervals)
    # dt = times[1]
    results = np.zeros(shape=[2, int(time_intervals)], dtype=np.float64)

    results[0, 0] = height
    results[1, 0] = velocity

    parachute = False

    for i in range(1, len(times)):
        dv = func(results[1, i - 1],times[i], dt)

        results[0, i] = results[0, i - 1] + results[1, i - 1] * dt
        results[1, i] = results[1, i - 1] + dv
        if not parachute and times[i] > 25:
            parachute = True
            print("Parašiutas išskleistas {0} m aukštyje".format(results[0,i]))
        if results[0, i] <= 0:
            print("Pasiekė žeme t = {0:.02f} s, v = {1:.02f} m/s".format(times[i], results[1,i]))
            break
    
    if plot_results:
        draw(times, results[0,:],"laikas, s", "aukštis, m", "Eulerio m. aukščio kitimas (dt={})".format(dt), "r")
        draw(times, results[1,:],"laikas, s", "greitis, m/s", "Eulerio m. greičio kitimas (dt={})".format(dt), "r")
    return [results,times]


def ivrk(height, velocity, max_time, dt, plot_results=True):
    time_intervals = int(max_time / dt)
    times = np.linspace(0, max_time, time_intervals)

    results = np.zeros(shape=[2, time_intervals], dtype=np.float64)
    
    results[:, 0] = [height, velocity]
    parachute = False
    for i in range(1, len(times)):
        dv1 = func(results[1, i - 1], times[i] - dt / 2, dt / 2)
      
        v1 = results[1, i - 1] + dv1
        
        dv2 = func(v1, times[i] - dt / 2, dt / 2)
        v2 = v1 + dv2

        dv3 = func(v2, times[i] - dt / 2, dt)
        v3 = v2 + dv3
        results[0, i] = results[0, i - 1] + results[1, i - 1] * dt
        results[1, i] = results[1, i - 1] + (dt / 6) * \
            (       func(results[1, i - 1],times[i], 1) + \
                2 * func(v1, times[i] - dt / 2, 1) + \
                2 * func(v2, times[i] - dt / 2, 1) + \
                    func(v3, times[i] - dt, 1))
        
        if not parachute and times[i] > 25:
            parachute = True
            print("Parašiutas išskleistas {0:.02f} m aukštyje".format(results[0,i]))
        if results[0, i] <= 0:
            print("Pasiekė žeme t = {0:.02f} s, v = {1:.02f} m/s".format(times[i], results[1,i]))
            break
    if plot_results:
        draw(times, results[0,:],"laikas, s", "aukštis, m", "IV RK m. aukščio kitimas", "r")
        draw(times, results[1,:],"laikas, s", "greitis, m/s", "IV RK m. greičio kitimas", "r")
    return [results, times]
  
def stability(function, height, velocity, max_time, dt):
    cur_dt = dt
    while True:
        res = function(height, velocity, max_time, cur_dt, plot_results=False)
        if res[0][1,-1] > 1e9:
            print("Gautas perpildymas, didžiausias žingsnis: {0:.02f}".format(cur_dt))
            break
        cur_dt += 0.01
        plt.plot(res[1], res[0][1])
        plt.title("Greičio kitimai su skirtingais žingsniais")
        plt.xlabel("laikas, s")
        plt.ylabel("greitis, m/s")
        plt.draw()
    plt.show()

def compare(height, velocity, dt, max_time):
    times = np.linspace(0, max_time, int(max_time / dt))
    result_eulers = eulers(height, velocity, max_time, dt, False)
    result_ivrk = ivrk(height, velocity, max_time, dt, False)

    result_integrate = scipy.integrate.solve_ivp(   \
        fun=_func,  \
        t_span=[0,300],     \
        t_eval=times,    \
        y0=[velocity])
    
    plt.plot(times, result_eulers[0][1], 'r')
    plt.plot(times, result_ivrk[0][1], 'g')
    plt.plot(times, result_integrate.y[0], 'b')
    plt.legend(["Eulerio", "IV RK", "SciPy.Integrate"])
    plt.title("Sprendinių palyginimas")
    plt.xlabel("laikas, s")
    plt.ylabel("greitis, m/s")
    plt.show()

if __name__ == "__main__":
    #pradiniai duomenys
    height = 3500.0
    velocity = 0.0
    max_time = 300 #240 s
    dt = 0.01
    # eulers(height, velocity, max_time, dt, plot_results=True)
    # ivrk(height, velocity, max_time,  dt, plot_results=True)
    # stability(eulers, height, velocity, max_time, dt)
    stability(ivrk, height, velocity, max_time, dt)

    # compare(height, velocity, dt, max_time)
    # for intervals in [10_000, 8_000, 6_000, 4_000, 3_000, 1500]:     
    #     a = eulers(height, velocity, max_time, intervals )
    #     b = ivrk(height, velocity, max_time, intervals)
   
    #     plt.plot(a[1], a[0][0], 'r')
    #     plt.plot(b[1], b[0][0], 'b--')
    #     plt.legend(["Euler", "IV RK"])
    #     plt.title("Aukščio grafikas, kai žingsnis= {0:0.02f} ".format(max_time / intervals))
    #     plt.show()

    
    #     plt.plot(a[1], a[0][1], 'r')
    #     plt.plot(b[1], b[0][1], 'b--')
    #     plt.title("Greičio grafikas, kai žingsnis= {0:0.02f} ".format(max_time / intervals))
    #     plt.legend(["Euler", "IV RK"])
        # plt.show())

 


    
    
