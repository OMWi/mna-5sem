import numpy as np
import math
from matplotlib import pyplot as plt
from scipy import integrate
import pylab

a, b = 0, 1
k = 1   
T = 0.05
phi = lambda x: abs(x - 0.5)
g_1 = lambda t: 0.5
g_2 = lambda t: 0.5
f = lambda x, t: 0

def main():
    N = 100
    h = (b - a) / N
    dt = h ** 2 / 6

    x_row, t_row, A = explicit_1(h, dt)
    show_plots(x_row, A)

    x_row, t_row, A = explicit_2(h, dt)
    print(len(x_row))
    show_plots(x_row, A)

    x_row, t_row, A = implicit_1(h, dt)
    show_plots(x_row, A)

    N = 100
    h = (b - a) / N
    dt = (((b - a) / N) ** 2) / 6
    x_row, t_row, A = implicit_2(h, dt)
    show_plots(x_row, A)


def explicit_1(h, dt):
    nt = int(T / dt) + 1
    nx = int((b - a) / h) + 1
    x_values = np.linspace(a, b, nx)
    t_values = np.linspace(0, dt, nt)
    A = np.zeros((nt, nx))
    A[0] = [phi(x) for x in x_values]
    A[:, 0] = [g_1(t) for t in t_values]
    
    for i in range(nt - 1):
        for j in range(1, nx - 1):            
            A[i + 1, j] = k * dt / h ** 2 * A[i, j-1] + (1 - 2 * k * dt / h ** 2) * A[i, j] + \
                          k * dt / h ** 2 * A[i, j + 1] + dt * f(x_values[j], t_values[i])
        A[i + 1, -1] = A[i + 1, -2] + h * g_2(t_values[i])
    return x_values, t_values, A

def show_plots(x, A):
    for i in range(0, len(A), 100):
        plt.plot(x, A[i])
    plt.grid()
    plt.xlabel('x', size=14)
    plt.ylabel('u(x, t)', size=14)
    plt.show()

def explicit_2(h, dt):
    nt = int(T / dt) + 1
    nx = int((b - a) / h) + 1
    x_values = np.linspace(a, b, nx)
    t_values = np.linspace(0, dt, nt)
    A = np.zeros((nt, nx))
    A[0] = [phi(x) for x in x_values]
    A[:, 0] = [g_1(t) for t in t_values]
    
    for i in range(nt - 1):
        for j in range(1, nx - 1):            
            A[i + 1, j] = k * dt / h ** 2 * A[i, j - 1] + (1 - 2 * k * dt / h ** 2) * A[i, j] + \
                          k * dt / h ** 2 * A[i, j + 1] + dt * f(x_values[j], t_values[i])
        A[i + 1, -1] = k * dt / h ** 2 * A[i, -2] + (1 - 2 * k * dt / h ** 2) * A[i, -1] + k * dt / h ** 2 * \
                        (2 * h * g_2(t_values[i]) + A[i + 1, -2]) + dt * f(x_values[j], t_values[i])
    return x_values, t_values, A

def implicit_1(h, dt):
    nt = int(T / dt) + 1
    nx = int((b - a) / h) + 1
    x_values = np.linspace(a, b, nx)
    t_values = np.linspace(0, dt, nt)
    A = np.zeros((nt, nx))
    A[0] = [phi(x) for x in x_values]
    A[:, 0] = [g_1(t) for t in t_values]
    
    matrix = np.zeros((nx, nx))
    for j in range(1, nx - 1):
            matrix[j, j-1] = - k * dt / h ** 2 
            matrix[j, j] = 1 + 2 * k * dt / h ** 2
            matrix[j, j+1] = - k * dt / h ** 2
    matrix[0, 0] = 1
    matrix[-1, -1], matrix[-1, -2] = 1, -1 
    
    for i in range(1, nt):
        free_members = np.array([dt * f(x, t_values[i] + dt) for x in x_values]) + A[i - 1]
        free_members[0] = g_1(t_values[i])
        free_members[-1] = h * g_2(t_values[i])
        A[i] = np.linalg.solve(matrix, free_members)
    return x_values, t_values, A

def implicit_2(h, dt):
    nt = int(T / dt) + 1
    nx = int((b - a) / h) + 1
    x_values = np.linspace(a, b, nx)
    t_values = np.linspace(0, dt, nt)
    A = np.zeros((nt, nx))
    A[0] = [phi(x) for x in x_values]
    A[:, 0] = [g_1(t) for t in t_values]
    
    matrix = np.zeros((nx + 1, nx + 1))
    for j in range(1, nx):
        matrix[j, j-1] = - k * dt / h ** 2 
        matrix[j, j] = 1 + 2 * k * dt / h ** 2
        matrix[j, j+1] = - k * dt / h ** 2
    matrix[0, 0] = 1
    matrix[-1, -1], matrix[-1, -3] = 1, -1 
    
    for i in range(1, nt):
        free_members = np.zeros(nx + 1)
        free_members[1:-1] = np.array([dt * f(x, t_values[i] + dt) for x in x_values[1:]]) + A[i - 1, 1:]
        free_members[0] = g_1(t_values[i])
        free_members[-1] = 2 * h * g_2(t_values[i])
        A[i] = np.linalg.solve(matrix, free_members)[:-1]
    return x_values, t_values, A

if __name__ == "__main__":
    main()
    