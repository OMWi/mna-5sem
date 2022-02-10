import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

def main():
    du = 0.2e-2
    L = 18e-2
    E = 120e9
    rho = 5.9e3
    T = 0.001
    phi = lambda x: (-4 * du / L ** 2) * x ** 2 + (4 * du / L) * x
    nx = 100
    nt = 10000
    
    A = task1(T, L, E, nx, nt, rho, phi)

    x_values = np.linspace(0, L, nx)
    for i in range(0, nt, int(nt / 20)):
        draw(x_values, A[i], 'x', 'u(x, t)')
    plt.grid()
    plt.show()



def task1(T, L, E, nx, nt, rho, phi):
    dt = T / nt
    dx = L / nx
    C = np.sqrt(E / rho) * dt / dx
    if C > 1:
        print('Не выполнено условие сходимости.')
        return
    
    x_values = np.linspace(0, L, nx)
    
    A = np.zeros((nt, nx))
    A[0] = [phi(x) for x in x_values]
    A[1, 1:-1] = [phi(x_values[i]) * (1 - dt ** 2 / 2) for i in range(1, nx - 1)]

    for i in range(1, nt-1):
        A[i + 1, 1:-1] = C ** 2 * (A[i, 2:] - 2 * A[i, 1:-1] + A[i, :-2]) + 2 * A[i, 1:-1] - A[i-1, 1:-1]
    return A

def draw(x, y, x_label, y_label):
    plt.plot(x, y)
    plt.xlabel(x_label, size=14)
    plt.ylabel(y_label, size=14)
    return plt

if __name__ == "__main__":
    main()