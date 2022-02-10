import numpy as np
import sympy as sp
from matplotlib import pyplot as plt
import math
from scipy import integrate

def main():
    task1()

    a, b = 1, 1.5
    ua, ub = 2, 1
    k2_values = [
        [1, 100],
        [100, 1],
    ]
    k3_values = [
        [1, 50, 100],
        [100, 50, 1],
        [50, 100, 50],
        [100, 5, 100],
    ]
    # 4a
    for k1, k2 in k2_values:
        x, y = task2(a, b, ua, ub, None, k1, k2)
        show_plot(x, y, f"k: {k1}, {k2}")
    # 4b
    for k1, k2, k3 in k3_values:
        x, y = task2(a, b, ua, ub, None, k1, k2, k3)
        show_plot(x, y, f"k: {k1}, {k2}, {k3}")
    # 5
    sources = [
        [(20, a + (b - a) / 2)],
        [(30, a + (b - a) / 4), (30, a + 3 * (b - a) / 4)],
        [(100, a + (b - a) / 4), (20, a + 3 * (b - a) / 4)],
    ]
    for k1, k2 in k2_values:
        for src in sources:
            x, y = task2(a, b, ua, ub, src, k1, k2)
            show_plot(x, y, f'k: {k1}, {k2}\nsources:{src}')
    for k1, k2, k3 in k3_values:
        for src in sources:
            x, y = task2(a, b, ua, ub, src, k1, k2, k3)
            show_plot(x, y, f'k: {k1}, {k2}, {k3}\nsources:{src}')
            

    a, b = 1.0, 1.5
    ua, ub = 2.0, 1.0
    T = 1

    fi = lambda x: x ** -3
    f = lambda x, t: (10 * math.cos(x)) * (1 - math.e ** (-t))
    k = lambda x: math.cos(x)
    dk = lambda x: 0

    n = 10
    h = 0.1
    tau = 0.001 
    m = int(T / tau)


    vector_x = np.linspace(a, b, n + 1)
    vector_t = np.linspace(0, T, m + 1)

    vector_x = np.linspace(a, b, n + 1)
    vector_t = np.linspace(0, T, m + 1)
    matrix_u = np.zeros((m + 1, n + 1))

    for i in range(m + 1):
        matrix_u[i][0] = ua
        matrix_u[i][n] = ub

    for i in range(n + 1):
        matrix_u[0][i] = fi(vector_x[i])
                
    for i in range(1, m + 1):
        for j in range(1, n):
            matrix_u[i][j] = (
            tau * dk(vector_x[j]) / h * (matrix_u[i - 1][j + 1] - matrix_u[i - 1][j]) +
            tau * k(vector_x[j]) / (h ** 2) * (matrix_u[i - 1][j + 1] - 2 * matrix_u[i - 1][j] + 
            matrix_u[i - 1][j - 1]) + tau * f(vector_x[j], vector_t[i - 1]) + matrix_u[i - 1][j]
            )
    for i in [0, 5, 20, 200]:
        plt.plot(vector_x, matrix_u[i], label=f'{i}t')
    plt.legend()
    plt.grid()
    plt.show()


    a, b = -1, 1
    ua, ub = 0, 5
    k = 2
    t = 0.1
    h = (b - a) / 10
    dt = 0.5 * h**2 / k
    x_row, t_row, A = task4(a, b, ua, ub, t, dt, k)
    show_plots_2(x_row, A, t, dt)


def task1():
    x, k, f = sp.symbols('x k f')
    f = sp.cos(x) * 10
    k = sp.cos(x)
    a, b = 1, 1.5
    ua, ub = 2, 1
    for var in range(1, 4):
        print('Набор #{}'.format(var))
        c, cur_k, cur_ua, cur_ub = get_params(k, ua, ub, var)
        u = solve(cur_k, a, b, cur_ua, cur_ub, f)
        plot_solution(u, a, b, 'k={}'.format(var))

    plt.legend()
    plt.show()

    for var in [1, 4]:
        print('Набор #{}'.format(var))
        c, cur_k, cur_ua, cur_ub = get_params(k, ua, ub, var)
        u = solve(cur_k, a, b, cur_ua, cur_ub, f)
        plot_solution(u, a, b, 'k={}'.format(var))

    plt.legend()
    plt.show()

    for var in [5, 6, 7]:
        print('Набор #{}'.format(var))
        c, cur_k, cur_ua, cur_ub = get_params(k, ua, ub, var)
        u = solve(cur_k, a, b, cur_ua, cur_ub, f)
        plot_solution(u, a, b, 'k={}'.format(var))

    plt.legend()
    plt.show()

def plot_solution(foo, x1, x2, label):
    x = sp.symbols("x")
    # TODO: change 200 param
    args = np.linspace(x1, x2, 100)
    y = [foo.subs({x: t}) for t in args]
    plt.plot(args, y, label=label)

# TODO change output
def solve(k, a, b, ua, ub, f):
    x, c1, c2 = sp.symbols("x,c1,c2")
    u = -sp.integrate(sp.integrate(f, x) / k, x) + c1 * x + c2
    sol = sp.solve([u.subs({x: a}) - ua, u.subs({x: b}) - ub], (c1, c2))
    print(f'k={k}\tf={f}\na={a}\tb={b}\tua={ua}\tub={ub}')
    print(u)
    print(f'Решение: c1 = {sol[c1]:.5f}\t c2 = {sol[c2]:.5f}')
    print()    
    return u.subs({c1: sol[c1], c2: sol[c2]})

def get_params(k, ua, ub, variant):
    params = [
        (1, k, ua, ub),
        (2, 2*k, ua, ub),
        (0.1, 0.1*k, ua, ub),
        (1, 1/k, ua, ub),
        (1, k, -ua, ub),
        (1, k, ua, -ub), 
        (1, k, -ua, -ub)
    ]
    return params[variant-1]


def task2(a, b, ua, ub, sources, *k_values):
    f = lambda x: 10 * math.cos(x)    
    h = (b - a) / 150
    n = int((b - a) / h) + 1    
    A = np.zeros((n, n))
    g = np.zeros((n, 1))
    
    A[0, 0] = A[-1, -1] = 1
    g[0] = ua
    g[-1] = ub
    x_values = np.linspace(a, b, n)

    for i in range(1, n - 1):
        A[i, i-1] = get_coeff(a, b, x_values[i-1], x_values[i], *k_values) 
        A[i, i] = (-get_coeff(a, b, x_values[i-1], x_values[i], *k_values) 
                   - get_coeff(a, b, x_values[i], x_values[i+1], *k_values))
        A[i, i+1] =  get_coeff(a, b, x_values[i], x_values[i+1],  *k_values)
        if sources:
            g[i] = -sum((get_delta_integral(x_values[i], x0, ci, h) for ci, x0 in sources))
        else:
            g[i] = -get_phi(f, x_values[i] - h/2, x_values[i] + h/2)
    return x_values, np.linalg.solve(A, g)

def show_plot(x, y, label):
    plt.plot(x, y, label=label)
    plt.grid()
    plt.legend()
    plt.xlabel('x', size=14)
    plt.ylabel('y(x)', size=14)
    plt.show()

def get_k(a, b, x, k1, k2, k3=None):
    if k3 is None:
        if a <= x <= (b + a) * 0.5:
            return k1
        return k2
    else:
        if a <= x <= (a + (b - a) / 3):
            return k1
        elif (a + 2 * (b - a) / 3) < x <= b:
            return k3
        return k2   

def get_coeff(a, b, x1, x2, *k_values):
    return (integrate.quad(lambda x: 1 / get_k(a, b, x, *k_values), x1, x2)[0]) ** (-1)

def get_phi(f, x1, x2):
    return integrate.quad(f, x1, x2)[0]

def get_delta_integral(x, x0, c, h):
    if abs(x - x0) - h / 2 < 1e-5:
        return c / 2
    elif abs(x - x0) < h/2:
        return c
    return 0


def task3(phi, a, b, ua, ub, t, dt, h):
    k = lambda x: math.cos(x)
    f = lambda x: 10 * math.cos(x)

    nt = int(t / dt) + 1
    nx = int((b - a) / h) + 1
    x_values = np.linspace(a, b, nx)
    t_values = np.linspace(0, t, nt)
    A = np.zeros((nt, nx))
    A[:, 0] = ua
    A[:, -1] = ub
    A[0] = [phi(x) for x in x_values]
    
    for i in range(nt - 1):
        for j in range(1, nx - 1):            
            A[i + 1, j] = dt * k(x_values[j] - h / 2) / (h ** 2) * A[i, j - 1] + \
                          (1 - dt * (k(x_values[j] - h / 2) + k(x_values[j] + h / 2)) / (h ** 2)) * A[i, j] + \
                          dt * k(x_values[j] + h / 2) / (h ** 2) * A[i, j + 1] + \
                          dt * f(x_values[j])*(1 - math.exp(-t_values[i]))
    return x_values, t_values, A

def task4(a, b, ua, ub, t, dt, k):
    f = lambda x, t: 0
    phi = lambda x: x**2
    h = (b - a) / 10

    nt = int(t / dt) + 1
    nx = int((b - a) / h) + 1
    x_values = np.linspace(a, b, nx)
    t_values = np.linspace(0, t, nt)
    A = np.zeros((nt, nx))
    A[0] = [phi(x) for x in x_values]
    A[:, 0] = ua
    A[:, -1] = ub
    
    for i in range(nt - 1):
        for j in range(1, nx - 1):            
            A[i + 1, j] = (k * dt / h ** 2 * A[i, j-1] + (1 - 2 * k * dt / h ** 2) * A[i, j] 
                           + k * dt / h ** 2* A[i, j+1] + dt * f(x_values[j], t_values[i]))
    return x_values, t_values, A

def show_plots_2(x, A, t, dt):
    for i in range(0, int(t / dt) + 1, 2):
        plt.plot(x, A[i], label=fr't={i}$\tau$')
    plt.grid()
    plt.legend()
    plt.xlabel('x', size=14)
    plt.ylabel('u(x, t)', size=14)
    plt.show()

if __name__ == "__main__":
    main()
