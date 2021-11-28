import numpy as np
from matplotlib import pyplot as plt
import math

def main():
    print()
    task1()

    print()
    task2()
    
    print()
    task3()

    print()
    task4()


def solve_diag(a, b, c, d):
    AlphaS = [-c[0] / b[0]]
    BetaS = [d[0] / b[0]]
    GammaS = [b[0]]
    n = len(d)
    result = [0 for i in range(n)]

    for i in range(1, n - 1):
        GammaS.append(b[i] + a[i] * AlphaS[i - 1])
        AlphaS.append(-c[i] / GammaS[i])
        BetaS.append((d[i] - a[i] * BetaS[i - 1]) / GammaS[i])

    GammaS.append(b[n - 1] + a[n - 1] * AlphaS[n - 2])
    BetaS.append((d[n - 1] - a[n - 1] * BetaS[n - 2]) / GammaS[n - 1])

    result[n - 1] = BetaS[n - 1]
    for i in reversed(range(n - 1)):
        result[i] = AlphaS[i] * result[i + 1] + BetaS[i]

    return result

# y'' + qy' + py = f
def get_coeffs(x_0, x_n, h, q, p, f):
    a, b, c, d = [], [], [], []
    for x_i in np.arange(x_0, x_n + h, h):
        a.append(1 / h**2 - q(x_i) / (2*h))        
        b.append(-2 / h**2 + p(x_i))
        c.append(1 / h**2 + q(x_i) / (2*h))
        d.append(f(x_i))
    return a, b, c, d


def task1():
    k = 11
    q = lambda x: 0
    p = lambda x: (1 + math.cos(k) * x**2) / math.sin(k)
    f = lambda x: -1 / math.sin(k)

    x_0, x_n = -1, 1
    y_0, y_n = 0, 0
    n = 10
    h = 1 / n
    e = 1e-3

    print("Задание 1")
    print(f"y'' + y(1 + cos({k})x^2)/sin({k}) = -1/sin({k})")
    print(f"x ∈ [{x_0}, {x_n}]\n")
    print(f"y({x_0}) = {y_0}\ty({x_n}) = {y_n}")

    a, b, c, d = get_coeffs(x_0, x_n, h * 2, q, p, f)
    y1 = solve_diag_2(a, b, c, d, y_n)
    cur_precision = e + 1
    while cur_precision > e:
        y_prev = np.copy(y1)
        a, b, c, d = get_coeffs(x_0, x_n, h, q, p, f)
        y1 = solve_diag_2(a, b, c, d, y_n)
        cur_precision = np.max(np.abs([y1[2*i] - y_prev[i] for i in range(len(y_prev))]))
        print(f"Точность: {cur_precision:.6f}, шаг: {h}")
        h /= 2
    h *= 2
    x1 = np.arange(x_0, x_n + h, h)
    plt.plot(x1, y1)
    plt.show()


def solve_diag_2(a, b, c, d, y_n):
    A, B = [-c[0]/b[0]], [d[0]/b[0]]
    size = len(a)

    for i in range(1, size-1):
        A.append(-c[i] / (b[i] + a[i]*A[i-1]))
        B.append((d[i] - a[i]*B[i-1]) / (b[i] + a[i]*A[i-1]))

    A.append(0)
    B.append((d[size-1] - a[size-1]*B[size-2]) / (b[size-1] +a[size-1]*A[size-2]))

    y = np.zeros(size)
    y[size-1] = y_n
    for i in range(size - 2, -1, -1):
        y[i] = A[i]*y[i+1] + B[i]
    return y

def task2():
    q = lambda x: 0.25 * (1 - x**2)
    p = lambda x: 5 * (1 + (math.cos(x))**2)
    f = lambda x: 15 * math.cos(x)
    x_0, x_n = 0, 2
    y_0, y_n = 0, 4
    n = 10
    h = 1 / n
    e = 2e-2

    print("Задание 2")
    print("u'' + 0.25(1-x^2)u' + 5(1+cos^2(x)) = 15cos(x)")
    print(f"x ∈ [{x_0}, {x_n}]\n")
    print(f"u({x_0}) = {y_0}\nu({x_n}) = {y_n}")

    a, b, c, d = get_coeffs(x_0, x_n, h * 2, q, p, f)
    y1 = solve_diag_2(a, b, c, d, y_n)
    cur_precision = e + 1
    while cur_precision > e:
        y_prev = np.copy(y1)
        a, b, c, d = get_coeffs(x_0, x_n, h, q, p, f)
        y1 = solve_diag_2(a, b, c, d, y_n)
        cur_precision = np.max(np.abs([y1[2*i] - y_prev[i] for i in range(len(y_prev))]))
        print(f"Точность: {cur_precision:.6f}, шаг: {h}")
        h /= 2
    h *= 2
    x1 = np.arange(x_0, x_n + h, h)
    plt.plot(x1, y1)
    plt.show()


def task3():
    A = 1.8
    B = 3.8
    
    n = 4
    eps = 0.1
    
    iteration_count = 0
    
    points = []
    current = []
    prev = []
    h = 0
    
    while True:
        h = (B - A) / n
        points = list(np.linspace(A, B, n + 1))

        del points[0]

        a = [(1 / h ** 2) + 8 * x / (2 * h) for x in points]
        b = [-(2 / h ** 2) -3 for x in points]
        c = [(1 / h ** 2) - (8 * x / (2 * h)) for x in points]
        d = [8 for x in points]

        d[0] = d[0] - a[0] * 1.2

        a[-1] = a[-1] - c[-1] / (h + 3)
        b[-1] = b[-1] + 4 * c[-1] / (h + 3)
        d[-1] = d[-1] - c[-1] * 3.2 * h / (h + 3)
        c[-1] = 0

        print(f'\nКоличество отрезков: {n}')
        print("Шаг равен: ", h)

        current = solve_diag(a, b, c, d)
        # current = solve_diag(a, b, c, d, B)
        if iteration_count != 0 and check_eps(current, prev, eps):
            break

        p = list(points)
        p.insert(0, A)
        c = list(current)
        c.insert(0, 1.2)
        plt.plot(p, c, label=n)

        prev = [item for item in current]
        iteration_count += 1
        n *= 2

    plt.plot(points, current, label=n)
    plt.legend()
    plt.grid()
    plt.show()

def check_eps(current, prev, eps):
    eps_t = max([math.fabs(current[i * 2] - prev[i]) for i in range(len(prev))])
    print(f'Погрешность: {eps_t}')
    if eps_t > eps:
        return False
    return True


def task4():
    A = 0
    B = 1.7
    C = 0.925
    k1 = 1.8
    k2 = 0.4
    q1 = 7.0
    q2 = 12.0

    f = lambda x: 8*x / (2 + x**3)
    
    eps = 1e-3
    n = 16
    iteration_count = 0
    prev = []
    points = []
    current = []
    while True:
        h = (B - A) / n
        points = list(np.linspace(A, B, n))[1:-1]

        a = [-k1 / h ** 2 if x < C else -k2 / h ** 2 for x in points]
        b = [2 * k1 / h ** 2 + q1 if x < C else 2 * k2 / h ** 2 + q2 for x in points]
        c = [-k1 / h ** 2 if x < C else -k2 / h ** 2 for x in points]
        d = [f(x) for x in points]

        current = solve_diag(a, b, c, d)
        
        if iteration_count != 0:
            print(f'\nКоличество отрезков: {n}')
            print("Шаг равен: ", h)

        if iteration_count != 0 and check_eps(current, prev, eps):
            break

        plt.plot(points, current, label=n)

        prev = [item for item in current]
        iteration_count += 1
        n *= 2

    plt.plot(points, current, label=n)
    plt.legend()
    plt.grid()
    plt.show()

if __name__ == "__main__":
    main()
