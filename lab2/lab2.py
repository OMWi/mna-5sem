import numpy as np
from matplotlib import pyplot as plt
import math

def main():
    print()
    # task1()
    # task2()
    # task3()


# метод прогонки
def solve_diag(a, b, c, d, y_n):
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


# y'' + qy' + py = f
def get_coeffs(x_0, x_n, h, q, p, f):
    a, b, c, d = [], [], [], []
    for x_i in np.arange(x_0, x_n + h, h):
        a.append(2 - h * q(x_i))        
        b.append(-4 + 2 * h**2 * p(x_i))
        c.append(2 + h * q(x_i))
        d.append(f(x_i) * 2 * h**2)
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

    a, b, c, d = get_coeffs(x_0, x_n, h * 2, q, p, f)
    y1 = solve_diag(a, b, c, d, y_n)
    cur_precision = e + 1
    while cur_precision > e:
        y_prev = np.copy(y1)
        a, b, c, d = get_coeffs(x_0, x_n, h, q, p, f)
        y1 = solve_diag(a, b, c, d, y_n)
        cur_precision = np.max(np.abs([y1[2*i] - y_prev[i] for i in range(len(y_prev))]))
        print(f"Точность: {cur_precision:.6f}, шаг: {h}")
        h /= 2

    h *= 2
    x1 = np.arange(x_0, x_n + h, h)
    plt.plot(x1, y1)
    plt.show()

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
    print(f"u({x_0}) = {y_0}\nu(x_2) = {y_n}")

    a, b, c, d = get_coeffs(x_0, x_n, h * 2, q, p, f)
    y1 = solve_diag(a, b, c, d, y_n)
    cur_precision = e + 1
    while cur_precision > e:
        y_prev = np.copy(y1)
        a, b, c, d = get_coeffs(x_0, x_n, h, q, p, f)
        y1 = solve_diag(a, b, c, d, y_n)
        cur_precision = np.max(np.abs([y1[2*i] - y_prev[i] for i in range(len(y_prev))]))
        print(f"Точность: {cur_precision:.6f}, шаг: {h}")
        h /= 2

    h *= 2
    x1 = np.arange(x_0, x_n + h, h)
    plt.plot(x1, y1)
    plt.show()


# def task3():
#     q = lambda x: -3
#     p = lambda x: 8*x
#     f = lambda x: 8
#     x_0, x_n = 
#     y_0, y_n = 
#     n = 10
#     h = 1 / n
#     e = 1e-1

#     print("Задание 3")
#     print("u'' - 3u' + 8xu = 8")
#     print("u(1.8) - 0.5u'(1.8) = 2")
#     print("u(3.8) = 5")
#     print(f"x ∈ [{x_0}, {x_n}]\n")
#     print(f"u({x_0}) = {y_0}\nu(x_2) = {y_n}")

#     a, b, c, d = get_coeffs(x_0, x_n, h * 2, q, p, f)
#     y1 = solve_diag(a, b, c, d, y_n)
#     cur_precision = e + 1
#     while cur_precision > e:
#         y_prev = np.copy(y1)
#         a, b, c, d = get_coeffs(x_0, x_n, h, q, p, f)
#         y1 = solve_diag(a, b, c, d, y_n)
#         cur_precision = np.max(np.abs([y1[2*i] - y_prev[i] for i in range(len(y_prev))]))
#         print(f"Точность: {cur_precision:.6f}, шаг: {h}")
#         h /= 2

#     h *= 2
#     x1 = np.arange(x_0, x_n + h, h)
#     plt.plot(x1, y1)
#     plt.show()








if __name__ == "__main__":
    main()