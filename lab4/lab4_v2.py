import sympy
import numpy
from matplotlib import pyplot


x = sympy.Symbol('x')
start = 0
end = 1
y_arrays = []
k = 1
T = 0.05
Y2 = []
g1 = lambda t: 0.5
g2 = lambda t: 0.5
fi = lambda x: 0
data_table = []


def solve(tau, h, t, last_layer = None):
    fx = 1
    A = []
    B = []
    Y = []
    i = 0
    if (True):
        pass
    xi = start
    while xi < end - 1e-5:
        if t != 0:
            a = -k * tau
            b = h**2 + 2 * k * tau
            c = -k * tau
            d = tau * h**2 * fx + h**2 * last_layer[i]

            if i == 0:
                Ai = - c / b
                Bi = (d - a * g1(t)) / b
            else:
                Ai = -c / (b + a * A[-1])
                Bi = (d - a * B[-1]) / (b + a * A[-1])

            A.append(Ai)
            B.append(Bi)
        else:
            Y.append(fi(xi))
        xi += h
        i += 1

    if t == 0:
        Y.append(end ** 2)
        return Y

    Y = [0] * 21

    Y[0] = g1(t)
    Y[-2] = (A[-2] * h * g2(t) + B[-2]) / (1 - A[-2])
    Y[-1] = Y[-2] + 4*h  # yb

    while i > 1:
        i -= 1
        Y[i] = A[i] * Y[i + 1] + B[i]

    return Y


def solve2(tau, h, t, j):
    Y = [0] * 21
    fx = 1
    xi = start

    if j == 0:
        i = 0
        while xi < end:
            Y[i] = fi(xi)
            i += 1
            xi += h
    else:
        Y[0] = g1(t)
        for i in range(1, 10):
            Y[i] = Y2[j - 1][i] + tau / h ** 2 * (Y2[j - 1][i + 1] - 2 * Y2[j - 1][i] + Y2[j - 1][i - 1]) #+ tau * fx
            xi += h

        Y[-2] = Y2[j - 1][-2] + tau / h ** 2 * (2 * h * g2(t) - Y2[j - 1][-2] + Y2[j - 1][-3]) #+ tau * fx
        Y[-1] = Y[-2] + 2 * h *g2(t)

    Y2.append(Y)


def implicit_function(h, tau):
    y_arrays.clear()
    t = 0
    for i in range(0, 10):
        if i == 0:
            y_arrays.append(solve(tau, h, t))
        else:
            y_arrays.append(solve(tau, h, t, y_arrays[-1]))
        if t == 0:
            t = 2 * tau
        else:
            t *= 2
        if t > T:
            break


def explicit_function(h, tau):
    t = 0
    Y2.clear()
    for i in range(0, 10):
        solve2(tau, h, t, i)
        if t == 0:
            t = 2 * tau
        else:
            t *= 2
        if t > T:
            break


def show_plot(array, h, method, t):
    x = numpy.arange(start, end + h, h)
    pyplot.plot(x, array, label=f't={t}')
    pyplot.legend()
    pyplot.title(method + " метод ")
  

def main():
    h = (end - start) / 20
    tau = 0.5 * h ** 2 / k
    implicit_function(h, tau)
    t = 0
    for array in y_arrays:
        show_plot(array, h, "Неявный", "{t:.5f}".format(t=t))
        if t == 0:
            t = 2 * tau
        else:
            t *= 2
    pyplot.show()

    explicit_function(h, tau)
    t = 0
    for array in Y2:
        show_plot(array, h, "Явный", "{t:.5f}".format(t=t))
        if t == 0:
            t = 2 * tau
        else:
            t *= 2
    pyplot.show()

    tau = h ** 2 / 6
    explicit_function(h, tau)
    t = 0
    for array in Y2:
        show_plot(array, h, "Явный t = (h^2)/6", "{t:.5f}".format(t=t))
        if t == 0:
            t = 2 * tau
        else:
            t *= 2
    pyplot.show()


if __name__ == "__main__":
    main()

