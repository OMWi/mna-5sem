import numpy as np
import sympy as sp
from matplotlib import pyplot as plt
from matplotlib import animation
import math
from sympy.functions.special.delta_functions import DiracDelta

x, k, f = sp.symbols('x k f')
f = sp.cos(x) * 10
k = sp.cos(x)
a, b = 1, 1.5
ua, ub = 2, 1

def main():    
    # task1()
    # task2()
    task3()


def task1():
    for var in range(1, 4):
        print('Solving problem set #{}'.format(var))
        c, cur_k, cur_ua, cur_ub = get_params(k, ua, ub, var)
        u = solve(cur_k, a, b, cur_ua, cur_ub, f)
        plot_solution(u, a, b, 'k={}'.format(var))

    plt.legend()
    plt.show()

    for var in [1, 4]:
        print('Solving problem set #{}'.format(var))
        c, cur_k, cur_ua, cur_ub = get_params(k, ua, ub, var)
        u = solve(cur_k, a, b, cur_ua, cur_ub, f)
        plot_solution(u, a, b, 'k={}'.format(var))

    plt.legend()
    plt.show()

    for var in [5, 6, 7]:
        print('Solving problem set #{}'.format(var))
        c, cur_k, cur_ua, cur_ub = get_params(k, ua, ub, var)
        u = solve(cur_k, a, b, cur_ua, cur_ub, f)
        plot_solution(u, a, b, 'k={}'.format(var))

    plt.legend()
    plt.show()


def task2():
    # TODO change k11, k12, k21, k22 ...
    task2_4aa()
    task2_4ab()
    task2_4ba()
    task2_4bb()
    task2_4bc()
    task2_4bd()
    f1 = DiracDelta(x - (a + b) / 2)
    task2_5a(f1)
    mid = (a + b) / 2
    x1 = (a + mid) / 2
    x2 = (mid + b) / 2
    f2 = DiracDelta(x - x1)
    f3 = DiracDelta(x - x2)
    task2_5b(f2, f3)
    task2_5c(f2, f3)
    task2_5d()
    
def task2_4aa():
    k11 = x**-2
    u = solve(k11, a, b, ua, ub, f)
    plot_solution(u, a, (b + a) / 2, 'k={}'.format(k11))
    k12 = x
    u = solve(k12, a, b, ua, ub, f)
    plot_solution(u, (b + a) / 2, b, 'k={}'.format(k12))
    plt.legend()
    plt.show()

def task2_4ab():
    k21 = sp.exp(2 * x)
    u = solve(k21, a, b, ua, ub, f)
    plot_solution(u, a, (b + a) / 2, 'k={}'.format(k21))
    k22 = sp.exp(-2 * x)
    u = solve(k22, a, b, ua, ub, f)
    plot_solution(u, (b + a) / 2, b, 'k={}'.format(k22))
    plt.legend()
    plt.show()

def task2_4ba():
    k31 = 0.34 * x**2
    u = solve(k31, a, b, ua, ub, f)
    plot_solution(u, a, a + (b - a) / 3, 'k={}'.format(k31))
    k32 = 2 * x**2
    u = solve(k32, a, b, ua, ub, f)
    plot_solution(u, a + (b - a) / 3, b - (b - a) / 3, 'k={}'.format(k32))
    k33 = sp.exp(3) * x**2
    u = solve(k33, a, b, ua, ub, f)
    plot_solution(u, b - (b - a) / 3, b, 'k={}'.format(k33))
    plt.legend()
    plt.show()

def task2_4bb():
    k31 = sp.exp(3) * x**2
    u = solve(k31, a, b, ua, ub, f)
    plot_solution(u, a, a + (b - a) / 3, 'k={}'.format(k31))
    k32 = 2 * x**2
    u = solve(k32, a, b, ua, ub, f)
    plot_solution(u, a + (b - a) / 3, b - (b - a) / 3, 'k={}'.format(k32))
    k33 = 0.5 * x**2
    u = solve(k33, a, b, ua, ub, f)
    plot_solution(u, b - (b - a) / 3, b, 'k={}'.format(k33))
    plt.legend()
    plt.show()

def task2_4bc():
    k31 = sp.exp(0.1 * x)
    u = solve(k31, a, b, ua, ub, f)
    plot_solution(u, a, a + (b - a) / 3, 'k={}'.format(k31))
    k32 = 2 * sp.exp(0.1 * x)
    u = solve(k32, a, b, ua, ub, f)
    plot_solution(u, a + (b - a) / 3, b - (b - a) / 3, 'k={}'.format(k32))
    k33 = sp.exp(0.1 * x)
    u = solve(k33, a, b, ua, ub, f)
    plot_solution(u, b - (b - a) / 3, b, 'k={}'.format(k33))
    plt.legend()
    plt.show()

def task2_4bd():
    k31 = 20 * sp.exp(0.1 * x)
    u = solve(k31, a, b, ua, ub, f)
    plot_solution(u, a, a + (b - a) / 3, 'k={}'.format(k31))
    k32 = sp.exp(0.1 * x)
    u = solve(k32, a, b, ua, ub, f)
    plot_solution(u, a + (b - a) / 3, b - (b - a) / 3, 'k={}'.format(k32))
    k33 = 20 * sp.exp(0.1 * x)
    u = solve(k33, a, b, ua, ub, f)
    plot_solution(u, b - (b - a) / 3, b, 'k={}'.format(k33))
    plt.legend()
    plt.show()

def task2_5a(f1):    
    u = solve(0.01 * k, a, b, ua, ub, f1)
    plot_solution(u, a, b, 'f = {}'.format(f1))
    plt.legend()
    plt.show()
    us = [small_solve(1 / i * k, a, b, ua, ub, f1) for i in range(1, 51)];
    # skipped some code 

def task2_5b(f2, f3):    
    u = solve(1/20/k, a, b, ua, ub, f2 + f3)
    plot_solution(u, a, b, 'f2 + f3 = {f2} + {f3}'.format(f2=f2, f3=f3))
    plt.legend()
    plt.show()

def task2_5c(f2, f3):
    u = solve(1/20/k, a, b, ua, ub, 6 * f2 + f3)
    plot_solution(u, a, b, 'f2 > f3 = {f2} + {f3}'.format(f2=6 * f2, f3=f3))
    plt.legend()
    plt.show()

def task2_5d():
    # 2 различных по мощности несимметрично(свой вариант)
    x1 = a + 0.1
    x2 = b - 0.1
    f2 = DiracDelta(x - x1)
    f3 = DiracDelta(x - x2) 
    u = solve(1/20/k, a, b, ua, ub, f2 + 6 * f3)
    plot_solution(u, a, b, 'f2 < f3 = {f2} + {f3}'.format(f2=6 * f2, f3=f3))
    plt.legend()
    plt.show()


def task3():
    pass

def task4():
    pass


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

def plot_solution(foo, x1, x2, label):
    x = sp.symbols("x")
    # TODO: change 200 param
    args = np.linspace(x1, x2, 200)
    y = [foo.subs({x: t}) for t in args]
    plt.plot(args, y, label=label)
# TODO change output
def solve(k, a, b, ua, ub, f):
    x, c1, c2 = sp.symbols("x,c1,c2")
    u = -sp.integrate(sp.integrate(f, x) / k, x) + c1 * x + c2
    sol = sp.solve([u.subs({x: a}) - ua, u.subs({x: b}) - ub], (c1, c2))
    print(f'Input data:\nk={k}\na={a}\nb={b}\nua={ua}\nub={ub}\nf={f}')
    print('Resulting equation:', u)
    print('Solution:', sol)
    print('-' * 30)
    print()    
    return u.subs({c1: sol[c1], c2: sol[c2]})
# TODO change output
def small_solve(k, a, b, ua, ub, f):
    x, c1, c2 = sp.symbols('x c1 c2')
    u = -sp.integrate(sp.integrate(f, x) / k, x) + c1 * x + c2
    sol = sp.solve([u.subs({x: a}) - ua, u.subs({x: b}) - ub], (c1, c2))
    # print('Input data:\nk={}\na={}\nb={}\nua={}\nub={}\nf={}'.format(k, a, b, ua, ub, f))
    # print('Resulting equation:', u)
    print('Solution:', sol)
    # print('-' * 30)
    # print()    
    return u.subs({c1: sol[c1], c2: sol[c2]})

if __name__ == "__main__":
    main()