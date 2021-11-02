import math
import numpy
from numpy import linspace
import pylab
import sympy
from sympy.solvers.solveset import linsolve
import functools

x = sympy.Symbol('x')

k = 11
num_var = 5

a = -1
b = 1

def task():
	dif_p = 5
	coeffs = [lambda x: 1 + math.cos(k) * x ** 2,
			  lambda x: 0,
			  lambda x: math.sin(k)]
	f = lambda x: -1
	print("")
	print( linspace(a + 0.2, b - 0.2, num_var))
	print("")
	print("Метод коллокаций: ", list(CollocationMethod.collocations_method(generate_basis_sequence(num_var),linspace(a + 0.2, b - 0.2, num_var))))
	print("Интегральный МНК: ", list(LSM.integral_least_square_method(generate_basis_sequence(num_var), a, b)))
	print("Дискретный МНК: ", list(LSM.discrete_least_square_method(generate_basis_sequence(num_var),num_var + dif_p, a, b)))
	print("Метод Галеркина: ", list(GalerkinMethod.galerkin_method(generate_basis_sequence(num_var), a, b)))
	collocation_points = numpy.linspace(a, b, 100)
	res_collocation = CollocationMethod.collocation_method(coeffs, f, collocation_points)
	res_integral = LSM.integral_LSM(coeffs, f, 50, a, b)
	points = numpy.linspace(a, b, 200)
	res_discrete = LSM.discrete_LSM(coeffs, f, 100, points)
	res_galerkin = GalerkinMethod.Galerkin_method(coeffs, f, 50, a, b)
	show_plots([[res_collocation]], a, b, 0.01, "Метод коллокаций")
	show_plots([[res_integral]], a, b, 0.01, "Интегральный МНК")
	show_plots([[res_discrete]], a, b, 0.01, "Дискретный МНК")
	show_plots([[res_galerkin]], a, b, 0.01, "Метод Галеркина")
	# show_plots([[res_collocation], [res_integral], [res_discrete],[res_galerkin]], a, b, 0.01, "All plots combined: ")

def func_for_substantiation(subs):
	func = sympy.sin(k) * sympy.diff(subs, x, x) + ((1 + sympy.cos(k) * x ** 2) * subs + 1)
	return func

def generate_basis_sequence(n):
	sequence = []
	for i in range(n):
		sequence.append((x ** i) * (1 - x ** 2))
	return sequence

def build_function_from_basis(basis):
	result = 0
	for i in range(len(basis)):
		current_a = sympy.Symbol('a' + str(i))
		result += current_a * basis[i]
	return result

def get_basis_function(n):
	if not n:
		return lambda x: 0
	return lambda x: x ** (n - 1) * (1 - x ** 2)


def get_basis_system(num_of_basis_functions):
	return [get_basis_function(n) for n in range(num_of_basis_functions)]

def numerical_integration(f, a, b):
	n = 100
	dx = (b - a) / n
	xlist = numpy.arange(a, b, dx)
	ylist = [f(p) for p in xlist]
	return numpy.trapz(ylist, dx=dx)


def numerical_diff(f, x, n):
	h = 0.001
	if not n:
		return f(x)
	elif n == 1:
		return (f(x + h) - f(x - h)) / (2 * h)
	elif n == 2:
		return (f(x - h) - 2 * f(x) + f(x + h)) / (h ** 2)
	else:
		raise NotImplementedError

def show_plots(functions, start_x, end_x, dx, title):
	for f in functions:
		for function in f:
			x_list = numpy.arange(start_x, end_x, dx)
			y_list = [function(p) for p in x_list]
			pylab.plot(x_list, y_list)
			pylab.title(title)
	pylab.grid(True)
	pylab.show()
	

class CollocationMethod:
	def collocations_method(basis, points):
		func = build_function_from_basis(basis)
		psi_func = func_for_substantiation(func)
		symbols = [sympy.Symbol('a' + str(i)) for i in range(len(points))]
		lin_system = []
		for point in points:
			lin_system.append(psi_func.subs(x, point).evalf())
		answer = linsolve(lin_system, *symbols)
		return answer

	def collocation_method(coefficients, f, collocation_points):
		num_of_basis_functions = len(collocation_points) + 1
		basis = get_basis_system(num_of_basis_functions)

		def resudial_part_diff(a):
			return lambda x: sum(coefficients[i](x) * numerical_diff(basis[a], x, i)
								 for i in range(len(coefficients)))

		matrix = [[resudial_part_diff(i)(point) for i in range(1, num_of_basis_functions)]
				  for point in collocation_points]
		right_side = [f(point) - resudial_part_diff(0)(point)
					  for point in collocation_points]
		answer = numpy.linalg.solve(numpy.matrix(matrix), numpy.array(right_side))
		return lambda x: basis[0](x) + sum(answer[i - 1] * basis[i](x) for i in range(1, num_of_basis_functions))

class GalerkinMethod:
	def Galerkin_method(coefficients, f, num_of_basis_functions, a, b):
		basis = get_basis_system(num_of_basis_functions)

		def resudial_part_diff(a):
			return lambda x: sum(coefficients[i](x) * numerical_diff(basis[a], x, i)
								 for i in range(len(coefficients)))

		part_diffs = [resudial_part_diff(i) for i in range(num_of_basis_functions)]

		def optimal_numerical_integration(f1, f2, points):
			dx = 2 / len(points);
			return sum([(f1[i] * f2[i] * dx) for i in range(0, len(points))])

		points = numpy.linspace(a, b, 100)
		points = [(points[i] + points[i - 1]) / 2 for i in range(1, len(points))]

		opt_part_diffs = [None] * num_of_basis_functions
		opt_basis = [None] * num_of_basis_functions
		for i in range(1, num_of_basis_functions):
			opt_part_diffs[i] = []
			opt_basis[i] = []
			for point in points:
				opt_part_diffs[i].append(part_diffs[i](point))
				opt_basis[i].append(basis[i](point))

		matrix = [[optimal_numerical_integration(opt_part_diffs[j], opt_basis[i], points)
				   for j in range(1, num_of_basis_functions)]
				  for i in range(1, num_of_basis_functions)]
		right_side = [numerical_integration(
			lambda x: (f(x) - part_diffs[0](x)) * basis[i](x), a, b)
			for i in range(1, num_of_basis_functions)]
		answer = numpy.linalg.solve(numpy.matrix(matrix), numpy.array(right_side))
		return lambda x: basis[0](x) + sum(answer[i - 1] * basis[i](x)
											   for i in range(1, num_of_basis_functions))
	def galerkin_method(basis, a, b):
		func = build_function_from_basis(basis)
		psi_func = func_for_substantiation(func)
		symbols = [sympy.Symbol('a' + str(i)) for i in range(len(basis))]
		lin_system = []
		for i in range(len(basis)):
			lin_system.append(sympy.integrate(
				psi_func * basis[i], (x, a, b)).evalf())
		answer = linsolve(lin_system, *symbols)
		return answer

class LSM:
	def integral_LSM(coefficients, f, num_of_basis_functions, a, b):
		basis = get_basis_system(num_of_basis_functions)

		def resudial_part_diff(a):
			return lambda x: sum(coefficients[i](x) * numerical_diff(basis[a], x, i)
								 for i in range(len(coefficients)))

		part_diffs = [resudial_part_diff(i) for i in range(num_of_basis_functions)]

		matrix = [[numerical_integration(lambda x: part_diffs[i](x) * part_diffs[j](x), a, b)
				   for j in range(1, num_of_basis_functions)]
				  for i in range(1, num_of_basis_functions)]
		right_side = [numerical_integration(
			lambda x: (f(x) - part_diffs[0](x)) * part_diffs[i](x), a, b)
			for i in range(1, num_of_basis_functions)]
		answer = numpy.linalg.solve(numpy.matrix(matrix), numpy.array(right_side))
		return lambda x: basis[0](x) + sum(answer[i - 1] * basis[i](x)
										   for i in range(1, num_of_basis_functions))

	def discrete_LSM(coefficients, f, num_of_basis_functions, points):
		basis = get_basis_system(num_of_basis_functions)

		def resudial_part_diff(a, x):
			return sum(coefficients[i](x) * numerical_diff(basis[a], x, i)
					   for i in range(len(coefficients)))

		part_diffs_x = {(i, x): resudial_part_diff(i, x) for i in range(num_of_basis_functions)
						for x in points}

		matrix = [[sum(part_diffs_x[(i, x)] * part_diffs_x[(j, x)] for x in points)
				   for j in range(1, num_of_basis_functions)]
				  for i in range(1, num_of_basis_functions)]
		right_side = [sum([(f(x) - part_diffs_x[(0, x)]) * part_diffs_x[(i, x)] for x in points])
					  for i in range(1, num_of_basis_functions)]
		answer = numpy.linalg.solve(numpy.matrix(matrix), numpy.array(right_side))
		return lambda x: basis[0](x) + sum(answer[i - 1] * basis[i](x)
										   for i in range(1, num_of_basis_functions))

	def integral_least_square_method(basis, a, b):
		func = build_function_from_basis(basis)
		psi_func = func_for_substantiation(func)
		symbols = [sympy.Symbol('a' + str(i)) for i in range(len(basis))]
		lin_system = []
		for i in range(len(basis)):
			lin_system.append(sympy.integrate(2 * sympy.diff(
				psi_func, symbols[i]) * psi_func, (x, a, b)).evalf())
		answer = linsolve(lin_system, *symbols)
		return answer


	def discrete_least_square_method(basis, points_num, a, b):
		func = build_function_from_basis(basis)
		psi_func = func_for_substantiation(func)
		seq = [psi_func.subs(x, point) ** 2 for point
			   in linspace(a + 0.05, b - 0.05, points_num)]
		psi_sqr_sum = functools.reduce((lambda a, b: a + b), seq)
		symbols = [sympy.Symbol('a' + str(i)) for i in range(len(basis))]
		lin_system = []
		for i in range(len(basis)):
			lin_system.append(sympy.diff(psi_sqr_sum, symbols[i]).evalf())
		answer = linsolve(lin_system, *symbols)
		return answer

if __name__=='__main__':
	task()
