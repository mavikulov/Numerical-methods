import numpy as np
import math
from functions import *
from numpy import linalg


def fill_bounds_v(v, x, y):
    for i in range(n + 1):
        v[i, 0] = u(x[i], 1)
        v[i, -1] = u(x[i], -1)
    for j in range(m + 1):
        v[0, j] = u(1, y[j])
        v[-1, j] = u(-1, y[j])


def fill_A():
    global n, m
    A = np.zeros(((n - 1) * (m - 1), (n - 1) * (m - 1)))
    np.fill_diagonal(A, coeff_A)

    for i, j in zip(range(0, A[0].size - (n - 1)), range(n - 1, np.size(A, 1))):
        A[i][j] = A[j][i] = coeff_k

    rows = [i for i in range((n - 1) * (m - 1))]
    columns = [i + 1 for i in range(0, (n - 1) * (m - 1) - 1)]
    for row, col in zip(rows, columns):
        if row % (n - 1) == n - 2:
            A[row][col] = 0
            A[col][row] = 0
        else:
            A[row][col] = coeff_h
            A[col][row] = coeff_h
    A = np.around(A, 2)
    return A


def generate_F(v):
    global coeff_k, coeff_A, coeff_h
    F = np.empty((n - 1) * (m - 1))
    k = 0
    for j in range(1, m):
        for i in range(1, n):
            if j == 1:
                if i == 1:
                    F[k] = -f(i, j) - coeff_h * v[i - 1, j] - coeff_k * v[i, j - 1]
                    k += 1
                elif i == n - 1:
                    F[k] = -f(i, j) - coeff_h * v[i + 1, j] - coeff_k * v[i, j - 1]
                    k += 1
                else:
                    F[k] = -f(i, j) - coeff_k * v[i, j - 1]
                    k += 1

            elif j == m - 1:
                if i == 1:
                    F[k] = -f(i, j) - coeff_h * v[i - 1, j] - coeff_k * v[i, j + 1]
                    k += 1
                elif i == n - 1:
                    F[k] = -f(i, j) - coeff_h * v[i + 1, j] - coeff_k * v[i, j + 1]
                    k += 1
                else:
                    F[k] = -f(i, j) - coeff_k * v[i, j + 1]
                    k += 1

            else:
                if i == 1:
                    F[k] = -f(i, j) - coeff_h * v[i - 1, j]
                    k += 1
                elif i == n - 1:
                    F[k] = -f(i, j) - coeff_h * v[i + 1, j]
                    k += 1
                else:
                    F[k] = -f(i, j)
                    k += 1

    return F


def seidel_method(v):
    global Nmax, eps, n, m, a, b, c, d, coeff_A, coeff_h, coeff_k
    fij = f()
    S = 0
    eps_max = 0
    eps_cur = 0
    flag = False
    h2 = -(n / (b - a)) ** 2
    k2 = -(m / (d - c)) ** 2
    a2 = -2 * (h2 + k2)

    while not flag:
        eps_max = 0
        for j in range(1, m):
            for i in range(1, n):
                v_old = v[i][j]
                v_new = -(h2 * (v[i + 1][j] + v[i - 1][j]) + k2 * (v[i][j + 1] + v[i][j - 1]))
                v_new = v_new + f(i, j)
                v_new = v_new / a2
                eps_cur = math.fabs(v_old - v_new)
                if eps_cur > eps_max:
                    eps_max = eps_cur

                v[i][j] = v_new
        S = S + 1

        if all([n < 20, m < 20]):
            if S == 1:
                print("Матрица V после 1-ой итерации метода Зейделя: \n")
                for row in v:
                    print(*row)
                print()
            if S == 2:
                print("Матрица V после 2-ой итерации метода Зейделя: \n")
                for row in v:
                    print(*row)
                print()

        if eps_max <= eps or S >= Nmax:
            flag = True

    return v, S, eps_max


print(
    "\t\tПрименение итерационного метода Зейделя для решения разностных схем\n\t\t    на примере задачи Дирихле для уравнения Пуассона")
a, b, c, d = -1, 1, -1, 1
print(f"{a} < x < {b},  {c} < y < {d}")
n, m = 4, 4
print(f"n = {n}, m = {m}")

h = (b - a) / n
k = (d - c) / m
x = [-1 + i * h for i in range(n + 1)]
y = [-1 + i * k for i in range(m + 1)]
coeff_h = 1 / h ** 2
coeff_k = 1 / k ** 2
coeff_A = -2 * (coeff_h + coeff_k)
print(f"Шаг h = {h}, шаг k = {k}")
print("Коэффициенты из разностной схемы: ")
print(f"1/h^2 = {coeff_h},  1/k^2 = {coeff_k},  A = {coeff_A}")
# print(f"Значения x узлов сетки", *x)
# print(f"Значения y узлов сетки", *y)
print()
print('Использовался итерационный метод Зейделя')
print("Тип начального приближения - нулевое")
print("Граничные условия и начальное приближение для решения СЛАУ")

v = np.zeros((n + 1, m + 1))
fill_bounds_v(v, x, y)
for row in np.around(v, 3):
    print(*row)
print('\n')

F = generate_F(v)
# print('\n\n Вектор F:', end=' ')
# print(F)
# print('\n\n')

print("Матрица A:")
A = fill_A()
for row in A:
    print(*row)

print('\n')
print('\t\tВвод параметров для критериев остановки метода Зейделя')
Nmax = int(input('Максимальное количество итераций метода Зейделя Nmax = '))
eps = float(input('Погрешность eps = '))
print('\n')

v, iterations, eps_max = seidel_method(v)

# print(f"\n\n, Матрица A * V =\n{A.dot(v[1:n, 1:m].flatten('F'))}")
# print('Вектор F: ', F)


print(f"При решении разностной схемы итерационным методом Зейделя с критериями остановки Nmax = {Nmax} и eps = {eps} "
      f"за {iterations} итераций достигнута точность {eps_max}")

print('Получено следующее численное решение:')
print('Матрица V, подсчитанная методом Зейделя:')
for row in v:
    print(*row)

print("\nМатрица U истинных значений:")

U = np.empty_like(v)
for i in range(n + 1):
    for j in range(m + 1):
        U[i, j] = u(x[i], y[j])
for row in np.around(U, 3):
     print(*row)

print("\nРазница между истинным и численным решением: ")
difference = np.abs(U - v)
for row in difference:
    print(*row)
print('\n\n')
print(f"Точность на выходе ε_N = {np.max(difference)}")

R_n = A.dot(v[1:n, 1:m].flatten('F')) - F
# print(f'R_n:')
# print(R_n)

R_n_2 = linalg.norm(R_n)

print(f"Невязка на выходе R_n = A * V_n - F имеет Евклидову норму = {R_n_2}")

inv_A = linalg.inv(A)
max_eigenval = max(map(abs, linalg.eigvals(inv_A)))

print(f"||Z_n||_oo <= ||Z_n||_2 <= {max_eigenval * R_n_2}")
print(f"|| A^(-1) ||_2 = {max_eigenval}")
print(f"||Z_n||_oo = {np.max(difference)}")
