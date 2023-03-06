import math


def f(x, u, a1, a2, m):
    return -(a1 * u + a2 * u**2) / m


def get_coefficients(x, v, h, a1, a2, m):
    K1, K2, K3 = [None for _ in range(3)]
    K1 = f(x, v, a1, a2, m)
    K2 = f(x + h / 2, v + h / 2 * K1, a1, a2, m)
    K3 = f(x + h, v + h * (-K1 + 2 * K2), a1, a2, m)
    return K1, K2, K3


def RK34_without_control(x, v, h, it, n_max, stop, a1, a2, m, u0, to_bound=None):
    u = [u0, ]
    while x[it] <= stop - to_bound:
        if abs(v[it] - u[it]) > 100:
            return x, v, h, u
        if it >= n_max:
            return x, v, h, u

        if abs(v[it]) < 1e-17 and it > 0:
            return x, v, h, u

        k1, k2, k3 = get_coefficients(x[it], v[it], h[it], a1, a2, m)
        if x[it] + h[it] > stop - to_bound:
            x[-1] = stop - to_bound
            h[-2] = h[0] - to_bound
            return x, v, h, u
        v.append(v[it] + h[it] / 6 * (k1 + 4 * k2 + k3))
        u.append(-a1 / (a2 - (math.exp(a1 * (x[it] + h[it]) / m) * (a2 * u0 + a1) / u0)))
        x.append(x[it] + h[it])
        h.append(h[it])
        it += 1
    return x, v, h, u


def RK34_with_control(x, v, h, it, eps, n_max, stop, a1, a2, m, u0, to_bound=None):
    u = [u0, ]
    count_of_divisions = 0
    OLP = [0]
    doublings = [0]
    divisions = doublings[:]
    V2i = v[:]
    while x[it] < stop - to_bound:
        if it >= n_max:
            return x, v, h, V2i, OLP, doublings, divisions, u

        if abs(v[it]) < 1e-17 and it > 0:
            return x, v, h, V2i, OLP, doublings, divisions, u

        if abs(v[it] - u[it]) > 100:
            return x, v, h, V2i, OLP, doublings, divisions, u

        if x[it] + h[it] >= stop - to_bound:
            h[it] = stop - x[it] - to_bound

        k1, k2, k3 = get_coefficients(x[it], v[it], h[it], a1, a2, m)
        Vi = v[it] + h[it] / 6 * (k1 + 4 * k2 + k3)
        half_step = h[it] / 2
        x_half = x[it] + half_step
        k1, k2, k3 = get_coefficients(x[it], v[it], half_step, a1, a2, m)
        v_with_first_half_step = v[it] + half_step / 6 * (k1 + 4 * k2 + k3)
        k1, k2, k3 = get_coefficients(x_half, v_with_first_half_step, half_step, a1, a2, m)
        v_with_double_step = (v_with_first_half_step + half_step / 6 * (k1 + 4 * k2 + k3))

        S = (v_with_double_step - Vi) / 7

        if eps / 16 <= math.fabs(S) <= eps:
            x.append(x[it] + h[it])
            v.append(Vi)
            u.append(-a1 / (a2 - (math.exp(a1 * (x[it] + h[it]) / m) * (a2 * u0 + a1) / u0)))
            V2i.append(v_with_double_step)
            h.append(h[it])
            doublings.append(0)
            OLP.append(math.fabs(S) * 8)
            divisions.append(count_of_divisions)
            count_of_divisions = 0
            it += 1
        elif math.fabs(S) < eps / 16:
            x.append(x[it] + h[it])
            v.append(Vi)
            u.append(-a1 / (a2 - (math.exp(a1 * (x[it] + h[it]) / m) * (a2 * u0 + a1) / u0)))
            V2i.append(v_with_double_step)
            h.append(h[it] * 2)
            doublings.append(1)
            divisions.append(count_of_divisions)
            count_of_divisions = 0
            OLP.append(math.fabs(S) * 8)
            it += 1
        else:
            count_of_divisions += 1
            h[it] /= 2
    return x, v, h, V2i, OLP, doublings, divisions, u
