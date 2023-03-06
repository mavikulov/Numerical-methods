import math
import numpy as np


def f_1(X: np.array, u_1: float, u_2: float) -> float:
    return u_2


def f_2(X: np.array, u_1: float, u_2: float, a, b) -> float:
    return -a * u_2**2 - b * math.sin(u_1)


def get_coefficients(x, v, h, F: tuple, it, a, b) -> tuple:
    K1, K2, K3, K4 = [np.zeros(2) for _ in range(4)]
    if len(v.shape) == 2:
        K1[0] = F[0](x[it], v[0, it], v[1, it])
        K1[1] = F[1](x[it], v[0, it], v[1, it], a, b)
        K2[0] = F[0](x[it] + h[it] / 2, v[0, it] + h[it] / 2 * K1[0], v[1, it] + h[it] / 2 * K1[1])
        K2[1] = F[1](x[it] + h[it] / 2, v[0, it] + h[it] / 2 * K1[0], v[1, it] + h[it] / 2 * K1[1], a, b)
        K3[0] = F[0](x[it] + h[it] / 2, v[0, it] + h[it] / 2 * K2[0], v[1, it] + h[it] / 2 * K2[1])
        K3[1] = F[1](x[it] + h[it] / 2, v[0, it] + h[it] / 2 * K2[0], v[1, it] + h[it] / 2 * K2[1], a, b)
        K4[0] = F[0](x[it] + h[it], v[0, it] + h[it] * K3[0], v[1, it] + h[it] * K3[1])
        K4[1] = F[1](x[it] + h[it], v[0, it] + h[it] * K3[0], v[1, it] + h[it] * K3[1], a, b)
    else:
        K1[0] = F[0](x[it], v[0], v[1])
        K1[1] = F[1](x[it], v[0], v[1], a, b)
        K2[0] = F[0](x[it] + h[it] / 2, v[0] + h[it] / 2 * K1[0], v[1] + h[it] / 2 * K1[1])
        K2[1] = F[1](x[it] + h[it] / 2, v[0] + h[it] / 2 * K1[0], v[1] + h[it] / 2 * K1[1], a, b)
        K3[0] = F[0](x[it] + h[it] / 2, v[0] + h[it] / 2 * K2[0], v[1] + h[it] / 2 * K2[1])
        K3[1] = F[1](x[it] + h[it] / 2, v[0] + h[it] / 2 * K2[0], v[1] + h[it] / 2 * K2[1], a, b)
        K4[0] = F[0](x[it] + h[it], v[0] + h[it] * K3[0], v[1] + h[it] * K3[1])
        K4[1] = F[1](x[it] + h[it], v[0] + h[it] * K3[0], v[1] + h[it] * K3[1], a, b)
    return K1, K2, K3, K4


def RK45(x, v, h, F: tuple, it, a, b, eps, n_max, stop, control=True):
    count_of_divisions = 0
    doublings = [0]
    divisions = doublings[:]
    V2i = [v[:, 0].tolist(), ]
    OLP = [0]
    final_it = None
    while x[it] <= stop:
        if it > n_max:
            return x, v, V2i, h, OLP, doublings, divisions
        if x[it] + h[it] >= stop:
            final_it = it
            x = np.append(x, stop)
            h = np.append(h, stop - x[it])
            h[it] = stop - x[it]
            tmp = np.array([0, 0])
            v = np.append(v, tmp[:, np.newaxis], axis=1)

        k1, k2, k3, k4 = get_coefficients(x, v, h, F, it, a, b)
        Vi = v[:, it] + h[it] / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        current_h = np.copy(h)
        current_h[it] /= 2
        k1, k2, k3, k4 = get_coefficients(x, v, current_h, F, it, a, b)
        v_with_first_half_step = v[:, it] + current_h[it] / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        x_first_half = x + current_h
        k1, k2, k3, k4 = get_coefficients(x_first_half, v_with_first_half_step, current_h, F, it, a, b)
        v_with_second_half_step = v_with_first_half_step + current_h[it] / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        V2i.append(v_with_second_half_step.tolist())

        S = (v_with_second_half_step - Vi) / 15

        if control:
            if eps / 32 <= np.max(S) <= eps:
                if it != final_it:
                    x = np.append(x, x[it] + h[it])
                    v = np.append(v, Vi[:, np.newaxis], axis=1)
                    doublings.append(0)
                    divisions.append(count_of_divisions)
                    OLP.append(np.max(S) * 16)
                    count_of_divisions = 0
                    h = np.append(h, h[it])
                    it += 1
                else:
                    v[:, it + 1] = Vi
                    OLP.append(np.max(S) * 16)
                    return x, v, V2i, h, OLP, doublings, divisions

            elif np.max(S) < eps / 32:
                if it != final_it:
                    x = np.append(x, x[it] + h[it])
                    v = np.append(v, Vi[:, np.newaxis], axis=1)
                    doublings.append(1)
                    divisions.append(count_of_divisions)
                    OLP.append(np.max(S) * 16)
                    count_of_divisions = 0
                    h = np.append(h, h[it] * 2)
                    it += 1
                else:
                    v[:, it + 1] = Vi
                    OLP.append(np.max(S) * 16)
                    return x, v, V2i, h, OLP, doublings, divisions
            else:
                if it != final_it:
                    count_of_divisions += 1
                    h[it] /= 2
                else:
                    v[:, it + 1] = Vi
                    doublings.append(0)
                    divisions.append(0)
                    OLP.append(np.max(S) * 16)
                    return x, v, V2i, h, OLP, doublings, divisions

        else:
            if it != final_it:
                x = np.append(x, x[it] + h[it])
                v = np.append(v, Vi[:, np.newaxis], axis=1)
                OLP.append(np.max(S) * 16)
                doublings.append(0)
                divisions.append(count_of_divisions)
                h = np.append(h, h[it])
                it += 1
            else:
                v[:, it + 1] = Vi
                OLP.append(np.max(S) * 16)
                return x, v, V2i, h, OLP, doublings, divisions

    return x, v, V2i, h, OLP, doublings, divisions
