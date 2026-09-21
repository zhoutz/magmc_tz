import numpy as np
from scipy.integrate import quad


def f(x, u, alpha):
    s2 = np.sin(alpha) ** 2
    c2 = 1 - s2
    x2 = x**2
    q = (2 - x2 - u * (1 - x2) ** 2 / (1 - u)) * s2
    tmp1 = np.sqrt(c2 + x2 * q)
    ret1 = x / tmp1
    return ret1


def cal_psi1(u, alpha):
    res, _err = quad(f, 0, 1, args=(u, alpha))
    s = np.sin(alpha)
    ret1 = 2 * s / np.sqrt(1 - u) * res
    return ret1


def cal_psi2(u, alpha):
    s = np.sin(alpha)
    tmp = (2 * s) / np.sqrt(3 * (1 - u))
    p_over_R = -tmp * np.cos((np.arccos(3 * u / tmp) + 2 * np.pi) / 3)
    psi_max = cal_psi1(u / p_over_R, np.pi / 2)
    psi1 = cal_psi1(u, np.pi - alpha)
    return 2 * psi_max - psi1


def cal_psi(u, alpha):
    if alpha <= np.pi / 2:
        return cal_psi1(u, alpha)
    else:
        return cal_psi2(u, alpha)


print(cal_psi(0.1, np.pi / 4))


def g(x, u, alpha):
    s2 = np.sin(alpha) ** 2
    c2 = 1 - s2
    x2 = x**2
    q = (2 - x2 - u * (1 - x2) ** 2 / (1 - u)) * s2
    tmp1 = np.sqrt(c2 + x2 * q)
    ret2 = x / (tmp1 * (1 + tmp1))
    return ret2


def cal_cdt_over_R1(u, alpha):
    res, _err = quad(g, 0, 1, args=(u, alpha))
    s = np.sin(alpha)
    ret2 = 2 * s**2 / (1 - u) * res
    return ret2


def cal_cdt_over_R2(u, alpha):
    s = np.sin(alpha)
    tmp = (2 * s) / np.sqrt(3 * (1 - u))
    p_over_R = -tmp * np.cos((np.arccos(3 * u / tmp) + 2 * np.pi) / 3)
    dt_max = cal_cdt_over_R1(u / p_over_R, np.pi / 2)
    dt1 = cal_cdt_over_R1(u, np.pi - alpha)
    return 2 * dt_max - dt1


def cal_cdt_over_R(u, alpha):
    if alpha <= np.pi / 2:
        return cal_cdt_over_R1(u, alpha)
    else:
        return cal_cdt_over_R2(u, alpha)
