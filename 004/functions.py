import numpy as np
from sympy import *
from scipy.signal import savgol_filter

def erro_mach(kT_plot, x_plot, sep, mask_range):
    windows = [9, 11, 13]
    mach_ws = []
    for w in windows:
        kT_s = savgol_filter(kT_plot, window_length=w, polyorder=3)

        grad_kT = np.zeros_like(kT_s)
        grad_kT[mask_range] = np.gradient(kT_s[mask_range], x_plot[mask_range])

        idx = np.argmin(grad_kT)
        idx_menor = idx + sep
        idx_maior = np.argmax(kT_s)

        T1 = kT_s[idx_menor]
        T2 = kT_s[idx_maior]

        if T1 <= 0:
            continue

        M = symbols('M')
        eq = Eq((5*M**4 + 14*M**2 - 3)/(16*M**2), T2/T1)
        mach_solutions = solve(eq)
        mach_real = max([m for m in mach_solutions if im(m) == 0])

        mach_ws.append(float(mach_real))

    return 0.5 * (max(mach_ws) - min(mach_ws))

def fill_nan_nearest(arr):
    arr = np.array(arr, dtype=float)
    nans = np.isnan(arr)
    if not np.any(nans):
        return arr

    idx = np.arange(len(arr))
    valid_idx = idx[~nans]
    valid_values = arr[~nans]

    nearest_idx = np.abs(valid_idx[:, None] - idx).argmin(axis=0)
    arr[nans] = valid_values[nearest_idx[nans]]
    return arr

def velocidade_som_keV(kT_keV, mi=0.6, Mh=1.67262192e-27, gamma=5/3):
    """
    Calcula a velocidade do som a partir de kT em keV.
    """
    kT_joule = kT_keV * 1.60218e-16  # keV -> J
    cs = np.sqrt(gamma * kT_joule / (mi * Mh))  # m/s
    return cs / 1000  # em km/s


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
