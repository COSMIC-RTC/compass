#
# This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
#
# COMPASS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# COMPASS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with COMPASS. If not, see <https://www.gnu.org/licenses/>.
#
# Copyright (C) 2011-2025 COSMIC Team
"""
STARLORD (SeT of Algorithms foR mOdified stRucture function computation)
Set of functions for structure function computation
"""

import numpy as np
from scipy.special import jv  # Bessel function


def dphi_highpass(r, x0, tabx, taby):
    """High-pass component of the phase structure function (von Kármán model).

    Computes the high-spatial-frequency contribution to the phase structure
    function (spatial frequencies above 1/(2*x0)). Result must be scaled by
    (1/r0)**(5/3) to obtain physical units.

    Args:
        r (np.ndarray): Separation distances [m].
        x0 (float): Inter-actuator pitch [m]; sets the high-pass cut-off frequency.
        tabx (np.ndarray): Abscissa of the tabulated Ij0 integral (from tabulateIj0).
        taby (np.ndarray): Values of the tabulated Ij0 integral (from tabulateIj0).

    Returns:
        np.ndarray: High-pass phase structure function, dimensionless (scale by (1/r0)**(5/3)).
    """
    return (
        (r ** (5.0 / 3.0))
        * (1.1183343328701949 - Ij0t83(r * (np.pi / x0), tabx, taby))
        * (2 * (2 * np.pi) ** (8 / 3.0) * 0.0228956)
    )


def dphi_lowpass(r, x0, L0, tabx, taby):
    """Low-pass component of the phase structure function (von Kármán model).

    Computes the low-spatial-frequency contribution to the phase structure
    function (spatial frequencies below 1/(2*x0)), including outer-scale effects.
    Result must be scaled by (1/r0)**(5/3) to obtain physical units.

    Args:
        r (np.ndarray): Separation distances [m].
        x0 (float): Inter-actuator pitch [m]; sets the low-pass cut-off frequency.
        L0 (float): Outer scale of atmospheric turbulence [m].
        tabx (np.ndarray): Abscissa of the tabulated Ij0 integral (from tabulateIj0).
        taby (np.ndarray): Values of the tabulated Ij0 integral (from tabulateIj0).

    Returns:
        np.ndarray: Low-pass phase structure function, dimensionless (scale by (1/r0)**(5/3)).
    """
    return rodconan(r, L0) - dphi_highpass(r, x0, tabx, taby)


def Ij0t83(x, tabx, taby):
    """Evaluate the tabulated integral of t^(-8/3) * (1 - J0(t)) from 0 to x.

    The integral is:

        I(x) = integral_0^x  t^(-8/3) * (1 - J0(t)) dt

    where J0 is the Bessel function of the first kind of order 0.
    Near the origin the approximation I(x) ≈ (3/4) x^(1/3) (1 - x²/112 + …)
    is used; for larger x the result is looked up from (tabx, taby).

    Args:
        x (np.ndarray): Evaluation points (must be ≥ 0).
        tabx (np.ndarray): Pre-computed abscissa grid (from tabulateIj0).
        taby (np.ndarray): Pre-computed integral values (from tabulateIj0).

    Returns:
        np.ndarray: Integral values, same shape as *x*.
    """
    res = x.copy()
    ismall = np.where(res < np.exp(-3.0))
    ilarge = np.where(res >= np.exp(-3.0))
    if ismall[0].size > 0:
        res[ismall] = 0.75 * x[ismall] ** (1.0 / 3) * (1 - x[ismall] ** 2 / 112.0)
    if ilarge[0].size > 0:
        res[ilarge] = np.interp(x[ilarge], tabx, taby)

    return res


def tabulateIj0():
    """Pre-compute the lookup table for the Ij0t83 integral.

    Tabulates the integral I(x) = integral_0^x t^(-8/3) (1 - J0(t)) dt over a
    logarithmic grid covering [exp(-4), exp(10)].  The result must be passed as
    (tabx, taby) to dphi_lowpass and dphi_highpass.

    Returns:
        tuple[np.ndarray, np.ndarray]: (tabx, taby) where *tabx* is the abscissa
        grid and *taby* the corresponding integral values.
    """
    n = 10000
    t = np.linspace(-4, 10, n)
    dt = (t[-1] - t[0]) / (n - 1)
    smallx = np.exp(-4.0)
    A = 0.75 * smallx ** (1.0 / 3) * (1 - smallx**2 / 112.0)
    X = np.exp(t)
    Y = np.exp(-t * (5.0 / 3.0)) * (1 - jv(0, X))
    Y[1:] = np.cumsum(Y[:-1] + np.diff(Y) / 2.0)
    Y[0] = 0.0
    Y = Y * dt + A

    return X, Y


def asymp_macdo(x):
    """Asymptotic expansion of the Macdo function for large arguments (x > 4.71).

    Used internally by rodconan to evaluate the von Kármán phase structure
    function.  The expansion is:

        asymp_macdo(x) ≈ k2 - k3 * exp(-x) * x^(1/3) * (1 + a1/x + a2/x² + a3/x³)

    where the coefficients reproduce the large-argument behaviour of the
    integral that defines the generalised structure function.

    Args:
        x (np.ndarray): Dimensionless separation, x = 2π r / L0.  Must satisfy x > 4.71.

    Returns:
        np.ndarray: Asymptotic approximation, same shape as *x*.
    """
    k2 = 1.00563491799858928388289314170833
    k3 = 1.25331413731550012081
    a1 = 0.22222222222222222222
    a2 = -0.08641975308641974829
    a3 = 0.08001828989483310284

    x_1 = 1.0 / x
    res = k2 - k3 * np.exp(-x) * x ** (1.0 / 3.0) * (1.0 + x_1 * (a1 + x_1 * (a2 + x_1 * a3)))
    return res


def macdo(x):
    """Evaluate the Macdo function using a power-series expansion for small arguments (x ≤ 4.71).

    The Macdo function arises in the von Kármán phase structure function and is
    defined through an integral involving modified Bessel functions.  For small
    arguments the series

        macdo(x) = sum_n [ Gma[n] * x^(5/3) + Ga[n] ] * (x²/4)^n

    converges rapidly.

    Args:
        x (float | np.ndarray): Dimensionless separation, x = 2π r / L0.  Must satisfy x ≤ 4.71.

    Returns:
        float | np.ndarray: Function value, same type/shape as *x*.
    """
    a = 5.0 / 6.0
    x2a = x ** (2.0 * a)
    x22 = x * x / 4.0
    s = 0.0

    Ga = [
        0,
        12.067619015983075,
        5.17183672113560444,
        0.795667187867016068,
        0.0628158306210802181,
        0.00301515986981185091,
        9.72632216068338833e-05,
        2.25320204494595251e-06,
        3.93000356676612095e-08,
        5.34694362825451923e-10,
        5.83302941264329804e-12,
    ]

    Gma = [
        -3.74878707653729304,
        -2.04479295083852408,
        -0.360845814853857083,
        -0.0313778969438136685,
        -0.001622994669507603,
        -5.56455315259749673e-05,
        -1.35720808599938951e-06,
        -2.47515152461894642e-08,
        -3.50257291219662472e-10,
        -3.95770950530691961e-12,
        -3.65327031259100284e-14,
    ]

    x2n = 0.5

    s = Gma[0] * x2a
    s *= x2n

    x2n *= x22

    for n in np.arange(10) + 1:
        s += (Gma[n] * x2a + Ga[n]) * x2n
        x2n *= x22

    return s


def rodconan(r, L0):
    """Von Kármán phase structure function with outer-scale correction.

    Evaluates the normalised (r0-independent) phase structure function for
    Kolmogorov turbulence modified by a finite outer scale L0:

        D_phi(r) = (r0)^(-5/3) * rodconan(r, L0)

    The piecewise implementation uses asymp_macdo for 2πr/L0 > 4.71 and
    macdo otherwise to ensure numerical accuracy across all separations.

    Args:
        r (np.ndarray): Separation distances [m].
        L0 (float): Outer scale of atmospheric turbulence [m].

    Returns:
        np.ndarray: Normalised structure function values, same shape as *r*.
        Multiply by (1/r0)**(5/3) to get physical phase variance [rad²].
    """
    res = r * 0.0
    k1 = 0.1716613621245709486
    dprf0 = (2 * np.pi / L0) * r
    ilarge = np.where(dprf0 > 4.71239)
    ismall = np.where(dprf0 <= 4.71239)
    if ilarge[0].size > 0:
        res[ilarge] = asymp_macdo(dprf0[ilarge])
    if ismall[0].size > 0:
        res[ismall] = -macdo(dprf0[ismall])

    res *= k1 * L0 ** (5.0 / 3.0)

    return res
