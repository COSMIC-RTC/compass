#!/usr/local/bin/python2.7
# encoding: utf-8
"""
Created on 5 juil. 2017

@author: vdeo

Numerical constants for shesha
Config enumerations for safe-typing

"""

import numpy as np


class CONST:
    RAD2ARCSEC = 3600. * 360. / (2 * np.pi)
    ARCSEC2RAD = 2. * np.pi / (360. * 3600.)
    RAD2DEG = 180. / np.pi
    DEG2RAD = np.pi / 180.


def check_enum(cls, name):
    """
        Create a safe-type enum instance from bytes contents
    """

    if isinstance(name, str):
        name = bytes(name.encode('UTF-8'))

    if not isinstance(name, bytes) or \
            not name in vars(cls).values():
        raise ValueError(
                "Invalid enumeration value for enum %s, value %s" %
                (cls, name))
    return name


class DmType:
    """
        Types of deformable mirrors
    """

    PZT = b'pzt'
    TT = b'tt'
    KL = b'kl'


class PatternType:
    """
        Types of Piezo DM patterns
    """

    SQUARE = b'square'
    HEXA = b'hexa'
    HEXAM4 = b'hexaM4'


class KLType:
    """
        Possible KLs for computations
    """

    KOLMO = b'kolmo'
    KARMAN = b'karman'


class InfluType:
    """
        Influence function types
    """

    DEFAULT = b'default'
    RADIALSCHWARTZ = b'radialSchwartz'
    SQUARESCHWARTZ = b'squareSchwartz'
    BLACKNUTT = b'blacknutt'
    GAUSSIAN = b'gaussian'
    BESSEL = b'bessel'


class ControllerType:
    """
        Controller types
    """

    GENERIC = b'generic'
    LS = b'ls'
    MV = b'mv'
    CURED = b'cured'
    GEO = b'geo'
    KALMAN_C = b'kalman_CPU'
    KALMAN_G = b'kalman_GPU'
    KALMAN_UN = b'kalman_uninitialized'


class CentroiderType:
    """
        Centroider types
    """

    COG = b'cog'
    TCOG = b'tcog'
    WCOG = b'wcog'
    BPCOG = b'bpcog'
    CORR = b'corr'
    PYR = b'pyr'


class CentroiderFctType:
    MODEL = b'model'
    GAUSS = b'gauss'


class WFSType:
    """
        WFS Types
    """
    SH = b'sh'
    PYRHR = b'pyrhr'


class TargetImageType:
    """
        Target Images
    """

    SE = b'se'
    LE = b'le'


class ApertureType:
    """
        Telescope apertures
    """
    GENERIC = b'Generic'
    EELT_NOMINAL = b'EELT-Nominal'
    EELT_BP1 = b'EELT-BP1'
    EELT_BP3 = b'EELT-BP3'
    EELT_BP5 = b'EELT-BP5'
    EELT_CUSTOM = b'EELT-Custom'
    VLT = b'VLT'


class SpiderType:
    """
        Spiders
    """
    FOUR = b'four'
    SIX = b'six'


class ProfType:
    """
        Sodium profile for LGS
    """
    GAUSS1 = b'Gauss1'
    GAUSS2 = b'Gauss2'
    GAUSS3 = b'Gauss3'
    EXP = b'Exp'
    FILES = dict({
            GAUSS1: b"allProfileNa_withAltitude_1Gaussian.npy",
            GAUSS2: b"allProfileNa_withAltitude_2Gaussian.npy",
            GAUSS3: b"allProfileNa_withAltitude_3Gaussian.npy",
            EXP: b"allProfileNa_withAltitude.npy"})


class FieldStopType:
    """
        WFS field stop
    """
    SQUARE = b'square'
    ROUND = b'round'