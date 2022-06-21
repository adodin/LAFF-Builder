"""
Implementation of Different Mixing Rules
Written By: Amro Dodin (Geissler Group UC Berkeley)
"""

import numpy as np

fpe0 =  0.000719756                     # (4 Pi eps0) in e^2/(kJ/mol A)


def arithmetic(x, y):
    return (x+y)/2


def geometric(x, y):
    return np.sqrt(x*y)


def detect_drude(mratio):
    assert 0 <= mratio <= 1
    if mratio == 1:
        return 'NP'
    elif mratio > 0.5:
        return'C'
    else:
        return 'D'


def drude_symmetrized(x, y, mratio_x, mratio_y, base_mixing_rule=geometric):
    # Compute Uncorrected Coefficients from Base Mixing Rule
    base = base_mixing_rule(x, y)
    # Autodetect if x or y are Drude Core, Shell or NonPolarizable (using mass ratio)
    xtype = detect_drude(mratio_x)
    ytype = detect_drude(mratio_y)
    # Calculate Correction Factor
    # Any particle is non-polarizable
    if 'NP' in (xtype, ytype):
        corr = mratio_x*mratio_y
    elif xtype == 'C' and ytype == 'C':
        corr = mratio_x + mratio_y - 1
    elif xtype == 'C' and ytype == 'D':
        corr = mratio_y
    elif xtype == 'D' and ytype == 'C':
        corr = mratio_x
    else:
        corr = 0
    return corr * base


def drude_asymmetric(x, y, mratio_x, mratio_y, base_mixing_rule=geometric):
    # Compute Uncorrected Coefficients from Base Mixing Rule
    base = base_mixing_rule(x, y)
    # Autodetect if x or y are Drude Core, Shell or NonPolarizable (using mass ratio)
    xtype = detect_drude(mratio_x)
    ytype = detect_drude(mratio_y)
    # Calculate Correction Factor
    # Any particle is non-polarizable
    if 'D' in (xtype, ytype):
        corr = 0
    else:
        corr = 1
    return corr * base


def roux_convention(alpha, drude_mass=0.8):
    k = 4184.0                              # Spring Constant in kJ/molA^2 (by Convention)
    return np.sqrt(fpe0*k*alpha)
