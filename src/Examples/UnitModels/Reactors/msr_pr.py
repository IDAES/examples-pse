#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Property package for Methane Steam Reforming in Syngas to Hydrogen
using Peng-Robinson equation of state.

Author: Brandon Paul
"""
# Import Python libraries
import logging

# Import Pyomo units
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import VaporPhase, Component

from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.ceos import Cubic, CubicType
from idaes.generic_models.properties.core.pure.NIST import NIST

# Set up logger
_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------
# Configuration dictionary for a simplified natural gas system

# Property Sources

# Source: NIST webbook
# Properties: Heat capacity coefficients, reference enthalpies and entropies.

# Source: The Properties of Gases and Liquids (1987)
# 4th edition, Chemical Engineering Series - Robert C. Reid
# Properties: Critical temperatures and pressures, Omega.

config_dict = {
    # Specifying components
    "components": {
        'H2': {'type': Component,
               "elemental_composition": {"H": 2},
               'enth_mol_ig_comp': NIST,
               'entr_mol_ig_comp': NIST,
               'parameter_data': {
                   'mw': (0.0020159, pyunits.kg/pyunits.mol),
                   'pressure_crit': (13e5, pyunits.Pa),
                   'temperature_crit': (33.2, pyunits.K),
                   'omega': -0.218,
                   'cp_mol_ig_comp_coeff': {'A': 33.066178,
                                            'B': -11.363417,
                                            'C': 11.432816,
                                            'D': -2.772874,
                                            'E': -0.158558,
                                            'F': -9.980797,
                                            'G': 172.707974,
                                            'H': 0.0},
                   'enth_mol_form_vap_comp_ref': (
                       0, pyunits.J/pyunits.mol),
                   'entr_mol_form_vap_comp_ref': (
                       130.7, pyunits.J/pyunits.mol/pyunits.K)}},

        'CO': {'type': Component,
               "elemental_composition": {"C": 1, "O": 1},
               'enth_mol_ig_comp': NIST,
               'entr_mol_ig_comp': NIST,
               'parameter_data': {
                   'mw': (0.0280101, pyunits.kg/pyunits.mol),
                   'pressure_crit': (35e5, pyunits.Pa),
                   'temperature_crit': (132.9, pyunits.K),
                   'omega': 0.066,
                   'cp_mol_ig_comp_coeff': {'A': 25.56759,
                                            'B': 6.09613,
                                            'C': 4.054656,
                                            'D': -2.671301,
                                            'E': 0.131021,
                                            'F': -118.0089,
                                            'G': 227.3665,
                                            'H': -110.5271},
                   'enth_mol_form_vap_comp_ref': (
                       -110530, pyunits.J/pyunits.mol),
                   'entr_mol_form_vap_comp_ref': (
                       197.7, pyunits.J/pyunits.mol/pyunits.K)}},

        'H2O': {'type': Component,
                "elemental_composition": {"H": 2, "O": 1},
                'enth_mol_ig_comp': NIST,
                'entr_mol_ig_comp': NIST,
                'parameter_data': {
                    'mw': (0.01801528, pyunits.kg/pyunits.mol),
                    'pressure_crit': (221.2e5, pyunits.Pa),
                    'temperature_crit': (647.3, pyunits.K),
                    'omega': 0.344,
                    'cp_mol_ig_comp_coeff': {'A': 30.092,
                                             'B': 6.832514,
                                             'C': 6.793435,
                                             'D': -2.53448,
                                             'E': 0.082139,
                                             'F': -250.881,
                                             'G': 223.3967,
                                             'H': -241.8264},
                    'enth_mol_form_vap_comp_ref': (
                        -241826, pyunits.J/pyunits.mol),
                    'entr_mol_form_vap_comp_ref': (
                        188.8, pyunits.J/pyunits.mol/pyunits.K)}},

        'CO2': {'type': Component,
                "elemental_composition": {"C": 1, "O": 2},
                'enth_mol_ig_comp': NIST,
                'entr_mol_ig_comp': NIST,
                'parameter_data': {
                    'mw': (0.04401, pyunits.kg/pyunits.mol),
                    'pressure_crit': (73.8e5, pyunits.Pa),
                    'temperature_crit': (304.1, pyunits.K),
                    'omega': 0.239,
                    'cp_mol_ig_comp_coeff': {'A': 24.99735,
                                             'B': 55.18696,
                                             'C': -33.69137,
                                             'D': 7.948387,
                                             'E': -0.136638,
                                             'F': -403.6075,
                                             'G': 228.2431,
                                             'H': -393.5224},
                    'enth_mol_form_vap_comp_ref': (
                        -393510, pyunits.J/pyunits.mol),
                    'entr_mol_form_vap_comp_ref': (
                        213.8, pyunits.J/pyunits.mol/pyunits.K)}},

        'CH4': {'type': Component,
                "elemental_composition": {"C": 1, "H": 4},
                'enth_mol_ig_comp': NIST,
                'entr_mol_ig_comp': NIST,
                'parameter_data': {
                    'mw': (0.0160425, pyunits.kg/pyunits.mol),
                    'pressure_crit': (46e5, pyunits.Pa),
                    'temperature_crit': (190.4, pyunits.K),
                    'omega': 0.011,
                    'cp_mol_ig_comp_coeff': {'A': -0.703029,
                                             'B': 108.4773,
                                             'C': -42.52157,
                                             'D': 5.862788,
                                             'E': 0.678565,
                                             'F': -76.84376,
                                             'G': 158.7163,
                                             'H': -74.8731},
                    'enth_mol_form_vap_comp_ref': (
                        -74870, pyunits.J/pyunits.mol),
                    'entr_mol_form_vap_comp_ref': (
                        -186.3, pyunits.J/pyunits.mol/pyunits.K)}}},

    # Specifying phases
    "phases":  {'Vap': {"type": VaporPhase,
                        "equation_of_state": Cubic,
                        "equation_of_state_options":
                            {"type": CubicType.PR}}},

    # Set base units of measurement
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},

    # Specifying state definition
    "state_definition": FTPx,
    "state_bounds": {"flow_mol": (0, 8, 40, pyunits.kmol/pyunits.s),
                     "temperature": (273.15, 500, 1500, pyunits.K),
                     "pressure": (50, 101, 3500, pyunits.kPa)},
    "pressure_ref": (101.325, pyunits.kPa),
    "temperature_ref": (298.15, pyunits.K),

    "parameter_data": {"PR_kappa":
                       {('H2', 'H2'): 0, ('H2', 'CO'): 0,
                        ('H2', 'H2O'): 0, ('H2', 'CO2'): 0,
                        ('H2', 'CH4'): 0, ('CO', 'H2'): 0,
                        ('CO', 'CO'): 0, ('CO', 'H2O'): 0,
                        ('CO', 'CO2'): 0, ('CO', 'CH4'): 0,
                        ('H2O', 'H2'): 0, ('H2O', 'CO'): 0,
                        ('H2O', 'H2O'): 0, ('H2O', 'CO2'): 0,
                        ('H2O', 'CH4'): 0, ('CO2', 'H2'): 0,
                        ('CO2', 'CO'): 0, ('CO2', 'H2O'): 0,
                        ('CO2', 'CO2'): 0, ('CO2', 'CH4'): 0,
                        ('CH4', 'H2'): 0, ('CH4', 'CO'): 0,
                        ('CH4', 'H2O'): 0, ('CH4', 'CO2'): 0,
                        ('CH4', 'CH4'): 0, }}}
