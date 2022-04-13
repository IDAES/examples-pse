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
from pyomo.environ import units as pyunits, exp

from idaes.models.properties.modular_properties.base.generic_reaction import (
        ConcentrationForm)
from idaes.models.properties.modular_properties.reactions.dh_rxn import \
    constant_dh_rxn

from msr_equilibrium_constant import empirical_1, empirical_2
# custom expressions for the equilibrium constants as functions of temperature
# Source: Int. J. Hydrogen Energy, 42 (2017), pp. 2889-2903

from idaes.models.properties.modular_properties.reactions.equilibrium_forms import \
    power_law_equil

#import msr_pr as thermo_props

# For this example, the thermophysical properties are imported from the
# corresponding PR property package, msr_pr.py
#thermo_configuration = thermo_props.configuration

# Next, create the reaction property definition which describes the system on
# reactions to be modeled.

config_dict = {
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},
    "equilibrium_reactions": {
        #  Steam Methane Reformation Reaction
        "E1": {"stoichiometry": {("Vap", "CH4"): -1,
                                 ("Vap", "H2O"): -1,
                                 ("Vap", "CO"): 1,
                                 ("Vap", "H2"): 3,
                                 ("Vap", "CO2"): 0},
               "heat_of_reaction": constant_dh_rxn,
               "equilibrium_constant": empirical_1,
               "equilibrium_form": power_law_equil,
               "concentration_form": ConcentrationForm.moleFraction,
               "parameter_data": {  # Heat of reaction, est. Keq at Tref
                                    # Int. J. Hydrogen Energy, 42 (2017),
                                    # pp. 2889-2903
                    "dh_rxn_ref": (206000, pyunits.J/pyunits.mol),
                    "k_eq_ref": (9.49001757E-27, None),
                    "T_eq_ref": (298, pyunits.K)}},

        #  Water Gas Shift Reaction
        "E2": {"stoichiometry": {("Vap", "CH4"): 0,
                                 ("Vap", "H2O"): -1,
                                 ("Vap", "CO"): -1,
                                 ("Vap", "H2"): 1,
                                 ("Vap", "CO2"): 1},
               "heat_of_reaction": constant_dh_rxn,
               "equilibrium_constant": empirical_2,
               "equilibrium_form": power_law_equil,
               "concentration_form": ConcentrationForm.moleFraction,
               "parameter_data": {  # Heat of reaction, est. Keq at Tref
                                    # Int. J. Hydrogen Energy, 42 (2017),
                                    # pp. 2889-2903
                    "dh_rxn_ref": (-41000, pyunits.J/pyunits.mol),
                    "k_eq_ref": (45665.6052, None),
                    "T_eq_ref": (298, pyunits.K)}}}}
