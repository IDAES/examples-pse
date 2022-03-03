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
Methods for calculating equilibrium constants
"""
from pyomo.environ import exp, log, Var, units as pyunits, value

from idaes.generic_models.properties.core.generic.utility import \
    ConcentrationForm
from idaes.generic_models.properties.core.reactions.dh_rxn import \
    constant_dh_rxn
from idaes.core.util.misc import set_param_from_config
from idaes.core.util.constants import Constants as c
from idaes.core.util.exceptions import BurntToast, ConfigurationError


# -----------------------------------------------------------------------------
# custom class from empirical correlations, modified from van_t_hoff (line 97)
# Source: Int. J. Hydrogen Energy, 42 (2017), pp. 2889-2903
class empirical_1():  # steam methane reforming, reaction 1

    @staticmethod
    def build_parameters(rblock, config):
        parent = rblock.parent_block()
        units = parent.get_metadata().derived_units

        c_form = config.concentration_form
        if c_form is None:
            raise ConfigurationError(
                "{} concentration_form configuration argument was not set. "
                "Please ensure that this argument is included in your "
                "configuration dict.".format(rblock.name))
        elif (c_form == ConcentrationForm.moleFraction or
              c_form == ConcentrationForm.massFraction or
              c_form == ConcentrationForm.activity):
            e_units = pyunits.dimensionless
        else:
            order = 0

            try:
                # This will work for Reaction Packages
                pc_set = parent.config.property_package._phase_component_set
            except AttributeError:
                # Need to allow for inherent reactions in Property Packages
                if not parent._electrolyte:
                    # In most cases ,should have _phase_component_set
                    pc_set = parent._phase_component_set
                else:
                    # However, for electrolytes need true species set
                    pc_set = parent.true_phase_component_set

            for p, j in pc_set:
                order += rblock.reaction_order[p, j].value

            if c_form == ConcentrationForm.molarity:
                c_units = units["density_mole"]
            elif c_form == ConcentrationForm.molality:
                c_units = units["amount"]*units["mass"]**-1
            elif c_form == ConcentrationForm.partialPressure:
                c_units = units["pressure"]
            else:
                raise BurntToast(
                    "{} received unrecognised ConcentrationForm ({}). "
                    "This should not happen - please contact the IDAES "
                    "developers with this bug."
                    .format(rblock.name, c_form))

            e_units = c_units**order

        rblock.k_eq_ref = Var(
                doc="Equilibrium constant at reference state",
                units=e_units)
        set_param_from_config(rblock, param="k_eq_ref", config=config)

        rblock.T_eq_ref = Var(
                doc="Reference temperature for equilibrium constant",
                units=units["temperature"])
        set_param_from_config(rblock, param="T_eq_ref", config=config)

    @staticmethod
    def return_expression(b, rblock, r_idx, T):
        units = pyunits.get_units(rblock.k_eq_ref)
        if units is None or units is pyunits.dimensionless:
            return exp(b.log_k_eq[r_idx])
        else:
            return exp(b.log_k_eq[r_idx]) * units

    @staticmethod
    def return_log_expression(b, rblock, r_idx, T):
        return (b.log_k_eq[r_idx] == -26830/T + 30.114)

    @staticmethod
    def calculate_scaling_factors(b, rblock):
        return 1/value(rblock.k_eq_ref)


# -----------------------------------------------------------------------------
# custom class from empirical correlations, modified from van_t_hoff (line 178)
# Source: Int. J. Hydrogen Energy, 42 (2017), pp. 2889-2903
class empirical_2():  # water gas shift, reaction 2

    @staticmethod
    def build_parameters(rblock, config):
        parent = rblock.parent_block()
        units = parent.get_metadata().derived_units

        c_form = config.concentration_form
        if c_form is None:
            raise ConfigurationError(
                "{} concentration_form configuration argument was not set. "
                "Please ensure that this argument is included in your "
                "configuration dict.".format(rblock.name))
        elif (c_form == ConcentrationForm.moleFraction or
              c_form == ConcentrationForm.massFraction or
              c_form == ConcentrationForm.activity):
            e_units = pyunits.dimensionless
        else:
            order = 0

            try:
                # This will work for Reaction Packages
                pc_set = parent.config.property_package._phase_component_set
            except AttributeError:
                # Need to allow for inherent reactions in Property Packages
                if not parent._electrolyte:
                    # In most cases ,should have _phase_component_set
                    pc_set = parent._phase_component_set
                else:
                    # However, for electrolytes need true species set
                    pc_set = parent.true_phase_component_set

            for p, j in pc_set:
                order += rblock.reaction_order[p, j].value

            if c_form == ConcentrationForm.molarity:
                c_units = units["density_mole"]
            elif c_form == ConcentrationForm.molality:
                c_units = units["amount"]*units["mass"]**-1
            elif c_form == ConcentrationForm.partialPressure:
                c_units = units["pressure"]
            else:
                raise BurntToast(
                    "{} received unrecognised ConcentrationForm ({}). "
                    "This should not happen - please contact the IDAES "
                    "developers with this bug."
                    .format(rblock.name, c_form))

            e_units = c_units**order

        rblock.k_eq_ref = Var(
                doc="Equilibrium constant at reference state",
                units=e_units)
        set_param_from_config(rblock, param="k_eq_ref", config=config)

        rblock.T_eq_ref = Var(
                doc="Reference temperature for equilibrium constant",
                units=units["temperature"])
        set_param_from_config(rblock, param="T_eq_ref", config=config)

    @staticmethod
    def return_expression(b, rblock, r_idx, T):
        units = pyunits.get_units(rblock.k_eq_ref)
        if units is None or units is pyunits.dimensionless:
            return exp(b.log_k_eq[r_idx])
        else:
            return exp(b.log_k_eq[r_idx]) * units

    @staticmethod
    def return_log_expression(b, rblock, r_idx, T):
        return (b.log_k_eq[r_idx] == 4400/T - 4.036)

    @staticmethod
    def calculate_scaling_factors(b, rblock):
        return 1/value(rblock.k_eq_ref)
