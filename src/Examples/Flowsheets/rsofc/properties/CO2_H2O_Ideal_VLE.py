##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Property package for MEA Carbon Capture System Condenser

CO2-H2O VLE system
Ideal VLE with non-condensables
"""
# Import Python libraries
import logging
import copy
import enum

from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import LiquidPhase, VaporPhase, Component, PhaseType

from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.forms import fugacity

from idaes.models.properties.modular_properties.pure import Perrys, NIST

# Set up logger
_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------
# Configuration dictionary for an ideal CO2-H2O system
# Assumptions
# 1. Ideal gas-liquid system
# 2. CO2 conc. negligible in liquid phase, at equilibrium
# 3. Properties are applicable at standard pressure (1 atm)

# Reference Temperature : 298.15 K
# Reference Pressure : 101325 Pa

# Data Sources:

# [1] NIST Webbook, https://webbook.nist.gov/
#     Retrieved 27th November, 2020. Converted from bar to Pa
# [2] Perry's Chemical Engineers' Handbook 7th Ed.

# Initial temperature and pressure reference:
# Measurement and Modeling of the Phase Behavior of the
# (Carbon Dioxide + Water) Mixture at Temperatures From 298.15 K to 448.15
# KShu-Xin Hou, Geoffrey C. Maitland, and J. P. Martin Trusler*
# Qatar Carbonates and Carbon Storage Research Centre
# Department of Chemical Engineering Imperial College London,
# South Kensington Campus, London SW7 2AZ.  U.K


class EosType(enum.Enum):
    PR = 1
    IDEAL = 2


_phase_dicts_ideal = {
    "Vap": {"type": VaporPhase, "equation_of_state": Ideal,},
    "Liq": {"type": LiquidPhase, "equation_of_state": Ideal,},
}

_phase_dicts_pr = {
    "Vap": {
        "type": VaporPhase,
        "equation_of_state": Cubic,
        "equation_of_state_options": {"type": CubicType.PR},
    },
    "Liq": {
        "type": LiquidPhase,
        "equation_of_state": Cubic,
        "equation_of_state_options": {"type": CubicType.PR},
    },
}

# Specifying components
_component_params = {
    "H2O": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase, PhaseType.liquidPhase],
        "dens_mol_liq_comp": Perrys,
        "enth_mol_liq_comp": Perrys,
        "enth_mol_ig_comp": NIST,
        "pressure_sat_comp": NIST,
        "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
        "parameter_data": {
            "mw": (18.0153e-3, pyunits.kg / pyunits.mol),  # [1]
            "pressure_crit": (220.64e5, pyunits.Pa),  # [1]
            "temperature_crit": (647, pyunits.K),  # [1]
            "dens_mol_liq_comp_coeff": {
                # [2] pg. 2-98, temperature range 273.16 K - 333.15 K
                "eqn_type": 1,
                "1": (5.459, pyunits.kmol * pyunits.m ** -3),
                "2": (0.30542, pyunits.dimensionless),
                "3": (647.13, pyunits.K),
                "4": (0.081, pyunits.dimensionless),
            },
            "cp_mol_ig_comp_coeff": {
                # [1] temperature range 500 K- 1700 K
                "A": (30.09200, pyunits.J / pyunits.mol / pyunits.K),
                "B": (
                    6.832514,
                    pyunits.J
                    * pyunits.mol ** -1
                    * pyunits.K ** -1
                    * pyunits.kiloK ** -1,
                ),
                "C": (
                    6.793435,
                    pyunits.J
                    * pyunits.mol ** -1
                    * pyunits.K ** -1
                    * pyunits.kiloK ** -2,
                ),
                "D": (
                    -2.534480,
                    pyunits.J
                    * pyunits.mol ** -1
                    * pyunits.K ** -1
                    * pyunits.kiloK ** -3,
                ),
                "E": (
                    0.082139,
                    pyunits.J
                    * pyunits.mol ** -1
                    * pyunits.K ** -1
                    * pyunits.kiloK ** 2,
                ),
                "F": (-250.8810, pyunits.kJ / pyunits.mol),
                "G": (223.3967, pyunits.J / pyunits.mol / pyunits.K),
                "H": (0, pyunits.kJ / pyunits.mol),
            },
            "cp_mol_liq_comp_coeff": {
                # [2] pg 2-174, temperature range 273.16 K - 533.15 K
                "1": (2.7637e5, pyunits.J / pyunits.kmol / pyunits.K),
                "2": (-2.0901e3, pyunits.J / pyunits.kmol / pyunits.K ** 2),
                "3": (8.125, pyunits.J / pyunits.kmol / pyunits.K ** 3),
                "4": (-1.4116e-2, pyunits.J / pyunits.kmol / pyunits.K ** 4),
                "5": (9.3701e-6, pyunits.J / pyunits.kmol / pyunits.K ** 5),
            },
            "enth_mol_form_liq_comp_ref": (-285.83e3, pyunits.J / pyunits.mol),  # [1]
            "enth_mol_form_vap_comp_ref": (0, pyunits.J / pyunits.mol),  # [1]
            "pressure_sat_comp_coeff": {
                "A": (4.6543, None),  # [1], temperature range 255.9 K - 373 K
                "B": (1435.264, pyunits.K),
                "C": (-64.848, pyunits.K),
            },
        },
    },
    "CO2": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        "enth_mol_ig_comp": NIST,
        "parameter_data": {
            "mw": (44.0095e-3, pyunits.kg / pyunits.mol),  # [1]
            "pressure_crit": (73.825e5, pyunits.Pa),  # [1]
            "temperature_crit": (304.23, pyunits.K),  # [1]
            "cp_mol_ig_comp_coeff": {  # [1], temperature range 298 K - 1200 K
                "A": (24.99735, pyunits.J / pyunits.mol / pyunits.K),
                "B": (
                    55.18696,
                    pyunits.J
                    * pyunits.mol ** -1
                    * pyunits.K ** -1
                    * pyunits.kiloK ** -1,
                ),
                "C": (
                    -33.69137,
                    pyunits.J
                    * pyunits.mol ** -1
                    * pyunits.K ** -1
                    * pyunits.kiloK ** -2,
                ),
                "D": (
                    7.948387,
                    pyunits.J
                    * pyunits.mol ** -1
                    * pyunits.K ** -1
                    * pyunits.kiloK ** -3,
                ),
                "E": (
                    -0.136638,
                    pyunits.J
                    * pyunits.mol ** -1
                    * pyunits.K ** -1
                    * pyunits.kiloK ** 2,
                ),
                "F": (-403.6075, pyunits.kJ / pyunits.mol),
                "G": (228.2431, pyunits.J / pyunits.mol / pyunits.K),
                "H": (0, pyunits.kJ / pyunits.mol),
            },
            "enth_mol_form_vap_comp_ref": (0, pyunits.J / pyunits.mol),  # [1]
        },
    },
    "N2": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        "enth_mol_ig_comp": NIST,
        "parameter_data": {
            "mw": (0.0280134, pyunits.kg / pyunits.mol),  # [1]
            "pressure_crit": (33.9e5, pyunits.Pa),  # [1]
            "temperature_crit": (126.2, pyunits.K),  # [1]
            "cp_mol_ig_comp_coeff": {  # [1], temperature range 298 K - 1200 K
                "A": (19.50583, pyunits.J / pyunits.mol / pyunits.K),
                "B": (
                    19.88705,
                    pyunits.J
                    * pyunits.mol ** -1
                    * pyunits.K ** -1
                    * pyunits.kiloK ** -1,
                ),
                "C": (
                    -8.598535,
                    pyunits.J
                    * pyunits.mol ** -1
                    * pyunits.K ** -1
                    * pyunits.kiloK ** -2,
                ),
                "D": (
                    1.369784,
                    pyunits.J
                    * pyunits.mol ** -1
                    * pyunits.K ** -1
                    * pyunits.kiloK ** -3,
                ),
                "E": (
                    0.527601,
                    pyunits.J
                    * pyunits.mol ** -1
                    * pyunits.K ** -1
                    * pyunits.kiloK ** 2,
                ),
                "F": (-4.935202, pyunits.kJ / pyunits.mol),
                "G": (212.39, pyunits.J / pyunits.mol / pyunits.K),
                "H": (0, pyunits.kJ / pyunits.mol),
            },
            "enth_mol_form_vap_comp_ref": (0, pyunits.J / pyunits.mol),  # [1]
        },
    },
    "O2": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        "enth_mol_ig_comp": NIST,
        "parameter_data": {
            "mw": (0.031998, pyunits.kg / pyunits.mol),  # [1]
            "pressure_crit": (50.4e5, pyunits.Pa),  # [1]
            "temperature_crit": (154.6, pyunits.K),  # [1]
            "cp_mol_ig_comp_coeff": {  # [1], temperature range 298 K - 1200 K
                "A": (30.03235, pyunits.J / pyunits.mol / pyunits.K),
                "B": (
                    8.772972,
                    pyunits.J
                    * pyunits.mol ** -1
                    * pyunits.K ** -1
                    * pyunits.kiloK ** -1,
                ),
                "C": (
                    -3.988133,
                    pyunits.J
                    * pyunits.mol ** -1
                    * pyunits.K ** -1
                    * pyunits.kiloK ** -2,
                ),
                "D": (
                    0.788313,
                    pyunits.J
                    * pyunits.mol ** -1
                    * pyunits.K ** -1
                    * pyunits.kiloK ** -3,
                ),
                "E": (
                    -0.741599,
                    pyunits.J
                    * pyunits.mol ** -1
                    * pyunits.K ** -1
                    * pyunits.kiloK ** 2,
                ),
                "F": (-11.32468, pyunits.kJ / pyunits.mol),
                "G": (236.1663, pyunits.J / pyunits.mol / pyunits.K),
                "H": (0, pyunits.kJ / pyunits.mol),
            },
            "enth_mol_form_vap_comp_ref": (0, pyunits.J / pyunits.mol),  # [1]
        },
    },
    "Ar": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        "enth_mol_ig_comp": NIST,
        "parameter_data": {
            "mw": (0.039948, pyunits.kg / pyunits.mol),  # [1]
            "pressure_crit": (48.7e5, pyunits.Pa),  # [1]
            "temperature_crit": (150.8, pyunits.K),  # [1]
            "cp_mol_ig_comp_coeff": {  # [1], temperature range 298 K - 1200 K
                "A": (20.786, pyunits.J / pyunits.mol / pyunits.K),
                "B": (
                    2.825911e-07,
                    pyunits.J
                    * pyunits.mol ** -1
                    * pyunits.K ** -1
                    * pyunits.kiloK ** -1,
                ),
                "C": (
                    -1.464191e-07,
                    pyunits.J
                    * pyunits.mol ** -1
                    * pyunits.K ** -1
                    * pyunits.kiloK ** -2,
                ),
                "D": (
                    1.092131e-08,
                    pyunits.J
                    * pyunits.mol ** -1
                    * pyunits.K ** -1
                    * pyunits.kiloK ** -3,
                ),
                "E": (
                    -3.661371e-08,
                    pyunits.J
                    * pyunits.mol ** -1
                    * pyunits.K ** -1
                    * pyunits.kiloK ** 2,
                ),
                "F": (-6.19735, pyunits.kJ / pyunits.mol),
                "G": (179.999, pyunits.J / pyunits.mol / pyunits.K),
                "H": (0, pyunits.kJ / pyunits.mol),
            },
            "enth_mol_form_vap_comp_ref": (0, pyunits.J / pyunits.mol),  # [1]
        },
    },
}


# returns a configuration dictionary for the list of specified components
def get_prop(components=None, phases="Vap", eos=EosType.IDEAL, scaled=False):
    if components is None:
        components = list(_component_params.keys())
    configuration = {
        "components": {},  # fill in later based on selected components
        "parameter_data": {},
        "phases": {},
        # Set base units of measurement
        "base_units": {
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
        # Specifying state definition
        "state_definition": FTPx,
        "state_bounds": {
            "flow_mol": (0, 8000, 1e6, pyunits.mol / pyunits.s),
            "temperature": (273.15, 500, 3000, pyunits.K),
            "pressure": (5e4, 1.3e5, 1e8, pyunits.Pa),
        },
        "pressure_ref": (101325, pyunits.Pa),
        "temperature_ref": (298.15, pyunits.K),
    }

    c = configuration["components"]
    for comp in components:
        c[comp] = copy.deepcopy(_component_params[comp])
    phases = phases
    for k in phases:
        if eos == EosType.PR:
            configuration["phases"][k] = copy.deepcopy(_phase_dicts_pr[k])
        elif eos == EosType.IDEAL:
            configuration["phases"][k] = copy.deepcopy(_phase_dicts_ideal[k])
        else:
            raise ValueError("Invalid EoS.")
    if len(phases) > 1:
        p = tuple(phases)
        configuration["phases_in_equilibrium"] = [p]
        configuration["phase_equilibrium_state"] = {p: SmoothVLE}

    # Fill the binary parameters with zeros.
    d = configuration["parameter_data"]
    d["PR_kappa"] = {(a, b): 0 for a in c for b in c}

    # # Change to scaled units if specified
    # if scaled:
    #     configuration["base_units"]["mass"] = pyunits.Mg
    #     configuration["base_units"]["amount"] = pyunits.kmol

    return configuration
