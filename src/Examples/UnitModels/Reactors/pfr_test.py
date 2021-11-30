# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 14:24:31 2021

@author: bpaul
"""

from pyomo.environ import (Constraint,
                           ConstraintList,
                           exp,
                           log,
                           Var,
                           ConcreteModel,
                           Expression,
                           Objective,
                           TransformationFactory,
                           value,
                           units as pyunits)
from pyomo.dae import DerivativeVar, ContinuousSet, Integral
from pyomo.network import Arc

from idaes.core import FlowsheetBlock
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)
from idaes.generic_models.properties.core.generic.generic_reaction import (
        GenericReactionParameterBlock)
from idaes.generic_models.unit_models import (Mixer,
                                              Heater,
                                              PFR)

from idaes.core.util.constants import Constants
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state

import idaes.logger as idaeslog

import egprod_ideal as thermo_props
import egprod_reaction as reaction_props

m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
m.fs.thermo_params = GenericParameterBlock(default=thermo_props.config_dict)
m.fs.reaction_params = GenericReactionParameterBlock(default={"property_package": m.fs.thermo_params,
                                                              **reaction_props.config_dict})
m.fs.M101 = Mixer(default={"property_package": m.fs.thermo_params,
                           "inlet_list": ["reagent_feed", "catalyst_feed"]})
m.fs.H101 = Heater(default={"property_package": m.fs.thermo_params,
                            "has_pressure_change": False,
                            "has_phase_equilibrium": False})
m.fs.R101 = PFR(
            default={"property_package": m.fs.thermo_params,
                     "reaction_package": m.fs.reaction_params,
                     "has_heat_of_reaction": True,
                     "has_heat_transfer": True,
                     "has_pressure_change": False,
                     "transformation_method": "dae.finite_difference",
                     "transformation_scheme": "BACKWARD",
                     "finite_elements": 1,
                     "length_domain_set": [0, 1]})
m.fs.s03 = Arc(source=m.fs.M101.outlet, destination=m.fs.H101.inlet)
m.fs.s04 = Arc(source=m.fs.H101.outlet, destination=m.fs.R101.inlet)
TransformationFactory("network.expand_arcs").apply_to(m)

m.fs.eg_prod = Expression(expr=pyunits.convert(m.fs.R101.outlet.flow_mol_phase_comp[0, "Liq", "ethylene_glycol"]
                                               *m.fs.thermo_params.ethylene_glycol.mw, # MW defined in properties as kg/mol
                                               to_units=pyunits.Mlb/pyunits.yr)) # converting kg/s to MM lb/year
m.fs.cooling_cost = Expression(expr=0.212e-7 * (-sum(m.fs.R101.heat_duty[0, x]
                                                     for x in m.fs.R101.control_volume.length_domain)))  # the reaction is exothermic, so R101 duty is negative
m.fs.heating_cost = Expression(expr=2.2e-7 * m.fs.H101.heat_duty[0])  # the stream must be heated to T_rxn, so H101 duty is positive
m.fs.operating_cost = Expression(expr=(3600 * 8000 *(m.fs.heating_cost + m.fs.cooling_cost)))

print('DOF total: ', str(degrees_of_freedom(m)))
# Check the degrees of freedom
assert degrees_of_freedom(m) == 16

m.fs.M101.reagent_feed.flow_mol_phase_comp[0, "Liq", "ethylene_oxide"].fix(58.0*pyunits.mol/pyunits.s)
m.fs.M101.reagent_feed.flow_mol_phase_comp[0, "Liq", "water"].fix(39.6*pyunits.mol/pyunits.s)  # calculated from 16.1 mol EO / cudm in stream
m.fs.M101.reagent_feed.flow_mol_phase_comp[0, "Liq", "sulfuric_acid"].fix(1e-5*pyunits.mol/pyunits.s)
m.fs.M101.reagent_feed.flow_mol_phase_comp[0, "Liq", "ethylene_glycol"].fix(1e-5*pyunits.mol/pyunits.s)
m.fs.M101.reagent_feed.temperature.fix(298.15*pyunits.K)
m.fs.M101.reagent_feed.pressure.fix(1e5*pyunits.Pa)

m.fs.M101.catalyst_feed.flow_mol_phase_comp[0, "Liq", "ethylene_oxide"].fix(1e-5*pyunits.mol/pyunits.s)
m.fs.M101.catalyst_feed.flow_mol_phase_comp[0, "Liq", "water"].fix(200*pyunits.mol/pyunits.s)
m.fs.M101.catalyst_feed.flow_mol_phase_comp[0, "Liq", "sulfuric_acid"].fix(0.334*pyunits.mol/pyunits.s)  # calculated from 0.9 wt% SA in stream
m.fs.M101.catalyst_feed.flow_mol_phase_comp[0, "Liq", "ethylene_glycol"].fix(1e-5*pyunits.mol/pyunits.s)
m.fs.M101.catalyst_feed.temperature.fix(298.15*pyunits.K)
m.fs.M101.catalyst_feed.pressure.fix(1e5*pyunits.Pa)

m.fs.H101.outlet.temperature.fix(328.15*pyunits.K)

m.fs.R101.conversion = Var(m.fs.R101.control_volume.length_domain, bounds=(0,1),
                           initialize=0.80, units=pyunits.dimensionless)  # fraction

def conv_rule(m, x):  # rule for calculating cumulative conversion at each finite element
    return (m.fs.R101.conversion[x] * m.fs.R101.inlet.flow_mol_phase_comp[0, "Liq", "ethylene_oxide"] ==
            (m.fs.R101.inlet.flow_mol_phase_comp[0, "Liq", "ethylene_oxide"] -
            m.fs.R101.control_volume.properties[0, x].flow_mol_phase_comp["Liq", "ethylene_oxide"]))
m.conv_constraint = Constraint(m.fs.R101.control_volume.length_domain, rule=conv_rule)

def rate_constant(m, x):  # rule for calculating reaction rate constant at every finite element
    return (m.fs.reaction_params.reaction_R1.arrhenius_const * exp(-m.fs.reaction_params.reaction_R1.energy_activation /
           (Constants.gas_constant*m.fs.R101.control_volume.properties[0, x].temperature)))  # r = k[EO], k = [/s]

#m.fs.R101.dconv_dlength = DerivativeVar(m.fs.R101.conversion, wrt=m.fs.R101.control_volume.length_domain)

#def forward_diff(m, idx):  # forward difference discretization for derivative variable, idx starts at 0 but x starts at 1
#    return (m.fs.R101.conversion[m.fs.R101.control_volume.length_domain[idx + 2]] ==
#            m.fs.R101.conversion[m.fs.R101.control_volume.length_domain[idx + 1]] +
#            m.fs.R101.dconv_dlength[m.fs.R101.control_volume.length_domain[idx + 1]] *
#            (m.fs.R101.control_volume.length_domain[idx + 2] - m.fs.R101.control_volume.length_domain[idx + 1]))
#m.forward_diff_approx = Constraint(range(len(m.fs.R101.control_volume.length_domain) - 1), rule=forward_diff)

def perf_eqn_rule(m, idx):  # rule for PFR performance equation at each finite element
    return (m.fs.R101.dconv_dlength[x] == m.fs.R101.control_volume.area * (rate_constant(m, x) /
                                                                           (2*3.62e-3)) * (1 - m.fs.R101.conversion[x]))
#m.perf_eqn_constraints = Constraint(m.fs.R101.control_volume.length_domain, rule=perf_eqn_rule)

m.fs.R101.control_volume.length.fix(1*pyunits.m)
m.fs.R101.volume.fix(5.538*pyunits.m**3)

print('DOF unspecified: ', str(degrees_of_freedom(m)))
# Check the degrees of freedom
assert degrees_of_freedom(m) == 0

# Initialize and solve each unit operation
m.fs.M101.initialize()
propagate_state(arc=m.fs.s03)

m.fs.H101.initialize()
propagate_state(arc=m.fs.s04)

m.fs.R101.initialize()

# set solver
solver = get_solver()

results = solver.solve(m, tee=True)