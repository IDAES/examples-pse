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
This is an example flowsheet for the charge cycle of a natural gas fuel cell
(NGFC) power plant with carbon capture integrated with compressed air energy
storage. The model uses a reduced order model created by PNNL to calculate
the performance of the solid oxide fuel cell.
The power plant model builds a simple steam cycle using the the efficiency
and energy balance.
For modeling compressed air energy storage, a single stage compressor and a
single cooler for interstage cooling is assumed. A compressed gas tank model
with fixed volume is used to model gas storage.
"""

# Import Pyomo libraries
from pyomo.environ import TransformationFactory, ConcreteModel, Param
from pyomo.network import Arc
import pyomo.environ as pyo

# IDAES Imports
import idaes
from idaes.core import FlowsheetBlock, MaterialBalanceType
from idaes.models.unit_models import Heater, PressureChanger
from idaes.models.unit_models.pressure_changer import \
    ThermodynamicAssumption
from idaes.core.solvers import use_idaes_solver_configuration_defaults
from idaes.core.util.model_statistics import degrees_of_freedom

import idaes.core.util.scaling as iscale
from idaes.core.util.initialization import propagate_state
from compressed_gas_tank import CompressedGasTank
import idaes.logger as idaeslog
import sofc as SOFC
from sofc_costing import get_capital_cost, get_fixed_costs, get_variable_costs

def build_model(m):
    
    # create unti model instances
    m.fs.air_compressor = PressureChanger(
            property_package= m.fs.air_props,
            compressor= True,
            material_balance_type= MaterialBalanceType.componentTotal,
            thermodynamic_assumption= ThermodynamicAssumption.isentropic,
    )

    m.fs.interstage_cooler = Heater(
            dynamic= False,
            property_package= m.fs.air_props,
            has_pressure_change= False,
    )

    m.fs.storage_tank = CompressedGasTank(
        property_package= m.fs.air_props,
        dynamic= False
    )

    # compressor to cooler arc
    m.fs.caes_air01 = Arc(
        source=m.fs.air_compressor.outlet,
        destination=m.fs.interstage_cooler.inlet
    )
    # cooler to tank arc
    m.fs.caes_air02 = Arc(
        source=m.fs.interstage_cooler.outlet,
        destination=m.fs.storage_tank.inlet
    )
    TransformationFactory("network.expand_arcs").apply_to(m)

def set_inputs(m):
    
    # *********** COMPRESSOR INPUTS ***********
    # Air inlet to compressor
    m.fs.air_compressor.inlet.pressure[0].fix(101.325)  # kPa
    m.fs.air_compressor.inlet.temperature[0].fix(303)  # K
    m.fs.air_compressor.inlet.flow_mol.fix(1)  # kmol/s
    m.fs.air_compressor.inlet.mole_frac_comp[0, "H2O"].fix(0.0104)
    m.fs.air_compressor.inlet.mole_frac_comp[0, "CO2"].fix(0.0003)
    m.fs.air_compressor.inlet.mole_frac_comp[0, "N2"].fix(0.7726)
    m.fs.air_compressor.inlet.mole_frac_comp[0, "O2"].fix(0.2077)
    m.fs.air_compressor.inlet.mole_frac_comp[0, "Ar"].fix(0.0094)

    m.fs.air_compressor.efficiency_isentropic.fix(0.9)
    m.fs.air_compressor.outlet.pressure[0].fix(7e3)  # kPa

    # *********** COOLER INPUTS ***********
    m.fs.interstage_cooler.outlet.temperature[0].fix(323)  # K

    # *********** TANK INPUTS ***********
    m.fs.storage_tank.volume_cons.deactivate()
    m.fs.storage_tank.control_volume.volume[:].fix(300000)
    
    # Fix initial state of tank
    m.fs.storage_tank.previous_state[0].temperature.fix(333)  # K
    m.fs.storage_tank.previous_state[0].pressure.fix(4e3)  # kPa
    m.fs.storage_tank.previous_state[0].mole_frac_comp["H2O"].fix(0.0104)
    m.fs.storage_tank.previous_state[0].mole_frac_comp["N2"].fix(0.7722)
    m.fs.storage_tank.previous_state[0].mole_frac_comp["O2"].fix(0.2077)
    m.fs.storage_tank.previous_state[0].mole_frac_comp["Ar"].fix(0.0094)
    
    # Fix Duration of Operation (Time Step =  1hr)
    m.fs.storage_tank.dt[0].fix(3600)
    
    # Fix the outlet flow to zero: during charge
    m.fs.storage_tank.control_volume.properties_out[0].flow_mol.fix(0)

def add_bounds(m):

    m.fs.air_compressor.control_volume.properties_out[0].\
        pressure.setub(1e5)
    m.fs.air_compressor.properties_isentropic[0].\
        pressure.setub(1e5)

    m.fs.interstage_cooler.control_volume.properties_in[0].\
        pressure.setub(1e5)
    m.fs.interstage_cooler.control_volume.properties_out[0].\
        pressure.setub(1e5)

    # Setting the bounds on the state variables
    m.fs.storage_tank.control_volume.properties_in[0].\
        pressure.setub(1e5)
    m.fs.storage_tank.control_volume.properties_out[0].\
        pressure.setub(1e8)
    m.fs.storage_tank.previous_state[0].pressure.setub(1e8)

def initialize(m, outlvl=idaeslog.NOTSET, solver=None, optarg=None):

    # m.fs.storage_tank.display()
    iscale.calculate_scaling_factors(m)

    # initialize compressor
    m.fs.air_compressor.initialize(outlvl=outlvl, optarg=solver.options)

    # initialize cooler
    propagate_state(m.fs.caes_air01)
    m.fs.interstage_cooler.initialize(outlvl=outlvl, optarg=solver.options)

    # initialize tank
    propagate_state(m.fs.caes_air02)
    m.fs.storage_tank.previous_state[0].\
        mole_frac_phase_comp["Vap", "H2O"].value = 0.0104
    m.fs.storage_tank.previous_state[0].\
        mole_frac_phase_comp["Vap", "CO2"].value = 0.0003
    m.fs.storage_tank.previous_state[0].\
        mole_frac_phase_comp["Vap", "N2"].value = 0.7722
    m.fs.storage_tank.previous_state[0].\
        mole_frac_phase_comp["Vap", "O2"].value = 0.2077
    m.fs.storage_tank.previous_state[0].\
        mole_frac_phase_comp["Vap", "Ar"].value = 0.0094
    m.fs.storage_tank.previous_state[0].\
        flow_mol_phase["Vap"].value = 0
    m.fs.storage_tank.initialize(outlvl=idaeslog.NOTSET, optarg=solver.options)

    solver.solve(m, tee=True)
    print("*************  Storage Unit Models Initialized   **************")

def build_costing(m, solver=None):

    # Chemical engineering cost index for 2019
    m.CE_index = 607.5
    # Number of years for annulaizing capital cost
    m.number_of_years = 15

    m.fs.ref_CAES_TOC = Param(
        initialize=1168,
        doc="Reference Total Owner's Cost for CAES [$/kW]")

    m.fs.ref_CAES_fixedOM = Param(
        initialize=16.12,
        doc="Reference Total Owner's Cost for CAES [$/kW]")

    m.fs.ref_CAES_variable_OM = Param(
        initialize=0.5125,
        doc="Reference variable O&M Cost for CAES [$/MWh]")

    m.fs.factor_TOC_to_TASC = Param(
        initialize=1.093,
        doc="Multiplier to convert TOC to TASC")

    m.fs.factor_TASC_to_Annual = Param(
        initialize=0.0707,
        doc="Multiplier to convert TOC to TASC")

    m.fs.storage_P_max = Param(
        initialize=100,
        doc="Pmax for storage [MW]")

    get_variable_costs(m)

    @m.fs.Expression()
    def caes_fixed_cost(fs):
        return (
            fs.ref_CAES_TOC * 1000 * fs.storage_P_max *
            fs.factor_TOC_to_TASC * fs.factor_TASC_to_Annual * 1e-6)

    @m.fs.Expression()
    def caes_fixed_om_cost(fs):
        return fs.ref_CAES_fixedOM * 1000 * fs.storage_P_max / 1e6

    @m.fs.Expression()
    def sofc_variable_fuel_cost(fs):
        return fs.costing.variable_operating_costs[0, "natural gas"]

    @m.fs.Expression()
    def sofc_variable_nonfuel_cost(fs):
        return (fs.costing.total_variable_OM_cost[0] -
                fs.costing.variable_operating_costs[0, "natural gas"])

    @m.fs.Expression()
    def caes_variable_cost(fs):
        return fs.ref_CAES_variable_OM

    @m.fs.Expression()
    def total_fixed_cost(fs):
        return fs.caes_fixed_cost + fs.caes_fixed_om_cost

    @m.fs.Expression()
    def total_variable_cost(fs):
        return (fs.caes_variable_cost +
                fs.sofc_variable_fuel_cost +
                fs.sofc_variable_nonfuel_cost)

    @m.fs.Expression()
    def total_variable_nonfuel_cost(fs):
        return fs.caes_variable_cost + fs.sofc_variable_nonfuel_cost

    solver.solve(m, tee=True)

    print("*************  Storage Costing Initialized   **************")

def integrate_storage(m):

    # Deactivating and writing new constraints to integrate storage
    del m.fs.steam_cycle_heat_duty
    del m.fs.net_power

    @m.fs.Expression(m.fs.time)
    def charge_rate_MW(fs, t):
        return m.fs.air_compressor.control_volume.work[t] * 1e-3

    @m.fs.Expression(m.fs.time)
    def steam_cycle_heat_duty(fs, t):
        return (fs.hrsg_heat_duty[t] -
                fs.asu_steam_duty[t] +
                fs.interstage_cooler.heat_duty[0] * 1e-3)

    @m.fs.Expression(m.fs.time)
    def net_power(fs, t):
        return (fs.gross_power[t] -
                fs.auxiliary_load[t] -
                m.fs.air_compressor.control_volume.work[t] * 1e-3)

    solver.solve(m, tee=True)

    print("************  Initialized after integrating storage  *************")

def build_caes_charge(m, solver=None):

    # build the storage model, set model inputs, and add bounds
    build_model(m)
    set_inputs(m)
    add_bounds(m)

    # check if the storage model is complete by asserting DOF = 0
    assert degrees_of_freedom(m) == 0

    # initialize the charge storage model
    initialize(m, solver=solver)

    return m

if __name__ == "__main__":

    use_idaes_solver_configuration_defaults()
    idaes.cfg.ipopt.options.nlp_scaling_method = "user-scaling"
    idaes.cfg.ipopt.options.linear_solver = "ma57"
    idaes.cfg.ipopt.options.OF_ma57_automatic_scaling = "yes"
    idaes.cfg.ipopt.options.ma57_pivtol = 1e-5
    idaes.cfg.ipopt.options.ma57_pivtolmax = 0.1
    idaes.cfg.ipopt.options.bound_push = 1e-20
    solver = pyo.SolverFactory("ipopt")

    # # create charge model
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic= False)

    # build SOFC model
    m = SOFC.get_model(m)
    SOFC.initialize(m)

    # build storage model
    m = build_caes_charge(m, solver=solver)

    # add costing
    build_costing(m, solver=solver)

    # initialize unit models in storage
    integrate_storage(m)
    

