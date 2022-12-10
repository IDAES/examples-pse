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
This is an example flowsheet for the discharge cycle of a solid oxide fuel cell
(SOFC) power plant with carbon capture integrated with a compressed air energy
storage. The model uses a reduced order model created by PNNL to calculate
the performance of the solid oxide fuel cell.
During the discharge cycle, the steam cycle is turned off.
A compressed gas tank model with fixed volume is used to model gas storage.
"""

import logging
import pandas as pd
import pyomo.environ as pyo
import idaes.core.util as iutil
import idaes.core.util.tables as tables

# Import Pyomo libraries
from pyomo.environ import TransformationFactory, ConcreteModel, Param
from pyomo.network import Arc
from pyomo.opt import TerminationCondition

# IDAES Imports
import idaes
from idaes.core import FlowsheetBlock, MaterialBalanceType
from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale
from idaes.core.solvers import use_idaes_solver_configuration_defaults
from idaes.models.unit_models import (
    Heater,
    HeatExchanger,
    PressureChanger,
    Valve,
    ValveFunctionType)
from idaes.models.unit_models.pressure_changer import \
    ThermodynamicAssumption
from idaes.models.unit_models.heat_exchanger import (
    delta_temperature_underwood_callback)
import idaes.logger as idaeslog
import sofc as SOFC
from compressed_gas_tank import CompressedGasTank
from sofc_costing import get_variable_costs
logging.basicConfig(level=logging.INFO)


def _build_model(m):

    # create unti model instances
    m.fs.air_turbine = PressureChanger(
            property_package= m.fs.air_props,
            compressor= False,
            material_balance_type= MaterialBalanceType.componentTotal,
            thermodynamic_assumption= ThermodynamicAssumption.isentropic,
    )

    m.fs.hrsg_heater = Heater(
            dynamic= False,
            property_package= m.fs.air_props,
            has_pressure_change= True
    )

    m.fs.air_preheater = HeatExchanger(
        hot_side_name="shell",
        cold_side_name="tube",
        shell={"property_package": m.fs.air_props,
               "material_balance_type": MaterialBalanceType.componentTotal,
               "has_pressure_change": True},
        tube={"property_package": m.fs.air_props,
              "material_balance_type": MaterialBalanceType.componentTotal,
              "has_pressure_change": True},
        delta_temperature_callback=delta_temperature_underwood_callback
    )

    m.fs.storage_tank = CompressedGasTank(
            property_package= m.fs.air_props,
            dynamic= False
    )

    m.fs.flow_valve = Valve(
            valve_function_callback= ValveFunctionType.linear,
            property_package= m.fs.air_props,
    )
    # deactivate pressure-flow correlation to control pressure and flow
    m.fs.flow_valve.pressure_flow_equation.deactivate()

    # create arcs

    # tank to valve
    m.fs.tank_to_valve = Arc(
        source=m.fs.storage_tank.outlet,
        destination=m.fs.flow_valve.inlet
    )
    # valve to preheater tube
    m.fs.caes_air03 = Arc(
        source=m.fs.flow_valve.outlet,
        destination=m.fs.air_preheater.tube_inlet
    )
    # preheater tube to HRSG
    m.fs.caes_air04 = Arc(
        source=m.fs.air_preheater.tube_outlet,
        destination=m.fs.hrsg_heater.inlet
    )
    # HRSG to turbine
    m.fs.caes_air05 = Arc(
        source=m.fs.hrsg_heater.outlet,
        destination=m.fs.air_turbine.inlet
    )
    # Turbine to preheater shell
    m.fs.caes_air06 = Arc(
        source=m.fs.air_turbine.outlet,
        destination=m.fs.air_preheater.shell_inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def _set_inputs(m):

    # *********** Preheater INPUTS ***********
    m.fs.air_preheater.area.fix(5000)  # m2
    m.fs.air_preheater.overall_heat_transfer_coefficient.fix(0.025)  # kW/m2/K
    m.fs.air_preheater.tube.deltaP[0].fix(-10)  # 0.1 bar = 100 kPa
    m.fs.air_preheater.shell.deltaP[0].fix(-10)  # 0.1 bar = 100 kPa

    # *********** HRSG INPUTS ***********
    m.fs.hrsg_heater.control_volume.properties_out[0].temperature.fix(1173)  # K
    m.fs.hrsg_heater.deltaP.fix(10)  # 0.1 bar = 10 kPa

    # *********** Turbine INPUTS ***********
    m.fs.air_turbine.control_volume.properties_out[0].pressure.fix(150)  # kPa
    m.fs.air_turbine.efficiency_isentropic.fix(0.9)

    # *********** Tank INPUTS ***********
    m.fs.storage_tank.volume_cons.deactivate()
    m.fs.storage_tank.control_volume.volume[:].fix(300000)
    m.fs.storage_tank.control_volume.properties_in[0].flow_mol.fix(0)
    m.fs.storage_tank.control_volume.properties_out[0].flow_mol.fix(1) # kmol/s

    # Fix initial state of tank
    m.fs.storage_tank.previous_state[0].temperature.fix(333)  # K
    m.fs.storage_tank.previous_state[0].pressure.fix(7e3)  # kPa
    m.fs.storage_tank.previous_state[0].mole_frac_comp['H2O'].fix(0.0104)
    m.fs.storage_tank.previous_state[0].mole_frac_comp['N2'].fix(0.7722)
    m.fs.storage_tank.previous_state[0].mole_frac_comp['O2'].fix(0.2077)
    m.fs.storage_tank.previous_state[0].mole_frac_comp['Ar'].fix(0.0094)

    # Fix initial state of tank
    m.fs.storage_tank.control_volume.properties_in[0].temperature.fix(333)  # K
    m.fs.storage_tank.control_volume.properties_in[0].pressure.fix(7e3)  # kPa
    m.fs.storage_tank.control_volume.properties_in[0].\
        mole_frac_comp['H2O'].fix(0.0104)
    m.fs.storage_tank.control_volume.properties_in[0].\
        mole_frac_comp['CO2'].fix(0.0003)
    m.fs.storage_tank.control_volume.properties_in[0].\
        mole_frac_comp['N2'].fix(0.7722)
    m.fs.storage_tank.control_volume.properties_in[0].\
        mole_frac_comp['O2'].fix(0.2077)
    m.fs.storage_tank.control_volume.properties_in[0].\
        mole_frac_comp['Ar'].fix(0.0094)

    # Fix Duration of Operation (Time Step =  1hr)
    m.fs.storage_tank.dt[0].fix(3600)

    # Fix outlet pressure for Valve
    m.fs.flow_valve.outlet.pressure[0].fix(4e3)


def _add_bounds(m):

    # Add required bounds to the variables
    m.fs.air_turbine.control_volume.properties_in[0].pressure.setub(1e5)
    m.fs.air_turbine.control_volume.properties_out[0].pressure.setub(1e5)

    m.fs.hrsg_heater.control_volume.properties_in[0].pressure.setub(1e5)
    m.fs.hrsg_heater.control_volume.properties_out[0].pressure.setub(1e5)

    m.fs.air_preheater.tube.properties_in[0].pressure.setub(1e5)
    m.fs.air_preheater.tube.properties_out[0].pressure.setub(1e5)

    m.fs.flow_valve.control_volume.properties_in[0].pressure.setub(1e5)
    m.fs.flow_valve.control_volume.properties_out[0].pressure.setub(1e5)

    # Setting the bounds on the state variables
    m.fs.storage_tank.control_volume.properties_in[0].pressure.setub(1e5)
    m.fs.storage_tank.control_volume.properties_out[0].pressure.setub(1e5)
    m.fs.storage_tank.previous_state[0].pressure.setub(1e5)


def _initialize_storage(m, outlvl=idaeslog.NOTSET, solver=None):

    iscale.calculate_scaling_factors(m)

    # initialize tank
    solver.solve(m.fs.storage_tank)
    m.fs.storage_tank.initialize(outlvl=outlvl, optarg=solver.options)

    # initialize valve
    propagate_state(m.fs.tank_to_valve)
    m.fs.flow_valve.control_volume.properties_in[0].flow_mol.value = 6
    m.fs.flow_valve.initialize(outlvl=outlvl, optarg=solver.options)

    # initialize preheater
    m.fs.air_preheater.shell.properties_in[0].\
        flow_mol.value = 1  # kmol/s
    m.fs.air_preheater.shell.properties_in[0].\
        temperature.value = 600  # K
    m.fs.air_preheater.shell.properties_in[0].\
        pressure.value = 1740  # kPa
    m.fs.air_preheater.shell.properties_in[0].\
        mole_frac_comp['H2O'].value = 0.0104
    m.fs.air_preheater.shell.properties_in[0].\
        mole_frac_comp['CO2'].value = 0.0003
    m.fs.air_preheater.shell.properties_in[0].\
        mole_frac_comp['N2'].value = 0.7722
    m.fs.air_preheater.shell.properties_in[0].\
        mole_frac_comp['O2'].value = 0.2077
    m.fs.air_preheater.shell.properties_in[0].\
        mole_frac_comp['Ar'].value = 0.0094
    propagate_state(m.fs.caes_air03)
    m.fs.air_preheater.initialize(outlvl=outlvl, optarg=solver.options)

    # initialize HRSG
    propagate_state(m.fs.caes_air04)
    m.fs.hrsg_heater.initialize(outlvl=outlvl, optarg=solver.options)

    # initialize turbine
    propagate_state(m.fs.caes_air05)
    m.fs.air_turbine.initialize(outlvl=outlvl, optarg=solver.options)

    # update conditions for discharge scenario
    m.fs.storage_tank.control_volume.properties_in[0].flow_mol.fix(0)  # kmol/s

    init_results = solver.solve(m, tee=False)
    print("CAES unit model initialization solver status:",
          init_results.solver.termination_condition)
    print("*************  Storage Unit Models Initialized   **************")


def _build_fs_costing(m, solver=None):

    # Chemical engineering cost index for 2019
    m.CE_index = 607.5
    # Number of years for annulaizing capital cost
    m.number_of_years = 15

    get_variable_costs(m)

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

    costing_results = solver.solve(m, tee=False)
    print("Model initialization with cost constraints, solver status:",
          costing_results.solver.termination_condition)
    print("*************  Storage Costing Initialized   **************")


def _integrate_storage_with_sofc(m, solver=None):

    # Deactivating constraints to turn off steams cycle and integrate storage
    del m.fs.steam_cycle_heat_duty
    del m.fs.gross_power

    @m.fs.Expression(m.fs.time)
    def discharge_rate_MW(fs, t):
        return (-1 * m.fs.air_turbine.control_volume.work[t] * 1e-3)

    @m.fs.Expression(m.fs.time)
    def gross_power(fs, t):
        return (fs.sofc_power_ac[t]
                + fs.air_turbine.control_volume.work[t]*-1e-3
                )

    # HRSG heat applied toward steam cycle in MMBtu/hr
    @m.fs.Expression(m.fs.time)
    def steam_cycle_heat_duty(fs, t):
        return (fs.hrsg_heat_duty[t]
                - fs.asu_steam_duty[t]
                - fs.hrsg_heater.heat_duty[0]*1e-3
                )

    ies_init_results = solver.solve(m, tee=False)
    print("Integrated SOFC + CAES model initialization, solver status:",
          ies_init_results.solver.termination_condition)
    print("************  Initialized after integrating storage  *************")


def build_model_sofc_caes_discharge(m, solver=None):

    # build SOFC model
    m = SOFC.get_model(m)
    SOFC.initialize(m)

    # build the storage model, set model inputs, and add bounds
    _build_model(m)
    _set_inputs(m)
    _add_bounds(m)

    # check if the storage model is complete by asserting DOF = 0
    assert degrees_of_freedom(m) == 0

    # initialize the charge storage model
    _initialize_storage(m, solver=solver)
    assert degrees_of_freedom(m) == 0

    # integrate storage with SOFC flowsheet
    _integrate_storage_with_sofc(m, solver=solver)

    # add costing
    _build_fs_costing(m, solver=solver)

    # setup discharge operational constraints
    m.fs.storage_tank.control_volume.properties_out[0].flow_mol.unfix()  # kmol/s
    m.fs.flow_valve.control_volume.properties_out[0].flow_mol.fix(6)  # kmol/s

    return m

def _add_tags(m):
    tag_group = iutil.ModelTagGroup()
    m._tags_streams = tag_group
    stream_states = tables.stream_states_dict(
        tables.arcs_to_stream_dict(
            m.fs,
            descend_into=False,
            additional={
                "ng00": m.fs.anode_mix.feed_inlet,
                "air_f00": m.fs.air_blower.inlet,
                "asu00": m.fs.air_compressor_s1.inlet,
                "n2": m.fs.ASU.N2_outlet,
                "air_f01": m.fs.cathode_HRSG.outlet,
                "caes_air07": m.fs.air_preheater.outlet_1,
                
            }
        )
    )
    for i, s in stream_states.items():  # create the tags for steam quantities
        tag_group[f"{i}_Fmass"] = iutil.ModelTag(
            doc=f"{i}: mass flow",
            expr=s.flow_mass,
            format_string="{:.3f}",
            display_units=pyo.units.kg / pyo.units.s,
        )
        tag_group[f"{i}_F"] = iutil.ModelTag(
            doc=f"{i}: mole flow",
            expr=s.flow_mol,
            format_string="{:.3f}",
            display_units=pyo.units.kmol / pyo.units.s,
        )
        for c in s.mole_frac_comp:
            tag_group[f"{i}_y{c}"] = iutil.ModelTag(
                doc=f"{i}: mole percent {c}",
                expr=s.mole_frac_comp[c] * 100,
                format_string="{:.3f}",
                display_units="%",
            )
        tag_group[f"{i}_P"] = iutil.ModelTag(
            doc=f"{i}: pressure",
            expr=s.pressure,
            format_string="{:.3f}",
            display_units=pyo.units.bar,
        )
        tag_group[f"{i}_T"] = iutil.ModelTag(
            doc=f"{i}: temperature",
            expr=s.temperature,
            format_string="{:.2f}",
            display_units=pyo.units.K,
        )
        tag_group[f"{i}_Fvol"] = iutil.ModelTag(
            doc=f"{i}: volumetric flow",
            expr=s.flow_vol,
            format_string="{:.3f}",
            display_units=pyo.units.m**3 / pyo.units.s,
        )
        if hasattr(s, "vapor_frac"):
            vf = 100 * s.vapor_frac
        else:
            vf = 100
        tag_group[f"{i}_vf"] = iutil.ModelTag(
            doc=f"{i}: vapor fraction",
            expr=vf,
            format_string="{:.2f}",
            display_units="%",
        )

    other_ports = {"wat01": m.fs.flash.liq_outlet,
                   "fg04": m.fs.flash.vap_outlet,
                   "wat02": m.fs.cpu.water,
                   "co2": m.fs.cpu.pureco2,
                   "vent": m.fs.cpu.vent}

    for i, p in other_ports.items():
        tag_group[f"{i}_F"] = iutil.ModelTag(
            doc=f"{i}: mole flow",
            expr=p.flow_mol[0],
            format_string="{:.3f}",
            display_units=pyo.units.kmol / pyo.units.s,
        )
        for c in p.mole_frac_comp:
            j = c[1]
            tag_group[f"{i}_y{j}"] = iutil.ModelTag(
                doc=f"{i}: mole percent {j}",
                expr=p.mole_frac_comp[c] * 100,
                format_string="{:.3f}",
                display_units="%",
            )
        tag_group[f"{i}_P"] = iutil.ModelTag(
            doc=f"{i}: pressure",
            expr=p.pressure[0],
            format_string="{:.3f}",
            display_units=pyo.units.bar,
        )
        tag_group[f"{i}_T"] = iutil.ModelTag(
            doc=f"{i}: temperature",
            expr=p.temperature[0],
            format_string="{:.2f}",
            display_units=pyo.units.K,
        )


def _stream_col_gen(tag_group):
    for tag in tag_group.values():
        spltstr = tag.doc.split(":")
        stream = spltstr[0].strip()
        col = f"{spltstr[1].strip()} ({tag.get_unit_str()})"
        yield tag, stream, col

def _stream_table(tag_group):
    rows = set()
    cols = set()
    tags = []
    for tag, stream, col in _stream_col_gen(tag_group):
        rows.add(stream)
        cols.add(col)
        tags.append((tag, stream, col))
    df = pd.DataFrame(index=sorted(rows), columns=sorted(cols))
    for tag, stream, col in tags:
        df.at[stream, col] = tag.get_display_value()
    return df

def streams_dataframe(m):
    df = _stream_table(m.tags_streams)
    return df


if __name__ == "__main__":

    use_idaes_solver_configuration_defaults()
    idaes.cfg.ipopt.options.nlp_scaling_method = "user-scaling"
    idaes.cfg.ipopt.options.linear_solver = "ma57"
    idaes.cfg.ipopt.options.OF_ma57_automatic_scaling = "yes"
    idaes.cfg.ipopt.options.ma57_pivtol = 1e-5
    idaes.cfg.ipopt.options.ma57_pivtolmax = 0.1
    idaes.cfg.ipopt.options.bound_push = 1e-20
    solver = pyo.SolverFactory("ipopt")

    # create charge model
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic= False)

    # build storage model
    m = build_model_sofc_caes_discharge(m, solver=solver)
    results = solver.solve(m, tee=True, symbolic_solver_labels=True)
