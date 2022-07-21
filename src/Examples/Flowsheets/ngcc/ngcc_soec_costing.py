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
Cost evaluation using IDAES Costing Framework for EMRE NGCC Model
Base Case : B31B - NETL Baseline Report Rev 4
Author: A. Deshpande, Alex Noring, M. Zamarripa
"""

# NGCC cost evaluation using IDAES Costing Framework

import pytest
import idaes
from idaes.models_extra.power_generation.costing.power_plant_costing import (
    get_PP_costing,
    get_total_TPC,
    costing_initialization,
    get_fixed_OM_costs,
    get_variable_OM_costs,
    initialize_fixed_OM_costs,
    initialize_variable_OM_costs,
)
from idaes.core.util.model_statistics import degrees_of_freedom
import pyomo.environ as pyo
from pyomo.environ import units as pyunits
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core.solvers import use_idaes_solver_configuration_defaults
import idaes.core.util.scaling as iscale

import ngcc_soec

from idaes.core.util.unit_costing import initialize as initialize_unit_costing


def add_total_plant_cost(b, installation_labor=1.25, eng_fee=1.2,
                         project_contingency=1.15, process_contingency=1.15):

    @b.costing.Expression()
    def total_plant_cost(c):
        return (c.purchase_cost * installation_labor * eng_fee *
                (project_contingency + process_contingency) / 1e6)


def add_results_for_costing(m):
    ###########################################################################
    # CO2 CAPTURE SYSTEM

    baseline_capture = 493487*pyunits.lb/pyunits.hr
    capture_percentage = 0.97

    @m.fs.Expression(m.fs.time)
    def CO2_captured(b, t):
        return pyunits.convert(
            (b.ngcc.gt.gts2.control_volume.properties_out[t].flow_mol_comp["CO2"] *
             0.04401*pyunits.kg/pyunits.mol * capture_percentage),
            pyunits.lb/pyunits.hr)

    # CO2 purification, drying, and compression cooling duty
    @m.fs.Expression(m.fs.time)
    def CO2_cooling_duty(b, t):
        baseline_duty = 1394.4*pyunits.MBtu/pyunits.hr
        return baseline_duty*(b.CO2_captured[t]/baseline_capture)

    # CO2 capture system makeup water
    @m.fs.Expression(m.fs.time)
    def CO2_water_makeup(b, t):
        baseline_makeup = 55.253*pyunits.gal/pyunits.min
        return baseline_makeup*(b.CO2_captured[t]/baseline_capture)

    #  CO2 capture system water discharge
    @m.fs.Expression(m.fs.time)
    def CO2_water_discharge(b, t):
        baseline_discharge = 604.73*pyunits.gal/pyunits.min
        return baseline_discharge*(b.CO2_captured[t]/baseline_capture)

    # CO2 capture auxiliary load
    @m.fs.Expression(m.fs.time)
    def CO2_aux_load(b, t):
        baseline_CO2_aux_load = 27690*pyunits.kW
        return baseline_CO2_aux_load*(b.CO2_captured[t]/baseline_capture)

    # CO2 compressor auxiliary load
    # this load is included in CO2_aux_load but needed seperatly for costing
    @m.fs.Expression(m.fs.time)
    def CO2_compressor_load(b, t):
        baseline_compressor_load = 17090*pyunits.kW
        return baseline_compressor_load*(b.CO2_captured[t]/baseline_capture)

    # CO2 aftercooler duty, MMBtu/hr
    @m.fs.Expression(m.fs.time)
    def CO2_aftercooler_duty(b, t):
        baseline_aftercooler_duty = 32.24*pyunits.MBtu/pyunits.hr
        return baseline_aftercooler_duty*(b.CO2_captured[t]/baseline_capture)

    ###########################################################################
    # COOLING WATER SYSTEM

    water_density = 8.333*pyunits.lb/pyunits.gal
    water_heat_capacity = 1*pyunits.Btu/pyunits.lb/pyunits.F
    CW_range = 20*pyunits.F
    cycles_of_concentration = 4
    baseline_cw_flow = 220784*pyunits.gal/pyunits.min

    # steam cycle condenser duty
    @m.fs.Expression(m.fs.time)
    def condenser_duty(b, t):
        return pyunits.convert(
            b.ngcc.st.main_condenser.heat_duty[t],
            pyunits.MBtu/pyunits.hr)

    # hydrogen compressor cooling duty
    @m.fs.Expression(m.fs.time)
    def h2_comp_cooling_duty(b, t):
        return -1*pyunits.convert(
            (b.soec.ic01.heat_duty[t] +
             b.soec.ic02.heat_duty[t] +
             b.soec.ic03.heat_duty[t]),
            pyunits.MBtu/pyunits.hr)

    # cooling water duty, includes 25 MMBtu of misc cooling loads
    @m.fs.Expression(m.fs.time)
    def cooling_water_duty(b, t):
        return (b.condenser_duty[t] +
                b.h2_comp_cooling_duty[t] +
                b.CO2_cooling_duty[t] +
                25*pyunits.MBtu/pyunits.hr)

    # cooling water circulation rate
    @m.fs.Expression(m.fs.time)
    def cooling_water_flow(b, t):
        return pyunits.convert(
            b.cooling_water_duty[t]/water_heat_capacity/CW_range/water_density,
            pyunits.gal/pyunits.min)

    # cooling tower evaporation losses
    @m.fs.Expression(m.fs.time)
    def evaporation_losses(b, t):
        return 0.008*20/10*b.cooling_water_flow[t]

    # cooling tower blowdown losses
    @m.fs.Expression(m.fs.time)
    def blowdown_losses(b, t):
        return b.evaporation_losses[t]/(cycles_of_concentration-1)

    # cooling water pump auxiliary load
    @m.fs.Expression(m.fs.time)
    def cooling_water_pump_load(b, t):
        baseline_pump_load = 4580*pyunits.kW
        return baseline_pump_load*(b.cooling_water_flow[t]/baseline_cw_flow)

    # cooling tower fans auxiliary load
    @m.fs.Expression(m.fs.time)
    def cooling_tower_fan_load(b, t):
        baseline_fan_load = 2370*pyunits.kW
        return baseline_fan_load*(b.cooling_water_flow[t]/baseline_cw_flow)

    ###########################################################################
    # SYSTEM WATER USAGE

    # BFW makeup
    @m.fs.Expression(m.fs.time)
    def BFW_makeup(b, t):
        return pyunits.convert(
            b.ngcc.st.hotwell.makeup_state[t].flow_vol,
            pyunits.gal/pyunits.min)

    # SOEC makeup
    @m.fs.Expression(m.fs.time)
    def SOEC_makeup(b, t):
        return pyunits.convert(
            b.soec.water_pump.control_volume.properties_in[t].flow_vol,
            pyunits.gal/pyunits.min)

    # water demand
    @m.fs.Expression(m.fs.time)
    def water_demand(b, t):
        return (b.evaporation_losses[t] +
                b.blowdown_losses[t] +
                b.BFW_makeup[t] +
                b.SOEC_makeup[t] +
                b.CO2_water_makeup[t])

    # internal recycle
    @m.fs.Expression(m.fs.time)
    def water_recycle(b, t):
        return 0*pyunits.gal/pyunits.min

    # raw water withdrawal
    @m.fs.Expression(m.fs.time)
    def raw_water_withdrawal(b, t):
        return b.water_demand[t] - b.water_recycle[t]

    # hydrogen dryer discharge
    @m.fs.Expression(m.fs.time)
    def dryer_discharge(b, t):
        return pyunits.convert(
            (b.soec.water_heater01.shell_outlet.flow_mol[t] *
             b.soec.water_heater01.shell_outlet.mole_frac_comp[t, "H2O"] *
             18.02*pyunits.g/pyunits.mol / water_density),
            pyunits.gal/pyunits.min)

    # process water discharge
    @m.fs.Expression(m.fs.time)
    def process_water_discharge(b, t):
        return (b.blowdown_losses[t]*0.9 +
                b.CO2_water_discharge[t] +
                b.BFW_makeup[t] +
                b.dryer_discharge[t])

    # raw water consumption
    @m.fs.Expression(m.fs.time)
    def raw_water_consumption(b, t):
        return b.raw_water_withdrawal[t] - b.process_water_discharge[t]

    # ground water pump auxiliary load
    @m.fs.Expression(m.fs.time)
    def ground_water_pump_load(b, t):
        baseline_pump_load = 430*pyunits.kW
        baseline_water_withdrawal = 4773.14*pyunits.gal/pyunits.min
        return baseline_pump_load*(b.raw_water_withdrawal[t]/baseline_water_withdrawal)

    ###########################################################################
    # AUXILIARY LOADS AND NET POWER

    baseline_gross_power = 689800*pyunits.kW

    # gas turbine gross power
    @m.fs.Expression(m.fs.time)
    def gas_turbine_gross_power(b, t):
        return pyunits.convert(
            -1*b.ngcc.gt.gt_power[t],
            pyunits.kW)

    # steam turbine gross power
    @m.fs.Expression(m.fs.time)
    def steam_turbine_gross_power(b, t):
        return pyunits.convert(
            -1*b.ngcc.st.steam_turbine.power[t],
            pyunits.kW)

    # total gross power
    @m.fs.Expression(m.fs.time)
    def gross_power(b, t):
        return b.gas_turbine_gross_power[t] + b.steam_turbine_gross_power[t]

    # gas turbine auxiliaries
    @m.fs.Expression(m.fs.time)
    def gas_turbine_aux_load(b, t):
        baseline_aux_load = 1020*pyunits.kW
        baseline_gt_power = 477300*pyunits.kW
        return baseline_aux_load*(b.gas_turbine_gross_power[t]/baseline_gt_power)

    # steam turbine auxiliaries
    @m.fs.Expression(m.fs.time)
    def steam_turbine_aux_load(b, t):
        baseline_aux_load = 200*pyunits.kW
        baseline_st_power = 212500*pyunits.kW
        return baseline_aux_load*(b.steam_turbine_gross_power[t]/baseline_st_power)

    # miscellaneous loads
    @m.fs.Expression(m.fs.time)
    def misc_aux_loads(b, t):
        baseline_misc_loads = 572*pyunits.kW
        return baseline_misc_loads*(b.gross_power[t]/baseline_gross_power)

    # transformer losses
    @m.fs.Expression(m.fs.time)
    def transformer_losses(b, t):
        baseline_transformer_losses = 2200*pyunits.kW
        return baseline_transformer_losses*(b.gross_power[t]/baseline_gross_power)

    # total auxiliary load
    @m.fs.Expression(m.fs.time)
    def aux_load(b, t):
        return (b.cooling_water_pump_load[t] +
                b.gas_turbine_aux_load[t] +
                b.cooling_tower_fan_load[t] +
                b.CO2_aux_load[t] +
                b.ground_water_pump_load[t] +
                b.misc_aux_loads[t] +
                b.steam_turbine_aux_load[t] +
                b.transformer_losses[t] +
                pyunits.convert(
                    (b.ngcc.st.cond_pump.work[t] +
                     b.ngcc.hrsg.pump_hp.work[t] +
                     b.ngcc.hrsg.pump_ip.work[t] +
                     b.soec.sweep_turbine.control_volume.work[t] +
                     b.soec.sweep_compressor.control_volume.work[t] +
                     b.soec.sweep_heater.control_volume.heat[t] +
                     b.soec.feed_heater.control_volume.heat[t] +
                     b.soec.water_pump.control_volume.work[t] +
                     b.soec.total_compressor_power[t]),
                    pyunits.kW))

    ###########################################################################
    # OTHER PROCESS PARAMETERS NEEDED FOR COSTING

    # fuel gas flow
    @m.fs.Expression(m.fs.time)
    def fuel_gas_flow(b, t):
        return pyunits.convert(b.ngcc.gt.feed_fuel1.properties[t].flow_mass,
                               pyunits.lb/pyunits.hr)

    # feedwater flow, lb/hr
    @m.fs.Expression(m.fs.time)
    def feedwater_flow(b, t):
        return pyunits.convert(
            b.ngcc.hrsg.econ_hp1.side_1.properties_in[t].flow_mass,
            pyunits.lb/pyunits.hr)

    # HRSG duty, MMBtu/hr
    @m.fs.Expression(m.fs.time)
    def hrsg_duty(b, t):
        return pyunits.convert(
            ((b.ngcc.hrsg.sh_hp4.side_2.properties_in[0.0].flow_mol *
             b.ngcc.hrsg.sh_hp4.side_2.properties_in[0.0].enth_mol) -
             (b.ngcc.hrsg.econ_lp.side_2.properties_out[0.0].flow_mol *
              b.ngcc.hrsg.econ_lp.side_2.properties_in[0.0].enth_mol)),
            pyunits.MBtu/pyunits.hr)

    # gas flow to HRSG, acfm
    @m.fs.Expression(m.fs.time)
    def hrsg_gas_flow(b, t):
        return pyunits.convert(
            b.ngcc.gt.mx3.mixed_state[0].flow_vol,
            pyunits.ft**3/pyunits.min)

    # CO2 capture absorber inlet flow rate, acfm
    @m.fs.Expression(m.fs.time)
    def absorber_gas_flow(b, t):
        return pyunits.convert(
            b.ngcc.hrsg.econ_lp.shell.properties_out[0].flow_vol,
            pyunits.ft**3/pyunits.min)

    # gas flow to stack, acfm
    # scale from baseline based on natural gas flow and adjust for CO2 capture
    @m.fs.Expression(m.fs.time)
    def stack_gas_flow(b, t):
        baseline_stack_flow = 1832949*pyunits.ft**3/pyunits.min
        baseline_natural_gas = 205630*pyunits.lb/pyunits.hr
        return (baseline_stack_flow *
                b.fuel_gas_flow[t]/baseline_natural_gas *
                (1-capture_percentage)/0.1)


def get_ngcc_soec_capital_cost(m):
    m.fs.get_costing(year="2018")

    ###########################################################################
    # FOSSIL ENERGY BASELINE CAPITAL COSTING, NGCC and COMBINED SYSTEMS

    # accounts with feedwater flow to HP section of HRSG, as the process parameter
    FW_accounts = ["3.1", "3.3", "8.4"]

    m.fs.b1 = pyo.Block()

    get_PP_costing(
        m.fs.b1, FW_accounts, m.fs.feedwater_flow[0], "lb/hr", 6)

    # accounts with raw water withdrawal as the process parameter
    water_withdraw_accounts = ["3.2", "3.4", "3.5", "9.5", "14.6"]

    m.fs.b2 = pyo.Block()

    get_PP_costing(
        m.fs.b2, water_withdraw_accounts, m.fs.raw_water_withdrawal[0], "gpm", 6)

    # accounts with fuel gas flowrate as the process parameter
    fuel_flow_accounts = ["3.6", "3.9", "6.1", "6.3", "6.4"]

    m.fs.b3 = pyo.Block()

    get_PP_costing(
        m.fs.b3, fuel_flow_accounts, m.fs.fuel_gas_flow[0], "lb/hr", 6)

    # accounts with process water discharge as the process parameter
    water_discharge_accounts = ["3.7"]

    m.fs.b4 = pyo.Block()

    get_PP_costing(
        m.fs.b4, water_discharge_accounts, m.fs.process_water_discharge[0], "gpm", 6)

    # accounts with CO2 product flow rate as the process parameter
    CO2_accounts_a = ["5.1.a.epri"]

    m.fs.b5a = pyo.Block()

    get_PP_costing(
        m.fs.b5a, CO2_accounts_a, m.fs.CO2_captured[0], "lb/hr", 6)

    CO2_accounts_b = ["5.12"]

    m.fs.b5b = pyo.Block()

    get_PP_costing(
        m.fs.b5b, CO2_accounts_b, m.fs.CO2_captured[0], "lb/hr", 6)

    # accounts with inlet flow rate to CO2 absorber as the process parameter
    CO2_absorber_accounts = ["5.1.b"]

    m.fs.b6 = pyo.Block()

    get_PP_costing(
        m.fs.b6, CO2_absorber_accounts, m.fs.absorber_gas_flow[0], "acfm", 6)

    # accounts with CO2 compressor auxiliary load as the process parameter
    CO2_comp_accounts = ["5.4"]

    m.fs.b7 = pyo.Block()

    get_PP_costing(
        m.fs.b7, CO2_comp_accounts, m.fs.CO2_compressor_load[0], "kW", 6)

    # accounts with CO2 aftercooler duty as the process parameter
    CO2_cooler_accounts = ["5.5"]

    m.fs.b8 = pyo.Block()

    get_PP_costing(
        m.fs.b8, CO2_cooler_accounts, m.fs.CO2_aftercooler_duty[0], "MMBtu/hr", 6)

    # accounts with gas turbine gross power as the process parameter
    gt_power_accounts = ["6.5", "14.1"]

    m.fs.b9 = pyo.Block()

    get_PP_costing(
        m.fs.b9, gt_power_accounts, m.fs.gas_turbine_gross_power[0], "kW", 6)

    # accounts with HRSG duty as the process parameter
    HRSG_duty_accounts = ["7.1", "7.2"]

    m.fs.b10 = pyo.Block()

    get_PP_costing(
        m.fs.b10, HRSG_duty_accounts, m.fs.hrsg_duty[0], "MMBtu/hr", 6)

    # accounts with gas flow to stack as the process parameter
    stack_flow_accounts = ["7.3", "7.4", "7.5"]

    m.fs.b11 = pyo.Block()

    get_PP_costing(
        m.fs.b11, stack_flow_accounts, m.fs.stack_gas_flow[0], "acfm", 6)

    # accounts with flue gas flowrate to HRSG as the process parameter
    hrsg_flow_accounts = ["7.6"]

    m.fs.b12 = pyo.Block()

    get_PP_costing(
        m.fs.b12, hrsg_flow_accounts, m.fs.hrsg_gas_flow[0], "acfm", 6)

    # accounts with steam turbine gross power as the process parameter
    st_power_accounts = ["8.1", "8.2", "8.5", "14.3"]

    m.fs.b13 = pyo.Block()

    get_PP_costing(
        m.fs.b13, st_power_accounts, m.fs.steam_turbine_gross_power[0], "kW", 6)

    # accounts with condenser duty as the process parameter
    condenser_duty_accounts = ["8.3"]

    m.fs.b14 = pyo.Block()

    get_PP_costing(
        m.fs.b14, condenser_duty_accounts, m.fs.condenser_duty[0], "MMBtu/hr", 6)

    # accounts with cooling tower duty as the process parameter
    cooling_tower_accounts = ["9.1"]

    m.fs.b15 = pyo.Block()

    get_PP_costing(
        m.fs.b15, cooling_tower_accounts, m.fs.cooling_water_duty[0], "MMBtu/hr", 6)

    # accounts with circulating water flowrate as the process parameter
    cooling_water_accounts = ["9.2", "9.3", "9.4", "9.6", "9.7", "14.5"]

    m.fs.b16 = pyo.Block()

    get_PP_costing(
        m.fs.b16, cooling_water_accounts, m.fs.cooling_water_flow[0], "gpm", 6)

    # accounts with total plant gross power as the process parameter
    gross_power_accounts = ["11.1", "11.7", "11.9",
                            "13.1", "13.2", "13.3",
                            "14.4", "14.7", "14.8", "14.9", "14.10"]

    m.fs.b17 = pyo.Block()

    get_PP_costing(
        m.fs.b17, gross_power_accounts, m.fs.gross_power[0], "kW", 6)

    # transformer costs
    gross_power_accounts_2 = ["11.8"]

    m.fs.b18 = pyo.Block()

    get_PP_costing(
        m.fs.b18, gross_power_accounts_2, m.fs.gross_power[0], "kW", 6)

    # accounts with auxilliary load as the process parameter
    aux_load_accounts = ["11.2", "11.3", "11.4", "11.5", "11.6",
                         "12.1", "12.2", "12.3", "12.4", "12.5", "12.6",
                         "12.7", "12.8", "12.9"]

    m.fs.b19 = pyo.Block()

    get_PP_costing(m.fs.b19, aux_load_accounts, m.fs.aux_load[0], "kW", 6)

    ###########################################################################
    # SOEC CAPITAL COSTING

    # soec modules
    m.fs.soec.soec_module.costing = pyo.Block()

    @m.fs.soec.soec_module.costing.Expression()
    def total_plant_cost(c):
        extra_installed_area = 1.10  # accounts for cell degradation
        return 60.87 * m.fs.soec.soec_module.number_cells * extra_installed_area / 1e6

    # high temperature soec trim heaters
    for heater in [m.fs.soec.feed_heater, m.fs.soec.sweep_heater]:
        heater.costing = pyo.Block()

        @heater.costing.Expression()
        def total_plant_cost(b):
            U = 56 * pyo.units.W / pyo.units.m ** 2 / pyo.units.K
            DT = 20 * pyo.units.K
            area = heater.heat_duty[0] / U / DT
            # Add factor of two to pathways cost to account for corrosion-resistant materials for trim heaters
            return 2*81.88*pyo.units.convert(area, pyo.units.ft**2) / 1e6

    # water heaters - U-tube HXs
    # costed with IDAES generic heat exchanger correlation
    m.fs.soec.water_heater01.get_costing(hx_type="U-tube")
    add_total_plant_cost(m.fs.soec.water_heater01)

    m.fs.soec.water_heater02.get_costing(hx_type="U-tube")
    add_total_plant_cost(m.fs.soec.water_heater02)

    # sweep hx and feed hx
    # costed with price/ft^2 from NGFC Pathways study
    for hx in [m.fs.soec.sweep_hx, m.fs.soec.feed_hx01]:
        hx.costing = pyo.Block()

        @hx.costing.Expression()
        def total_plant_cost(c):
            return 81.88*pyo.units.convert(hx.area, pyo.units.ft**2) / 1e6

    # sweep compressor
    m.fs.soec.sweep_compressor.costing = pyo.Block()

    # air is coming in at standard conditions
    sweep_compressor_scfm = pyo.units.convert(
        m.fs.soec.sweep_compressor.control_volume.properties_in[0].flow_vol,
        pyo.units.ft**3/pyo.units.min)

    n_sections = 11

    @m.fs.soec.sweep_compressor.costing.Expression()
    def total_plant_cost(c):  # in 2018 $
        return n_sections * 10.9995 * pyo.exp(
            0.4692 +
            0.1203*pyo.log(sweep_compressor_scfm/1000/n_sections) +
            0.0931*pyo.log(sweep_compressor_scfm/1000/n_sections)**2) / 1e3

    # sweep turbine
    # costed with IDAES generic turbine correlation
    m.fs.soec.sweep_turbine.get_costing()
    add_total_plant_cost(m.fs.soec.sweep_turbine)

    # H2 compressor
    m.fs.soec.cmp01.costing = pyo.Block()

    h2_comp_process_param = pyo.units.convert(
        m.fs.soec.cmp01.control_volume.properties_out[0].flow_mass,
        pyo.units.lb/pyo.units.hr)

    @m.fs.soec.cmp01.costing.Expression()
    def total_plant_cost(c):
        # This is for 4 stages.  The original was for 2 hence the multiply by 2
        ref_cost = 11.408  # MM$ 2018
        ref_param = 44369*pyo.units.lb/pyo.units.hr
        alpha = 0.7
        return 2*ref_cost*(h2_comp_process_param/ref_param)**alpha

    # build cost constraints
    get_total_TPC(m.fs)

    # costing initialization
    costing_initialization(m.fs)

    initialize_unit_costing(m.fs.soec.water_heater01.costing)
    initialize_unit_costing(m.fs.soec.water_heater02.costing)
    initialize_unit_costing(m.fs.soec.sweep_turbine.costing)

    calculate_variable_from_constraint(
        m.fs.costing.total_TPC,
        m.fs.costing.total_TPC_eq)


def display_capital_costs(m):

    @m.fs.costing.Expression()
    def feedwater_and_miscellaneous_bop_cost(b):
        return (m.fs.b1.costing.total_plant_cost["3.1"] +
                m.fs.b2.costing.total_plant_cost["3.2"] +
                m.fs.b1.costing.total_plant_cost["3.3"] +
                m.fs.b2.costing.total_plant_cost["3.4"] +
                m.fs.b2.costing.total_plant_cost["3.5"] +
                m.fs.b3.costing.total_plant_cost["3.6"] +
                m.fs.b4.costing.total_plant_cost["3.7"] +
                m.fs.b3.costing.total_plant_cost["3.9"])

    @m.fs.costing.Expression()
    def carbon_capture_system_cost(b):
        return (m.fs.b5a.costing.total_plant_cost["5.1.a.epri"] +
                m.fs.b5b.costing.total_plant_cost["5.12"] +
                m.fs.b6.costing.total_plant_cost["5.1.b"] +
                m.fs.b7.costing.total_plant_cost["5.4"] +
                m.fs.b8.costing.total_plant_cost["5.5"])

    @m.fs.costing.Expression()
    def combustion_turbine_cost(b):
        return (m.fs.b3.costing.total_plant_cost["6.1"] +
                m.fs.b3.costing.total_plant_cost["6.3"] +
                m.fs.b3.costing.total_plant_cost["6.4"] +
                m.fs.b9.costing.total_plant_cost["6.5"])

    @m.fs.costing.Expression()
    def hrsg_ductwork_and_stack_cost(b):
        return (m.fs.b10.costing.total_plant_cost["7.1"] +
                m.fs.b10.costing.total_plant_cost["7.2"] +
                m.fs.b11.costing.total_plant_cost["7.3"] +
                m.fs.b11.costing.total_plant_cost["7.4"] +
                m.fs.b11.costing.total_plant_cost["7.5"] +
                m.fs.b12.costing.total_plant_cost["7.6"])

    @m.fs.costing.Expression()
    def steam_turbine_cost(b):
        return (m.fs.b13.costing.total_plant_cost["8.1"] +
                m.fs.b13.costing.total_plant_cost["8.2"] +
                m.fs.b14.costing.total_plant_cost["8.3"] +
                m.fs.b1.costing.total_plant_cost["8.4"] +
                m.fs.b13.costing.total_plant_cost["8.5"])

    @m.fs.costing.Expression()
    def cooling_water_system_cost(b):
        return (m.fs.b15.costing.total_plant_cost["9.1"] +
                m.fs.b16.costing.total_plant_cost["9.2"] +
                m.fs.b16.costing.total_plant_cost["9.3"] +
                m.fs.b16.costing.total_plant_cost["9.4"] +
                m.fs.b2.costing.total_plant_cost["9.5"] +
                m.fs.b16.costing.total_plant_cost["9.6"] +
                m.fs.b16.costing.total_plant_cost["9.7"])

    @m.fs.costing.Expression()
    def accessory_electric_plant_cost(b):
        return (m.fs.b17.costing.total_plant_cost["11.1"] +
                m.fs.b19.costing.total_plant_cost["11.2"] +
                m.fs.b19.costing.total_plant_cost["11.3"] +
                m.fs.b19.costing.total_plant_cost["11.4"] +
                m.fs.b19.costing.total_plant_cost["11.5"] +
                m.fs.b19.costing.total_plant_cost["11.6"] +
                m.fs.b17.costing.total_plant_cost["11.7"] +
                m.fs.b18.costing.total_plant_cost["11.8"] +
                m.fs.b17.costing.total_plant_cost["11.9"])

    @m.fs.costing.Expression()
    def instrumentation_and_control_cost(b):
        return (m.fs.b19.costing.total_plant_cost["12.1"] +
                m.fs.b19.costing.total_plant_cost["12.2"] +
                m.fs.b19.costing.total_plant_cost["12.3"] +
                m.fs.b19.costing.total_plant_cost["12.4"] +
                m.fs.b19.costing.total_plant_cost["12.5"] +
                m.fs.b19.costing.total_plant_cost["12.6"] +
                m.fs.b19.costing.total_plant_cost["12.7"] +
                m.fs.b19.costing.total_plant_cost["12.8"] +
                m.fs.b19.costing.total_plant_cost["12.9"])

    @m.fs.costing.Expression()
    def improvements_to_site_cost(b):
        return (m.fs.b17.costing.total_plant_cost["13.1"] +
                m.fs.b17.costing.total_plant_cost["13.2"] +
                m.fs.b17.costing.total_plant_cost["13.3"])

    @m.fs.costing.Expression()
    def buildings_and_structures_cost(b):
        return (m.fs.b9.costing.total_plant_cost["14.1"] +
                m.fs.b13.costing.total_plant_cost["14.3"] +
                m.fs.b17.costing.total_plant_cost["14.4"] +
                m.fs.b16.costing.total_plant_cost["14.5"] +
                m.fs.b2.costing.total_plant_cost["14.6"] +
                m.fs.b17.costing.total_plant_cost["14.7"] +
                m.fs.b17.costing.total_plant_cost["14.8"] +
                m.fs.b17.costing.total_plant_cost["14.9"] +
                m.fs.b17.costing.total_plant_cost["14.10"])

    @m.fs.costing.Expression()
    def soec_module_cost(b):
        return m.fs.soec.soec_module.costing.total_plant_cost

    @m.fs.costing.Expression()
    def trim_heater_cost(b):
        return (m.fs.soec.feed_heater.costing.total_plant_cost +
                m.fs.soec.sweep_heater.costing.total_plant_cost)

    @m.fs.costing.Expression()
    def hydrogen_compressor_cost(b):
        return (m.fs.soec.cmp01.costing.total_plant_cost)

    @m.fs.costing.Expression()
    def sweep_compressor_and_turbine_cost(b):
        return (m.fs.soec.sweep_compressor.costing.total_plant_cost +
                m.fs.soec.sweep_turbine.costing.total_plant_cost)

    @m.fs.costing.Expression()
    def heat_exchanger_costs(b):
        return (m.fs.soec.sweep_hx.costing.total_plant_cost +
                m.fs.soec.feed_hx01.costing.total_plant_cost +
                m.fs.soec.water_heater01.costing.total_plant_cost +
                m.fs.soec.water_heater02.costing.total_plant_cost)

    print("Total plant cost: ${:.2f}M".format(
        pyo.value(m.fs.costing.total_TPC)))
    print("Feedwater & misc. BOP cost: ${:.2f}M".format(
        pyo.value(m.fs.costing.feedwater_and_miscellaneous_bop_cost)))
    print("Carbon capture system cost: ${:.2f}M".format(
        pyo.value(m.fs.costing.carbon_capture_system_cost)))
    print("Combustion turbine cost: ${:.2f}M".format(
        pyo.value(m.fs.costing.combustion_turbine_cost)))
    print("HRSG, ductwork, and stack cost: ${:.2f}M".format(
        pyo.value(m.fs.costing.hrsg_ductwork_and_stack_cost)))
    print("Steam turbine cost: ${:.2f}M".format(
        pyo.value(m.fs.costing.steam_turbine_cost)))
    print("Cooling water system cost: ${:.2f}M".format(
        pyo.value(m.fs.costing.cooling_water_system_cost)))
    print("Accessory plant electric cost: ${:.2f}M".format(
        pyo.value(m.fs.costing.accessory_electric_plant_cost)))
    print("Instrumentation and control cost: ${:.2f}M".format(
        pyo.value(m.fs.costing.instrumentation_and_control_cost)))
    print("Improvements to site cost: ${:.2f}M".format(
        pyo.value(m.fs.costing.improvements_to_site_cost)))
    print("Buildings and structures cost: ${:.2f}M".format(
        pyo.value(m.fs.costing.buildings_and_structures_cost)))
    print("SOEC module cost: ${:.2f}M".format(
        pyo.value(m.fs.costing.soec_module_cost)))
    print("Trim heater cost: ${:.2f}M".format(
        pyo.value(m.fs.costing.trim_heater_cost)))
    print("Hydrogen compressor cost: ${:.2f}M".format(
        pyo.value(m.fs.costing.hydrogen_compressor_cost)))
    print("Sweep turbomachinery cost: ${:.2f}M".format(
        pyo.value(m.fs.costing.sweep_compressor_and_turbine_cost)))
    print("Heat exchanger cost: ${:.2f}M".format(
        pyo.value(m.fs.costing.heat_exchanger_costs)))


# if __name__ == "__main__":
#     # set up solver
#     use_idaes_solver_configuration_defaults()
#     idaes.cfg.ipopt.options.nlp_scaling_method = "user-scaling"
#     solver = pyo.SolverFactory("ipopt")

#     # create the ngcc model
#     m = pyo.ConcreteModel()
#     m.fs = ngcc_soec.NgccSoecFlowsheet(default={"dynamic": False})
#     iscale.calculate_scaling_factors(m)
#     m.fs.initialize(
#         load_from="ngcc_soec_init.json.gz",
#         save_to="ngcc_soec_init.json.gz")
#     # res = solver.solve(m, tee=True)

#     # add capital costing
#     add_results_for_costing(m)
#     get_ngcc_soec_capital_cost(m)
#     # solver.solve(m, tee=True)

#     display_capital_costs(m)
