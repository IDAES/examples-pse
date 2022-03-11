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

__author__ = "John Eslick"

import os
import pyomo.environ as pyo
from pyomo.network import Arc
import idaes.generic_models.unit_models as um # um = unit models
from idaes.core import FlowsheetBlockData, declare_process_block_class
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
import gas_turbine
import hrsg
import steam_turbine
import idaes.core.util as iutil
from idaes.core.util.initialization import propagate_state
from idaes.power_generation.costing.power_plant_costing import (
    get_fixed_OM_costs,
    get_variable_OM_costs,
    initialize_fixed_OM_costs,
    initialize_variable_OM_costs,
)

def _build_ngcc_OM_costs(m):
    """Cost evaluation using IDAES Costing Framework for EMRE NGCC Model
    Base Case : B31B - NETL Baseline Report Rev 4
    Author: A. Deshpande, Alex Noring, M. Zamarripa
    """
    get_fixed_OM_costs(
        m, 650, operators_per_shift=5, tech=6, fixed_TPC=555.35  # MW
    )  # MW
    initialize_fixed_OM_costs(m)

    m.net_power_cost = pyo.Var(m.time, initialize=650, units=pyo.units.MW)
    m.net_power_cost.fix()
    m.natural_gas = pyo.Var(
        m.time, initialize={0: 7239}, units=pyo.units.MBtu / pyo.units.day
    )
    m.natural_gas.fix(110955)
    m.water_use = pyo.Var(
        m.time, initialize={0: 2090 * 1000}, units=pyo.units.gallon / pyo.units.day
    )
    m.water_use.fix()
    #
    get_variable_OM_costs(
        m,
        m.net_power_cost,
        ["natural gas", "water"],
        [m.natural_gas, m.water_use],
    )
    initialize_variable_OM_costs(m)

@declare_process_block_class(
    "NgccFlowsheet",
    doc=(
        "The NGCC flowsheet is base on NETL report 'Cost and Performance "
        "Baseline for Fossil Energy Plants Volume 1: Bituminous Coal and Natural "
        "Gas to Electricity.' Sept 2019, Case B31B."
    ),
)
class NgccFlowsheetData(FlowsheetBlockData):
    def build(self):
        super().build()
        self._add_flowsheets()
        self._add_units()
        self._add_arcs()
        self._add_constraints()
        self._scaling_guess()
        self._add_tags()

    def _add_tags(self):
        tag_group = iutil.ModelTagGroup()
        self.tags_output = tag_group
        tag_group["st_power"] = iutil.ModelTag(
            doc=f"Steam turbine electric power output",
            expr=-self.st.steam_turbine.power[0],
            format_string="{:.2f}",
            display_units=pyo.units.MW,
        )
        tag_group["gt_power"] = iutil.ModelTag(
            doc=f"Gas turbine electric power output",
            expr=-self.gt.gt_power[0],
            format_string="{:.2f}",
            display_units=pyo.units.MW,
        )
        tag_group["gross_power"] = iutil.ModelTag(
            doc=f"Gross electric power output",
            expr=-self.gross_power[0],
            format_string="{:.2f}",
            display_units=pyo.units.MW,
        )
        tag_group["net_power"] = iutil.ModelTag(
            doc=f"Net electric power output",
            expr=self.net_power_mw[0],
            format_string="{:.2f}",
            display_units=pyo.units.MW,
        )
        tag_group["lhv_efficiency"] = iutil.ModelTag(
            doc=f"Overall LHV efficiency",
            expr=100*self.lhv_efficiency[0],
            format_string="{:.2f}",
            display_units="%",
        )
        tag_group["combustor_temperature"] = iutil.ModelTag(
            doc=f"Overall LHV efficiency",
            expr=self.gt.cmb1.control_volume.properties_out[0].temperature,
            format_string="{:.2f}",
            display_units=pyo.units.K,
        )
        tag_group["fuel_flow"] = iutil.ModelTag(
            doc=f"Fuel mass flow",
            expr=self.gt.feed_fuel1.properties[0].flow_mass,
            format_string="{:.3f}",
            display_units=pyo.units.kg / pyo.units.s,
        )
        st = self.st.steam_turbine
        tag_group["st_throttle_delta_pressure"] = iutil.ModelTag(
            doc=f"Pressure change in the steam turbine throttle valve",
            expr=st.throttle_valve[1].deltaP[0],
            format_string="{:.2f}",
            display_units=pyo.units.bar,
        )
        tag_group["st_condenser_pressure"] = iutil.ModelTag(
            doc=f"Steam turbine comdenser pressure",
            expr=st.outlet_stage.control_volume.properties_out[0].pressure,
            format_string="{:.2f}",
            display_units=pyo.units.bar,
        )
        tag_group["st_condenser_pressure"] = iutil.ModelTag(
            doc=f"Steam turbine condenser pressure",
            expr=st.outlet_stage.control_volume.properties_out[0].pressure,
            format_string="{:.2f}",
            display_units=pyo.units.bar,
        )
        tag_group["st_throttle_inlet_temperature"] = iutil.ModelTag(
            doc=f"Steam turbine throttle valve outlet temperature",
            expr=st.throttle_valve[1].control_volume.properties_in[0].temperature,
            format_string="{:.2f}",
            display_units=pyo.units.K,
        )
        tag_group["st_throttle_outlet_temperature"] = iutil.ModelTag(
            doc=f"Steam turbine throttle valve outlet temperature",
            expr=st.throttle_valve[1].control_volume.properties_out[0].temperature,
            format_string="{:.2f}",
            display_units=pyo.units.K,
        )


    def _add_flowsheets(self):
        self.gt = gas_turbine.GasTurbineFlowsheet(
            default={
                "dynamic":self.config.dynamic,
                "time":self.time,
                "time_units": self.config.time_units
            }
        )
        self.hrsg = hrsg.HrsgFlowsheet(
            default={
                "dynamic":self.config.dynamic,
                "time":self.time,
                "time_units": self.config.time_units
            }
        )
        self.st = steam_turbine.SteamTurbineFlowsheet(
            default={
                "dynamic":self.config.dynamic,
                "time":self.time,
                "time_units": self.config.time_units
            }
        )

    def _add_units(self):
        self.fg_translate = um.Translator(
            doc="Translate from generic to flue gas prop and drop Ar.",
            default={
                "inlet_property_package": self.gt.flue_prop_params,
                "outlet_property_package": self.hrsg.prop_gas
            }
        )

    def _add_arcs(self):
        self.g08a = Arc(
            source=self.gt.mx3.outlet,
            destination=self.fg_translate.inlet,
        )
        self.g08b = Arc(
            source=self.fg_translate.outlet,
            destination=self.hrsg.sh_hp4.shell_inlet,
        )
        self.t05a = Arc(
            source=self.hrsg.sh_lp.tube_outlet,
            destination=self.st.steam_turbine_lp_mix.hrsg
        )
        self.t02a = Arc(
            source=self.st.steam_turbine.hp_stages[7].outlet,
            destination=self.hrsg.splitter_ip2.inlet
        )
        self.t03a = Arc(
            source=self.hrsg.sh_ip3.tube_outlet,
            destination=self.st.steam_turbine.ip_stages[1].inlet
        )
        self.t11a = Arc(
            source=self.st.return_mix.outlet,
            destination=self.hrsg.econ_lp.tube_inlet
        )
        self.t01a = Arc(
            source=self.hrsg.sh_hp4.tube_outlet,
            destination=self.st.steam_turbine.inlet_split.inlet
        )
        self.st01a = Arc(
            source=self.hrsg.splitter_ip1.toNGPH,
            destination=self.gt.ng_preheat.shell_inlet
        )
        self.st02a = Arc(
            source=self.gt.ng_preheat.shell_outlet,
            destination=self.hrsg.mixer1.Preheater
        )
        pyo.TransformationFactory("network.expand_arcs").apply_to(self)

    def _add_constraints(self):
        self.net_power_mw = pyo.Var(
            self.time,
            initialize=600,
            units=pyo.units.MW
        )
        # default to 90% to compare to baseline
        self.cap_fraction = pyo.Var(
            initialize=0.90,
            units=pyo.units.dimensionless
        )
        # Specific Reboiler Duty (SRD)
        # From baseline with 90% capture 2.7 MJ/kg
        # For PZ and 97% capture 3.6 MJ/kg
        # For PZAS as 97% capture 2.4 MJ/kg
        self.cap_specific_reboiler_duty = pyo.Var(
            initialize=2.7e6,
            units=pyo.units.J/pyo.units.kg
        )
        self.cap_addtional_co2 = pyo.Var(
            self.config.time,
            initialize=0.0,
            units=pyo.units.kg/pyo.units.s
        )
        self.cap_specific_compression_power = pyo.Var(
            initialize=0.2748e6,
            units=pyo.units.J/pyo.units.kg
        )
        self.cap_additional_reboiler_duty = pyo.Var(
            self.config.time,
            initialize=0.0,
            units=pyo.units.W
        )
        self.fuel_lhv = pyo.Var(
            initialize=47.2e6,
            units=pyo.units.J/pyo.units.kg
        )
        self.fuel_hhv = pyo.Var(
            initialize=52.3e6,
            units=pyo.units.J/pyo.units.kg
        )
        self.lp_steam_temperature = pyo.Var(
            self.config.time,
            initialize=554.0,
            units=pyo.units.K
        )

        @self.fg_translate.Constraint(self.time, self.hrsg.prop_gas.component_list)
        def mol_frac_eqn(b, t, i):
            return (
                b.outlet.flow_mol_comp[t, i] ==
                b.inlet.flow_mol[t] * b.inlet.mole_frac_comp[t, i]
            )
        @self.fg_translate.Constraint(self.time)
        def temperature_eqn(b, t):
            return b.outlet.temperature[t] == b.inlet.temperature[t]
        @self.fg_translate.Constraint(self.time)
        def pressure_eqn(b, t):
            return b.outlet.pressure[t] == b.inlet.pressure[t]

        @self.Expression(self.config.time)
        def gross_power(b, t):
            return b.gt.gt_power[t] + b.st.steam_turbine.power[t]

        @self.Expression(self.config.time)
        def aux_cooling_fans(b, t):
            return 1e3 * (2370) * pyo.units.MW

        @self.Expression(self.config.time)
        def aux_cooling_pumps(b, t):
            return 1e3 * (4580) * pyo.units.MW

        # Aux power expressions
        @self.Expression(self.config.time)
        def aux_combustion(b, t):
            return 1e3 * 1020.0 * pyo.units.MW

        @self.Expression(self.config.time)
        def aux_capture(b, t): #scale to flue gas flow
            return (
                10600000*pyo.units.W/62.1/pyo.units.kg*pyo.units.s * b.cap_fraction *
                b.gt.gts2.control_volume.properties_out[t].flow_mol_comp["CO2"] *
                0.04401 * pyo.units.kg / pyo.units.mol
            )

        @self.Expression(self.config.time)
        def aux_compression(b, t):
            return (
                b.cap_specific_compression_power * b.cap_fraction *
                b.gt.gts2.control_volume.properties_out[t].flow_mol_comp["CO2"] *
                0.04401 * pyo.units.kg / pyo.units.mol
            )

        @self.Expression(self.config.time)
        def aux_transformer(b, t): # scale to gross power
            return - 1e3 * 2200 * b.gross_power[t] / 687.0e6

        @self.Expression(self.config.time)
        def aux_misc(b, t):
            # turbine, scr, gw pumps, bop
            return 1e3 * 1202.0 * pyo.units.MW

        @self.Expression(self.config.time)
        def net_power(b, t):
            return (
                b.gt.gt_power[t] +
                b.st.steam_turbine.power[t] +
                b.st.cond_pump.work[t] +
                b.hrsg.pump_hp.work[t] +
                b.hrsg.pump_ip.work[t] +
                b.aux_cooling_pumps[t] +
                b.aux_cooling_fans[t] +
                b.aux_combustion[t] +
                b.aux_capture[t] +
                b.aux_compression[t] +
                b.aux_transformer[t] +
                b.aux_misc[t]
            )

        @self.Expression(self.config.time)
        def fuel_thermal_in_mbtu(b, t):
            return pyo.units.convert(
                b.fuel_lhv * b.gt.inject1.gas_state[t].flow_mass,
                pyo.units.MBtu/pyo.units.hr
            )

        @self.Expression(self.config.time)
        def lhv_efficiency(b, t):
            return -b.net_power[t] / b.gt.inject1.gas_state[t].flow_mass / b.fuel_lhv

        @self.Expression(self.config.time)
        def reboiler_duty_expr(b, t): #scale to flue gas flow
            return (
                -b.cap_specific_reboiler_duty * b.cap_fraction *
                (
                    b.gt.gts2.control_volume.properties_out[0].flow_mol_comp["CO2"] *
                    0.04401 * pyo.units.kg/pyo.units.mol + b.cap_addtional_co2[t]
                ) +
                b.cap_additional_reboiler_duty[t]
            )

        @self.Constraint(self.config.time)
        def reboiler_duty_eqn(b, t):
            return b.reboiler_duty_expr[t] == b.st.reboiler.control_volume.heat[t]

        @self.Constraint(self.config.time)
        def net_power_constraint(b, t):
            return b.net_power_mw[t]/100.0 == -b.net_power[t]/1e6/100.0

        @self.Constraint(self.config.time)
        def lp_steam_temperature_eqn(b, t):
            return (b.lp_steam_temperature[t] ==
                b.hrsg.sh_lp.tube.properties_out[t].temperature)

    def _scaling_guess(self):
        for i, c in self.fg_translate.mol_frac_eqn.items():
            iscale.constraint_scaling_transform(c, 1e-2)
        for t, c in self.fg_translate.temperature_eqn.items():
            iscale.constraint_scaling_transform(c, 1e-2)
        for t, c in self.fg_translate.pressure_eqn.items():
            iscale.constraint_scaling_transform(c, 1e-5)
        for t, c in self.reboiler_duty_eqn.items():
            iscale.constraint_scaling_transform(c, 1e-7)
        iscale.set_scaling_factor(
            self.st.steam_turbine.throttle_valve[1].control_volume.deltaP, 1e-6)
        iscale.set_scaling_factor(
            self.st.steam_turbine.hp_stages[1].control_volume.deltaP, 1e-6)
        iscale.set_scaling_factor(self.st.main_condenser.tube.heat, 1e-8)
        iscale.set_scaling_factor(self.st.main_condenser.shell.heat, 1e-8)

    def initialize(
        self,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
        load_from="ngcc_init.json.gz",
        save_to="ngcc_init.json.gz",
    ):

        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="flowsheet")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="flowsheet")
        solver_obj = iutil.get_solver(solver, optarg)

        if load_from is not None and os.path.exists(load_from):
            init_log.info(f"NGCC load initial from {load_from}")
            # here suffix=False avoids loading scaling factors
            iutil.from_json(
                self, fname=load_from, wts=iutil.StoreSpec(suffix=False))
        else:
            self.cap_addtional_co2.fix()
            self.cap_fraction.fix()
            self.cap_specific_reboiler_duty.fix()
            self.cap_specific_compression_power.fix()
            self.cap_additional_reboiler_duty.fix()
            self.fuel_lhv.fix()
            self.fuel_hhv.fix()

            self.gt.initialize(
                load_from="data_init/gas_turbine_init.json.gz",
                save_to="data_init/gas_turbine_init.json.gz",
            )
            propagate_state(self.g08a)
            self.fg_translate.initialize()
            propagate_state(self.g08b, overwrite_fixed=True)
            self.hrsg.initialize(
                load_from="data_init/hrsg_init.json.gz",
                save_to="data_init/hrsg_init.json.gz",
            )
            self.hrsg.sh_hp4.shell_inlet.unfix()
            propagate_state(self.t05a, overwrite_fixed=True)
            self.st.initialize(
                load_from="data_init/steam_turbine_init.json.gz",
                save_to="data_init/steam_turbine_init.json.gz",
            )

            init_log.info(f"Open tears")
            self.st02a_expanded.deactivate() # steam from ng preheat
            self.st01a_expanded.deactivate() # steam to ng preheat
            self.t01a_expanded.deactivate() # main steam to turbine
            self.t02a_expanded.deactivate() # cold reheat
            self.t03a_expanded.deactivate() # hot reheat
            self.st.steam_turbine_lp_mix.hrsg.unfix()
            self.t11a_expanded.deactivate()
            #self.hrsg.econ_lp.tube_inlet.unfix()
            self.st.steam_turbine_lp_split.reboiler.flow_mol.unfix()
            solver_obj.solve(self, tee=True)

            init_log.info(f"HRSG flow constraint active")
            self.hrsg.evap_hp.hp_sat_vap_eqn.activate()
            self.hrsg.econ_lp.tube_inlet.flow_mol.unfix()
            solver_obj.solve(self, tee=True)

            # hook in preheater
            init_log.info(f"Connect preheater and reheater")
            self.st02a_expanded.activate() # steam from ng preheat
            self.st01a_expanded.activate() # steam to ng preheat
            self.gt.ng_preheat.shell_inlet.unfix()
            self.hrsg.mixer1.Preheater.unfix()
            #propagate_state(self.t01a, overwrite_fixed=True)
            propagate_state(self.t02a, overwrite_fixed=True)
            propagate_state(self.t03a, overwrite_fixed=True)
            self.hrsg.splitter_ip2.inlet.unfix()
            self.st.t02_dummy.deactivate()
            self.st.t03_dummy.deactivate()
            self.st.dummy_reheat.deactivate()
            self.t02a_expanded.activate() # hot reheat
            self.t03a_expanded.activate() # cold reheat
            solver_obj.solve(self, tee=True)

            # hook main steam
            init_log.info(f"Finish turbine sizing/connect main steam")
            self.t01a_expanded.activate()
            #propagate_state(self.t01a, overwrite_fixed=True)
            self.st.steam_turbine.outlet_stage.flow_coeff.unfix()
            self.st.steam_turbine.inlet_split.inlet.unfix()
            solver_obj.solve(self, tee=True)

            init_log.info(f"Fix flow coefficent and free throttle")
            self.st.steam_turbine.throttle_valve[1].pressure_flow_equation.deactivate()
            self.st.steam_turbine.outlet_stage.flow_coeff.fix()
            solver_obj.solve(self, tee=True)

            init_log.info(f"Connect feedwater")
            self.t11a_expanded.activate()
            self.st.hotwell.makeup.flow_mol.unfix()
            self.hrsg.econ_lp.tube_inlet.unfix()
            solver_obj.solve(self, tee=True)

            init_log.info("Set estimated parameters")
            self.hrsg.evap_hp.area.set_value(15000)
            self.hrsg.evap_ip.area.set_value(11000)
            self.hrsg.evap_lp.area.set_value(14000)
            self.hrsg.split_fg_lp.split_fraction[:, "toLP_SH"].set_value(0.5)
            self.hrsg.sh_hp4.fcorrection_htc.set_value(0.65939)
            self.hrsg.sh_ip3.fcorrection_htc.set_value(1.372068)
            self.hrsg.econ_hp1.fcorrection_htc.set_value(7.863405)
            self.hrsg.econ_ip1.fcorrection_htc.set_value(0.6)
            self.hrsg.econ_lp.fcorrection_htc.set_value(20)
            self.gt.cmp1.efficiency_isentropic[:].set_value(0.840187)
            self.st.steam_turbine.outlet_stage.flow_coeff.set_value(0.1099345)
            self.st.main_condenser.tube_inlet.enth_mol.fix(1260)
            self.st.main_condenser.area.fix(1000)
            self.hrsg.evap_hp_valve.deltaP.fix(-7.0e6)

            init_log.info("Set net power to 646 MW")
            self.net_power_mw.fix(646)
            self.gt.gt_power.unfix()
            solver_obj.solve(self, tee=True)

            init_log.info("Unfix GT exhaust pressure and fix stack pressure")
            self.gt.exhaust_1.pressure.unfix()
            self.hrsg.econ_lp.shell.properties_out[:].pressure.fix(1.01e5)
            solver_obj.solve(self, tee=True)

            if save_to is not None:
                iutil.to_json(self)
                if save_to is not None:
                    iutil.to_json(self, fname=save_to)
                    init_log.info(f"Initialization saved to {save_to}")

    def check_scaling(self):
        jac, nlp = iscale.get_jacobian(self, scaled=True)
        print("Extreme Jacobian entries:")
        for i in iscale.extreme_jacobian_entries(jac=jac, nlp=nlp, large=100):
            print(f"    {i[0]:.2e}, [{i[1]}, {i[2]}]")
        print("Badly scaled variables:")
        for v, sv in iscale.badly_scaled_var_generator(m, large=1e2, small=1e-2, zero=1e-12):
            print(f"    {v} -- {sv} -- {iscale.get_scaling_factor(v)}")
        print(f"Jacobian Condition Number: {iscale.jacobian_cond(jac=jac):.2e}")

    def add_costing(self):
        _build_ngcc_OM_costs(self)
        self.tags_output["variable_om_cost_per_time"] = iutil.ModelTag(
            doc=f"Variable OM cost currency per time",
            expr=self.costing.total_variable_OM_cost[0]*(-self.net_power[0]),
            format_string="{:.2f}",
            display_units=pyo.units.USD/pyo.units.hour,
        )
        self.tags_output["variable_om_cost_per_energy"] = iutil.ModelTag(
            doc=f"Variable OM cost currency per time",
            expr=self.costing.total_variable_OM_cost[0],
            format_string="{:.2f}",
            display_units=pyo.units.USD/pyo.units.MWh,
        )
        @self.Constraint(self.config.time)
        def cost_eq1(c, t):
            return self.natural_gas[t] == pyo.units.convert(
                self.gt.inject1.gas_state[t].flow_mass*self.fuel_hhv,
                pyo.units.MBtu/pyo.units.day)

        self.natural_gas.unfix()

        @self.Constraint(self.config.time)
        def cost_eq2(c, t):
            return self.net_power_cost[t] == self.net_power_mw[t]

        self.net_power_cost.unfix()

        for t in self.cost_eq1:
            iscale.constraint_scaling_transform(self.cost_eq1[t], 1e-4)
        for t in self.cost_eq2:
            iscale.constraint_scaling_transform(self.cost_eq2[t], 1e-2)
