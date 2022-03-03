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

    def _scaling_guess(self):
        for i, c in self.fg_translate.mol_frac_eqn.items():
            iscale.constraint_scaling_transform(c, 1e-2)
        for t, c in self.fg_translate.temperature_eqn.items():
            iscale.constraint_scaling_transform(c, 1e-2)
        for t, c in self.fg_translate.pressure_eqn.items():
            iscale.constraint_scaling_transform(c, 1e-5)
        for t, c in self.reboiler_duty_eqn.items():
            iscale.constraint_scaling_transform(c, 1e-7)

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

        if load_from is not None:
            if os.path.exists(load_from):
                init_log.info_high(f"NGCC load initial from {load_from}")
                # here suffix=False avoids loading scaling factors
                iutil.from_json(
                    self, fname=load_from, wts=iutil.StoreSpec(suffix=False)
                )
                return

        self.cap_addtional_co2.fix()
        self.cap_fraction.fix()
        self.cap_specific_reboiler_duty.fix()
        self.cap_specific_compression_power.fix()
        self.cap_additional_reboiler_duty.fix()
        self.fuel_lhv.fix()
        self.fuel_hhv.fix()

        solver_obj = iutil.get_solver(solver, optarg)

        self.gt.initialize()
        propagate_state(self.g08a)
        self.fg_translate.initialize()
        propagate_state(self.g08b, overwrite_fixed=True)
        self.hrsg.initialize()
        self.hrsg.sh_hp4.shell_inlet.unfix()
        propagate_state(self.t05a, overwrite_fixed=True)
        self.st.initialize()

        # solver sheet with some disconnected streams
        #self.st02a_expanded.deactivate() # steam from ng preheat
        #self.st01a_expanded.deactivate() # steam to ng preheat
        self.t01a_expanded.deactivate() # main steam to turbine
        #self.t02a_expanded.deactivate() # cold reheat
        #self.t03a_expanded.deactivate() # hot reheat
        self.st.steam_turbine_lp_mix.hrsg.unfix()
        self.hrsg.econ_lp.tube_inlet.unfix()
        self.st.steam_turbine_lp_split.reboiler.flow_mol.unfix()
        #solver_obj.solve(self, tee=True)

        self.t01a.source.display()

        # hook in preheater
        print("preheater and reheater")
        #self.st02a_expanded.activate() # steam from ng preheat
        #self.st01a_expanded.activate() # steam to ng preheat
        self.gt.ng_preheat.shell_inlet.unfix()
        self.hrsg.mixer1.Preheater.unfix()
        #propagate_state(self.t01a, overwrite_fixed=True)
        propagate_state(self.t02a, overwrite_fixed=True)
        propagate_state(self.t03a, overwrite_fixed=True)
        self.hrsg.splitter_ip2.inlet.unfix()
        self.st.t02_dummy.deactivate()
        self.st.t03_dummy.deactivate()
        self.st.dummy_reheat.deactivate()
        #self.t02a_expanded.activate() # hot reheat
        #self.t03a_expanded.activate() # cold reheat
        solver_obj.solve(self, tee=True)

        #self.t01a.source.display()

        # hook main steam
        print("main steam")
        self.t01a_expanded.activate()
        self.st.hotwell.makeup.flow_mol[:].unfix()
        self.st.steam_turbine.inlet_split.inlet.unfix()
        #propagate_state(self.t01a)
        #self.hrsg.pump_hp.outlet.pressure[0].unfix()
        #self.st.steam_turbine.inlet_split.inlet.flow_mol.fix()
        solver_obj.solve(self, tee=True)

        print("fix steam temperature")
        self.st.hp_steam_temperature.display()
        self.st.hp_steam_temperature.fix()
        self.hrsg.pump_hp.outlet.pressure[0].unfix()
        solver_obj.solve(self, tee=True)

        """
        print("set steam temperature 855")
        self.st.hp_steam_temperature.fix(855)
        solver_obj.solve(self, tee=True)
        self.gt.cmp1.efficiency_isentropic.fix(0.85)
        solver_obj.solve(self, tee=True)
        """
        """
        #power to 646
        self.net_power_mw.fix(646)
        self.gt.gt_power.unfix()
        solver_obj.solve(self, tee=True)

        #steam temp to 858
        self.st.hp_steam_temperature.fix(858)
        solver_obj.solve(self, tee=True)
        """
        if save_to is not None:
            iutil.to_json(self)
            if save_to is not None:
                iutil.to_json(self, fname=save_to)
                init_log.info_high(f"Initialization saved to {save_to}")
