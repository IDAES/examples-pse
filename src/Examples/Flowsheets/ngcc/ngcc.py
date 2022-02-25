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
        self.ng_preheat = um.HeatExchanger(
            default={
                "shell": {"property_package": self.hrsg.prop_water},
                "tube": {"property_package": self.gt.cmb_prop_params}
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
            source=self.hrsg.sh_hp4.side_1_outlet,
            destination=self.st.steam_turbine.inlet_split.inlet
        )
        self.fuel01dupe = Arc(
            source=self.ng_preheat.tube_outlet,
            destination=self.gt.inject1.gas
        )
        self.st01 = Arc(
            source=self.hrsg.splitter_ip1.toNGPH,
            destination=self.ng_preheat.shell_inlet
        )
        self.st02 = Arc(
            source=self.ng_preheat.shell_outlet,
            destination=self.hrsg.mixer1.Preheater
        )
        pyo.TransformationFactory("network.expand_arcs").apply_to(self)

    def _add_constraints(self):
        self.net_power_mw = pyo.Var(self.time, initialize=600, units=pyo.units.MW)

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


    def _scaling_guess(self):
        for (t, i) in self.fg_translate.mol_frac_eqn:
            iscale.constraint_scaling_transform(
                self.fg_translate.mol_frac_eqn[t, i], 1e-2)
        for t in self.fg_translate.temperature_eqn:
            iscale.constraint_scaling_transform(
                self.fg_translate.temperature_eqn[t], 1e-2)
        for t in self.fg_translate.pressure_eqn:
            iscale.constraint_scaling_transform(
                self.fg_translate.pressure_eqn[t], 1e-5)
        iscale.set_scaling_factor(self.ng_preheat.shell.heat, 1e-5)
        iscale.set_scaling_factor(self.ng_preheat.tube.heat, 1e-5)
        iscale.set_scaling_factor(
            self.ng_preheat.overall_heat_transfer_coefficient, 1e-2)
        iscale.set_scaling_factor(self.ng_preheat.area, 1e-3)

    def initialize(
        self,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
        load_from="ngcc_init.json.gz",
        save_to="ngcc_init.json.gz",
    ):

        solver_obj = iutil.get_solver(solver, optarg)
        self.t01a_expanded.deactivate()
        self.gt.initialize()

        propagate_state(self.g08a)
        self.fg_translate.initialize()
        propagate_state(self.g08b, overwrite_fixed=True)

        self.hrsg.initialize()
        solver_obj.solve(self.hrsg, tee=True)
        self.hrsg.sh_hp4.shell_inlet.unfix()
        propagate_state(self.t05a, overwrite_fixed=True)

        self.st.initialize()

        self.ng_preheat.area.fix(5000)
        self.ng_preheat.overall_heat_transfer_coefficient.fix(100)
        self.ng_preheat.tube_inlet.fix()
        propagate_state(self.st01)
        propagate_state(
            source=self.gt.feed_fuel1.outlet,
            destination=self.ng_preheat.tube_inlet,
            overwrite_fixed=True
        )
        self.ng_preheat.tube_inlet.temperature.fix(290)
        self.ng_preheat.initialize()
        #self.hrsg.mixer1.Preheater.unfix()
        self.st02_expanded.deactivate()
        self.ng_preheat.tube_inlet.flow_mol.unfix()
        self.gt.feed_fuel1.flow_mol.unfix()
        self.gt.feed_fuel1.pressure.unfix()
        self.gt.feed_fuel1.temperature.unfix()
        self.gt.feed_fuel1.mole_frac_comp.unfix()

        self.st.steam_turbine_lp_mix.hrsg.unfix()

        self.st.t02_dummy.deactivate()
        self.st.t03_dummy.deactivate()
        self.st.dummy_reheat.deactivate()
        self.hrsg.splitter_ip2.inlet.unfix()
        propagate_state(arc=self.t02a)
        propagate_state(arc=self.t03a)
        self.st.steam_turbine.throttle_valve[1].deltaP[0].fix(-5e5)
        for i, c in self.st.steam_turbine.throttle_valve.items():
            c.pressure_flow_equation.deactivate()
        @self.st.steam_turbine.Constraint([2,3,4], self.time)
        def same_throttle_pos_eqn(b, i, t):
            return b.throttle_valve[i].deltaP[t] == b.throttle_valve[1].deltaP[t]
        for i in self.st.steam_turbine.same_throttle_pos_eqn:
            iscale.constraint_scaling_transform(
                self.st.steam_turbine.same_throttle_pos_eqn[i], 1e-5)
        self.hrsg.econ_lp.tube_inlet.unfix()
        solver_obj.solve(self, tee=True)

        self.hrsg.mixer1.Preheater.unfix()
        self.st02_expanded.activate()
        solver_obj.solve(self, tee=True)


        self.t01a_expanded.activate()
        self.st.cond_pump.control_volume.properties_out[0].pressure.fix(
            pyo.value(self.hrsg.econ_lp.tube_inlet.pressure[0]))
        self.st.cond_pump.deltaP.unfix()
        self.st.hp_steam_temperature.unfix()
        self.st.hotwell.makeup.flow_mol[:].unfix()
        self.st.steam_turbine.inlet_split.inlet.unfix()
        self.hrsg.evap_lp.vapor_frac_control.fix(0.12)
        solver_obj.solve(self, tee=True)

        self.gt.gt_power[0].fix(-477e6)
        self.gt.exhaust_1.temperature.fix(898)
        self.gt.gts3.control_volume.properties_out[0].pressure.fix(103421)
        solver_obj.solve(self, tee=True)
