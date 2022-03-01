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
NGCC HRSG Subsystem for a 690MWe
"""
import os

import pandas as pd

from pyomo.common.fileutils import this_file_dir
import pyomo.environ as pyo
from pyomo.network import Arc

import idaes.logger as idaeslog
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state
from idaes.core import FlowsheetBlockData, declare_process_block_class
import idaes.core.util as iutil
import idaes.core.util.scaling as iscale
from idaes.generic_models.unit_models import Mixer, Separator as Splitter
from idaes.generic_models.properties import iapws95
from idaes.power_generation.properties import FlueGasParameterBlock
from idaes.power_generation.unit_models.helm import (
    HelmMixer,
    MomentumMixingType,
    HelmSplitter,
    HelmIsentropicCompressor as WaterPump,
)
from idaes.generic_models.unit_models.heat_exchanger import (
    HeatExchanger,
    HeatExchangerFlowPattern,
    delta_temperature_underwood_callback,
)
from idaes.power_generation.unit_models.boiler_heat_exchanger import (
    BoilerHeatExchanger,
    TubeArrangement,
)
from idaes.power_generation.unit_models.helm.phase_separator import HelmPhaseSeparator
import idaes.core.util.tables as ta
from idaes.core.util.tags import svg_tag

__author__ = "M. Zamarripa"


@declare_process_block_class(
    "HrsgFlowsheet",
    doc=(
        "The HRSG Flowsheet base on NETL report 'Cost and Performance Baseline "
        "for Fossil Energy Plants Volume 1: Bituminous Coal and Natural Gas to "
        "Electricity.' Sept 2019, Case B31B."
    ),
)
class HrsgFlowsheetData(FlowsheetBlockData):
    def build(self):
        super().build()
        self.prop_water = iapws95.Iapws95ParameterBlock()
        self.prop_gas = FlueGasParameterBlock(
            default={"components": ["N2", "O2", "CO2", "H2O"]}
        )
        self._add_unit_models()
        self._add_flowsheet_constraints()
        self._add_arcs()
        self._initial_design_variables()
        self._set_scaling_factors()
        self._stream_tags()

    def _add_unit_models(self):
        prop_water = self.prop_water
        prop_gas = self.prop_gas

        ######### LP Section ###########
        self.econ_lp = BoilerHeatExchanger(
            doc="Economizer",
            default={
                "delta_temperature_callback": delta_temperature_underwood_callback,
                "tube": {"property_package": prop_water},
                "shell": {"property_package": prop_gas},
                "has_pressure_change": False,
                "has_holdup": False,
                "flow_pattern": HeatExchangerFlowPattern.countercurrent,
                "tube_arrangement": TubeArrangement.inLine,
                "side_1_water_phase": "Liq",
                "has_radiation": False,
            },
        )
        self.mixer1 = HelmMixer(
            doc="Mixer for econ_lp outlet and preheater streams",
            default={
                "dynamic": False,
                "property_package": prop_water,
                "momentum_mixing_type": MomentumMixingType.none,
                "inlet_list": ["econ_lp", "Preheater"],
            },
        )
        self.drum_lp = HelmPhaseSeparator(
            doc="LP drum, flash seperator",
            default={"property_package": prop_water},
        )
        self.evap_lp = HeatExchanger(
            doc="LP evaporator",
            default={
                "shell": {"property_package": prop_gas},
                "tube": {"property_package": prop_water},
                "delta_temperature_callback": delta_temperature_underwood_callback,
                "flow_pattern": HeatExchangerFlowPattern.countercurrent,
            },
        )
        self.mixer_soec = HelmMixer(
            doc="Mixer for econ_lp outlet and preheater streams",
            default={
                "dynamic": False,
                "property_package": prop_water,
                "momentum_mixing_type": MomentumMixingType.none,
                "inlet_list": ["main", "soec_makeup"],
            },
        )
        self.split_fg_lp = Splitter(
            doc="LP section bypass flue gas splitter",
            default={
                "property_package": prop_gas,
                "ideal_separation": False,
                "outlet_list": ["toLP_SH", "toMixer"],
            },
        )
        self.mixer_lp2 = Mixer(
            doc="LP section flue gas bypass mixer",
            default={
                "property_package": prop_gas,
                "inlet_list": ["fromLP_SH", "bypass"],
            },
        )
        self.sh_lp = BoilerHeatExchanger(
            doc="LP superheater",
            default={
                "delta_temperature_callback": delta_temperature_underwood_callback,
                "tube": {"property_package": prop_water},
                "shell": {"property_package": prop_gas},
                "has_pressure_change": False,
                "has_holdup": False,
                "flow_pattern": HeatExchangerFlowPattern.countercurrent,
                "tube_arrangement": TubeArrangement.inLine,
                "side_1_water_phase": "Vap",
                "has_radiation": False,
            },
        )
        self.splitter1 = HelmSplitter(
            doc="LP liquid split to IP and HP pumps",
            default={"property_package": prop_water, "outlet_list": ["toIP", "toHP"]},
        )
        self.pump_ip = WaterPump(
            doc="Intermediate pressure pump",
            default={
                "property_package": prop_water,
            },
        )
        self.pump_hp = WaterPump(
            doc="High pressure pump",
            default={
                "property_package": prop_water,
            },
        )
        ######### IP Section ###########
        self.econ_ip1 = BoilerHeatExchanger(
            doc="IP ecomonmizer part 1",
            default={
                "delta_temperature_callback": delta_temperature_underwood_callback,
                "tube": {"property_package": prop_water},
                "shell": {"property_package": prop_gas},
                "has_pressure_change": True,
                "has_holdup": False,
                "flow_pattern": HeatExchangerFlowPattern.countercurrent,
                "tube_arrangement": TubeArrangement.inLine,
                "side_1_water_phase": "Liq",
                "has_radiation": False,
            },
        )
        self.splitter_ip1 = HelmSplitter(
            doc="IP economizer split for natrual gas preheater",
            default={
                "property_package": prop_water,
                "outlet_list": ["toIP_ECON2", "toNGPH"],
            },
        )
        self.econ_ip2 = BoilerHeatExchanger(
            doc="IP ecomonmizer part 2",
            default={
                "delta_temperature_callback": delta_temperature_underwood_callback,
                "tube": {"property_package": prop_water},
                "shell": {"property_package": prop_gas},
                "has_pressure_change": True,
                "has_holdup": False,
                "flow_pattern": HeatExchangerFlowPattern.countercurrent,
                "tube_arrangement": TubeArrangement.inLine,
                "side_1_water_phase": "Liq",
                "has_radiation": False,
            },
        )
        self.evap_ip = HeatExchanger(
            doc="IP evaporator",
            default={
                "shell": {"property_package": prop_gas},
                "tube": {"property_package": prop_water, "has_pressure_change": True},
                "delta_temperature_callback": delta_temperature_underwood_callback,
                "flow_pattern": HeatExchangerFlowPattern.countercurrent,
            },
        )
        self.sh_ip1 = BoilerHeatExchanger(
            doc="IP superheater 1",
            default={
                "delta_temperature_callback": delta_temperature_underwood_callback,
                "tube": {"property_package": prop_water},
                "shell": {"property_package": prop_gas},
                "has_pressure_change": True,
                "has_holdup": False,
                "flow_pattern": HeatExchangerFlowPattern.countercurrent,
                "tube_arrangement": TubeArrangement.inLine,
                "side_1_water_phase": "Vap",
                "has_radiation": False,
            },
        )
        self.mixer_ip1 = HelmMixer(
            doc="Mixer for econ_lp outlet and Preheater streams",
            default={
                "dynamic": False,
                "property_package": prop_water,
                "momentum_mixing_type": MomentumMixingType.none,
                "inlet_list": ["sh_ip1", "Cold_reheat"],
            },
        )
        self.splitter_ip2 = HelmSplitter(
            doc="IP Splitter 2, for ejector, reclaimer and dryer",
            default={
                "property_package": prop_water,
                "outlet_list": ["Cold_reheat", "toEjector", "toReclaimer", "toDryer"],
            },
        )
        self.sh_ip2 = BoilerHeatExchanger(
            doc="IP superheater 2",
            default={
                "delta_temperature_callback": delta_temperature_underwood_callback,
                "tube": {"property_package": prop_water},
                "shell": {"property_package": prop_gas},
                "has_pressure_change": True,
                "has_holdup": False,
                "flow_pattern": HeatExchangerFlowPattern.countercurrent,
                "tube_arrangement": TubeArrangement.inLine,
                "side_1_water_phase": "Vap",
                "has_radiation": False,
            },
        )
        self.sh_ip3 = BoilerHeatExchanger(
            doc="IP superheater 3",
            default={
                "delta_temperature_callback": delta_temperature_underwood_callback,
                "tube": {"property_package": prop_water},
                "shell": {"property_package": prop_gas},
                "has_pressure_change": True,
                "has_holdup": False,
                "flow_pattern": HeatExchangerFlowPattern.countercurrent,
                "tube_arrangement": TubeArrangement.inLine,
                "side_1_water_phase": "Vap",
                "has_radiation": False,
            },
        )
        self.econ_hp1 = BoilerHeatExchanger(
            doc="HP economizer 1",
            default={
                "delta_temperature_callback": delta_temperature_underwood_callback,
                "tube": {"property_package": prop_water},
                "shell": {"property_package": prop_gas},
                "has_pressure_change": True,
                "has_holdup": False,
                "flow_pattern": HeatExchangerFlowPattern.countercurrent,
                "tube_arrangement": TubeArrangement.inLine,
                "side_1_water_phase": "Liq",
                "has_radiation": False,
            },
        )
        self.econ_hp2 = BoilerHeatExchanger(
            doc="HP economizer 2",
            default={
                "delta_temperature_callback": delta_temperature_underwood_callback,
                "tube": {"property_package": prop_water},
                "shell": {"property_package": prop_gas},
                "has_pressure_change": True,
                "has_holdup": False,
                "flow_pattern": HeatExchangerFlowPattern.countercurrent,
                "tube_arrangement": TubeArrangement.inLine,
                "side_1_water_phase": "Liq",
                "has_radiation": False,
            },
        )
        self.econ_hp3 = BoilerHeatExchanger(
            doc="HP economizer 3",
            default={
                "delta_temperature_callback": delta_temperature_underwood_callback,
                "tube": {"property_package": prop_water},
                "shell": {"property_package": prop_gas},
                "has_pressure_change": True,
                "has_holdup": False,
                "flow_pattern": HeatExchangerFlowPattern.countercurrent,
                "tube_arrangement": TubeArrangement.inLine,
                "side_1_water_phase": "Liq",
                "has_radiation": False,
            },
        )
        self.econ_hp4 = BoilerHeatExchanger(
            doc="HP economizer 4",
            default={
                "delta_temperature_callback": delta_temperature_underwood_callback,
                "tube": {"property_package": prop_water},
                "shell": {"property_package": prop_gas},
                "has_pressure_change": True,
                "has_holdup": False,
                "flow_pattern": HeatExchangerFlowPattern.countercurrent,
                "tube_arrangement": TubeArrangement.inLine,
                "side_1_water_phase": "Liq",
                "has_radiation": False,
            },
        )
        self.econ_hp5 = BoilerHeatExchanger(
            doc="HP economizer 5",
            default={
                "delta_temperature_callback": delta_temperature_underwood_callback,
                "tube": {"property_package": prop_water},
                "shell": {"property_package": prop_gas},
                "has_pressure_change": True,
                "has_holdup": False,
                "flow_pattern": HeatExchangerFlowPattern.countercurrent,
                "tube_arrangement": TubeArrangement.inLine,
                "side_1_water_phase": "Liq",
                "has_radiation": False,
            },
        )
        self.evap_hp = HeatExchanger(
            doc="HP evaporator",
            default={
                "shell": {"property_package": prop_gas},
                "tube": {"property_package": prop_water},
                "delta_temperature_callback": delta_temperature_underwood_callback,
                "flow_pattern": HeatExchangerFlowPattern.countercurrent,
            },
        )
        self.sh_hp1 = BoilerHeatExchanger(
            doc="HP superheater 1",
            default={
                "delta_temperature_callback": delta_temperature_underwood_callback,
                "tube": {"property_package": prop_water},
                "shell": {"property_package": prop_gas},
                "has_pressure_change": True,
                "has_holdup": False,
                "flow_pattern": HeatExchangerFlowPattern.countercurrent,
                "tube_arrangement": TubeArrangement.inLine,
                "side_1_water_phase": "Vap",
                "has_radiation": False,
            },
        )
        self.sh_hp2 = BoilerHeatExchanger(
            doc="HP superheater 2",
            default={
                "delta_temperature_callback": delta_temperature_underwood_callback,
                "tube": {"property_package": prop_water},
                "shell": {"property_package": prop_gas},
                "has_pressure_change": True,
                "has_holdup": False,
                "flow_pattern": HeatExchangerFlowPattern.countercurrent,
                "tube_arrangement": TubeArrangement.inLine,
                "side_1_water_phase": "Vap",
                "has_radiation": False,
            },
        )
        self.sh_hp3 = BoilerHeatExchanger(
            doc="HP superheater 3",
            default={
                "delta_temperature_callback": delta_temperature_underwood_callback,
                "tube": {"property_package": prop_water},
                "shell": {"property_package": prop_gas},
                "has_pressure_change": True,
                "has_holdup": False,
                "flow_pattern": HeatExchangerFlowPattern.countercurrent,
                "tube_arrangement": TubeArrangement.inLine,
                "side_1_water_phase": "Vap",
                "has_radiation": False,
            },
        )
        self.sh_hp4 = BoilerHeatExchanger(
            doc="HP superheater 4",
            default={
                "delta_temperature_callback": delta_temperature_underwood_callback,
                "tube": {"property_package": prop_water},
                "shell": {"property_package": prop_gas},
                "has_pressure_change": True,
                "has_holdup": False,
                "flow_pattern": HeatExchangerFlowPattern.countercurrent,
                "tube_arrangement": TubeArrangement.inLine,
                "side_1_water_phase": "Vap",
                "has_radiation": False,
            },
        )

    def _add_flowsheet_constraints(self):
        """Add additional flowsheet constraints."""
        #
        # LP evaporator performance fixed to vaporize ~18 % of water inlet
        #
        self.evap_lp.vapor_frac_control = pyo.Var(
            initialize=0.12, doc="parameter to determine vapor flowrate"
        )
        self.evap_lp.vapor_frac_control.fix(0.12)

        @self.evap_lp.Constraint(self.time)
        def lp_vap_frac_eqn(b, t):
            return (
                b.tube.properties_out[0].vapor_frac == self.evap_lp.vapor_frac_control
            )

        # heat transfer determined by above
        self.evap_lp.heat_transfer_equation.deactivate()

        #
        # Everything evaporates in IP evaporator
        #
        @self.evap_ip.Constraint(self.config.time)
        def ip_sat_vap_eqn(b, t):
            return (
                b.tube.properties_out[t].enth_mol
                == b.tube.properties_out[t].enth_mol_sat_phase["Vap"] + 30
            )

        # heat transfer determined by above
        self.evap_ip.heat_transfer_equation.deactivate()

        #
        # HP Inlet must be vaporized
        #
        @self.evap_hp.Constraint(self.config.time)
        def hp_sat_vap_eqn(b, t):
            return (
                b.tube.properties_out[t].enth_mol
                == b.tube.properties_out[t].enth_mol_sat_phase["Vap"] + 30
            )

        # heat transfer determined by above
        self.evap_hp.heat_transfer_equation.deactivate()

        @self.mixer1.Constraint(self.config.time)
        def mixer1_pressure_eqn(b, t):
            return b.mixed_state[t].pressure == b.econ_lp_state[t].pressure

        @self.mixer_soec.Constraint(self.config.time)
        def mixer_soec_pressure_eqn(b, t):
            return b.mixed_state[t].pressure == b.main_state[t].pressure

        @self.mixer_ip1.Constraint(self.config.time)
        def mixer_ip1_pressure_eqn(b, t):
            return b.mixed_state[t].pressure == b.sh_ip1_state[t].pressure

    def _add_arcs(self):
        ######### LP Section ###########
        self.lp02 = Arc(
            doc="econ_lp to mixer1",
            source=self.econ_lp.side_1_outlet,
            destination=self.mixer1.econ_lp,
        )
        self.lp04 = Arc(
            doc="mixer1 to evap_lp",
            source=self.mixer1.outlet,
            destination=self.evap_lp.tube_inlet,
        )
        self.lp05 = Arc(
            doc="evap_lp to drum_lp",
            source=self.evap_lp.tube_outlet,
            destination=self.mixer_soec.main,
        )
        self.lp13 = Arc(
            doc="evap_lp to drum_lp",
            source=self.mixer_soec.outlet,
            destination=self.drum_lp.inlet,
        )
        self.lp10 = Arc(
            doc="drum_lp vapor to sh_lp",
            source=self.drum_lp.vap_outlet,
            destination=self.sh_lp.side_1_inlet,
        )
        self.lp06 = Arc(
            doc="drum_lp liquid to splitter1",
            source=self.drum_lp.liq_outlet,
            destination=self.splitter1.inlet,
        )
        self.lp08 = Arc(
            doc="splitter1 to pump_ip",
            source=self.splitter1.toIP,
            destination=self.pump_ip.inlet,
        )
        self.lp09 = Arc(
            doc="splitter1 to pump_hp",
            source=self.splitter1.toHP,
            destination=self.pump_hp.inlet,
        )
        ######### IP Section ###########
        self.ip01 = Arc(
            doc="pump_ip to econ_ip1",
            source=self.pump_ip.outlet,
            destination=self.econ_ip1.side_1_inlet,
        )
        self.ip02 = Arc(
            doc="econ_ip1 to splitter_ip1",
            source=self.econ_ip1.side_1_outlet,
            destination=self.splitter_ip1.inlet,
        )
        self.ip03 = Arc(
            doc="splitter_ip1 to econ_ip2",
            source=self.splitter_ip1.toIP_ECON2,
            destination=self.econ_ip2.side_1_inlet,
        )
        self.ip05 = Arc(
            doc="econ_ip2 to evap_ip",
            source=self.econ_ip2.side_1_outlet,
            destination=self.evap_ip.tube_inlet,
        )
        self.ip06 = Arc(
            doc="evap_ip to IPSH1",
            source=self.evap_ip.tube_outlet,
            destination=self.sh_ip1.side_1_inlet,
        )
        self.ip07 = Arc(
            doc="sh_ip1 to mixer_ip1",
            source=self.sh_ip1.side_1_outlet,
            destination=self.mixer_ip1.sh_ip1,
        )
        self.ip15 = Arc(
            doc="splitter_ip2 to Cold_reheat",
            source=self.splitter_ip2.Cold_reheat,
            destination=self.mixer_ip1.Cold_reheat,
        )
        self.ip08 = Arc(
            doc="mixer_ip1 to sh_ip2",
            source=self.mixer_ip1.outlet,
            destination=self.sh_ip2.side_1_inlet,
        )
        self.ip09 = Arc(
            doc="sh_ip2 to IPSH3",
            source=self.sh_ip2.side_1_outlet,
            destination=self.sh_ip3.side_1_inlet,
        )
        ######### HP Section ###########
        self.hp01 = Arc(
            doc="pump_hp to HP_ECON",
            source=self.pump_hp.outlet,
            destination=self.econ_hp1.side_1_inlet,
        )
        self.hp02 = Arc(
            doc="HP_ECON to econ_hp2",
            source=self.econ_hp1.side_1_outlet,
            destination=self.econ_hp2.side_1_inlet,
        )
        self.hp03 = Arc(
            doc="econ_hp3 to econ_hp4",
            source=self.econ_hp2.side_1_outlet,
            destination=self.econ_hp3.side_1_inlet,
        )
        self.hp04 = Arc(
            doc="econ_hp3 to econ_hp4",
            source=self.econ_hp3.side_1_outlet,
            destination=self.econ_hp4.side_1_inlet,
        )
        self.hp05 = Arc(
            doc="econ_hp4 to econ_hp5",
            source=self.econ_hp4.side_1_outlet,
            destination=self.econ_hp5.side_1_inlet,
        )
        self.hp06 = Arc(
            doc="econ_hp5 to evap_hp",
            source=self.econ_hp5.side_1_outlet,
            destination=self.evap_hp.tube_inlet,
        )
        self.hp07 = Arc(
            doc="evap_hp to sh_hp1",
            source=self.evap_hp.tube_outlet,
            destination=self.sh_hp1.side_1_inlet,
        )
        self.hp08 = Arc(
            doc="sh_hp1 to sh_hp2",
            source=self.sh_hp1.side_1_outlet,
            destination=self.sh_hp2.side_1_inlet,
        )
        self.hp09 = Arc(
            doc="sh_hp2 to sh_hp3",
            source=self.sh_hp2.side_1_outlet,
            destination=self.sh_hp3.side_1_inlet,
        )
        self.hp10 = Arc(
            doc="Flue gas from sh_hp3 to sh_hp4",
            source=self.sh_hp3.side_1_outlet,
            destination=self.sh_hp4.side_1_inlet,
        )
        self.g09 = Arc(
            doc="Flue gas from sh_hp4 to sh_ip3",
            source=self.sh_hp4.side_2_outlet,
            destination=self.sh_ip3.side_2_inlet,
        )
        self.g10 = Arc(
            doc="Flue gas from sh_ip3 to sh_hp3",
            source=self.sh_ip3.side_2_outlet,
            destination=self.sh_hp3.side_2_inlet,
        )
        self.g11 = Arc(
            doc="Flue gas from sh_hp3 to sh_hp2",
            source=self.sh_hp3.side_2_outlet,
            destination=self.sh_hp2.side_2_inlet,
        )
        self.g12 = Arc(
            doc="Flue gas from sh_hp2 to sh_ip2",
            source=self.sh_hp2.side_2_outlet,
            destination=self.sh_ip2.side_2_inlet,
        )
        self.g13 = Arc(
            doc="Flue gas from sh_ip2 to sh_hp1",
            source=self.sh_ip2.side_2_outlet,
            destination=self.sh_hp1.side_2_inlet,
        )
        self.g14 = Arc(
            doc="Flue gas from sh_hp1 to evap_hp",
            source=self.sh_hp1.side_2_outlet,
            destination=self.evap_hp.shell_inlet,
        )
        self.g15 = Arc(
            doc="Flue gas from evap_hp to econ_hp5",
            source=self.evap_hp.shell_outlet,
            destination=self.econ_hp5.side_2_inlet,
        )
        self.g16 = Arc(
            doc="Flue gas from econ_hp5 to sh_ip1",
            source=self.econ_hp5.side_2_outlet,
            destination=self.sh_ip1.side_2_inlet,
        )
        self.g17 = Arc(
            doc="Flue gas from sh_ip1 to econ_hp4",
            source=self.sh_ip1.side_2_outlet,
            destination=self.econ_hp4.side_2_inlet,
        )
        self.g18 = Arc(
            doc="Flue gas from econ_hp4 to econ_hp3",
            source=self.econ_hp4.side_2_outlet,
            destination=self.econ_hp3.side_2_inlet,
        )
        self.g19 = Arc(
            doc="Flue gas from econ_hp3 to split_fg_lp",
            source=self.econ_hp3.side_2_outlet,
            destination=self.split_fg_lp.inlet,
        )
        self.g20 = Arc(
            doc="Flue gas from split_fg_lp to sh_lp",
            source=self.split_fg_lp.toLP_SH,
            destination=self.sh_lp.side_2_inlet,
        )
        self.g21 = Arc(
            doc="Flue gas from sh_lp to mixer_lp2",
            source=self.sh_lp.side_2_outlet,
            destination=self.mixer_lp2.fromLP_SH,
        )
        self.g22 = Arc(
            doc="Flue gas from split_fg_lp to mixer_lp2",
            source=self.split_fg_lp.toMixer,
            destination=self.mixer_lp2.bypass,
        )
        self.g23 = Arc(
            doc="Flue gas from mixer_lp2 to evap_ip",
            source=self.mixer_lp2.outlet,
            destination=self.evap_ip.shell_inlet,
        )
        self.g24 = Arc(
            doc="Flue gas from evap_ip to econ_ip2",
            source=self.evap_ip.shell_outlet,
            destination=self.econ_ip2.side_2_inlet,
        )
        self.g25 = Arc(
            doc="Flue gas from econ_ip2 to econ_hp2",
            source=self.econ_ip2.side_2_outlet,
            destination=self.econ_hp2.side_2_inlet,
        )
        self.g26 = Arc(
            doc="Flue gas from econ_hp2 to econ_ip1",
            source=self.econ_hp2.side_2_outlet,
            destination=self.econ_ip1.side_2_inlet,
        )
        self.g27 = Arc(
            doc="Flue gas from econ_ip1 to econ_hp1",
            source=self.econ_ip1.side_2_outlet,
            destination=self.econ_hp1.side_2_inlet,
        )
        self.g28 = Arc(
            doc="Flue gas from econ_hp1 to evap_lp",
            source=self.econ_hp1.side_2_outlet,
            destination=self.evap_lp.shell_inlet,
        )
        self.g29 = Arc(
            doc="Flue gas from evap_lp to econ_lp",
            source=self.evap_lp.shell_outlet,
            destination=self.econ_lp.side_2_inlet,
        )
        pyo.TransformationFactory("network.expand_arcs").apply_to(self)

    def _initial_design_variables(self):
        """Set the basic inputs for the unit models including geometry, design
        and operating variables.
        """
        ######### LP Section ###########
        self.econ_lp.tube_di.fix(0.038)
        self.econ_lp.tube_thickness.fix(0.003)
        self.econ_lp.pitch_x.fix(0.09)
        self.econ_lp.pitch_y.fix(0.09)
        self.econ_lp.tube_length.fix(94.08 / 2)
        self.econ_lp.tube_nrow.fix(4)
        self.econ_lp.tube_ncol.fix(78)
        self.econ_lp.nrow_inlet.fix(2)
        self.econ_lp.delta_elevation.fix(1)
        self.econ_lp.tube_r_fouling = 0.000176
        self.econ_lp.shell_r_fouling = 0.00088
        if self.econ_lp.config.has_radiation:
            self.econ_lp.emissivity_wall.fix(0.7)
        self.econ_lp.fcorrection_htc.fix(1.5)
        self.econ_lp.fcorrection_dp_tube.fix(1.0)
        self.econ_lp.fcorrection_dp_shell.fix(1.0)
        self.splitter1.split_fraction[0, "toIP"].fix(0.20626)

        self.pump_ip.efficiency_isentropic.fix(0.80)
        self.pump_hp.efficiency_isentropic.fix(0.80)

        self.evap_lp.area.fix(42032.0 * 0.5)

        self.mixer_soec.soec_makeup.flow_mol.fix(0)
        self.mixer_soec.soec_makeup.pressure.fix(8*pyo.units.bar)
        self.mixer_soec.soec_makeup.enth_mol.fix(
            iapws95.htpx(P=8*pyo.units.bar, x=0.2))

        self.sh_lp.tube_di.fix(0.051)
        self.sh_lp.tube_thickness.fix(0.003)
        self.sh_lp.pitch_x.fix(0.11)
        self.sh_lp.pitch_y.fix(0.11)
        self.sh_lp.tube_length.fix(7 * 9)
        self.sh_lp.tube_nrow.fix(30)
        self.sh_lp.tube_ncol.fix(10 * 3)
        self.sh_lp.nrow_inlet.fix(10)
        self.sh_lp.delta_elevation.fix(1)
        self.sh_lp.tube_r_fouling = 0.000176
        self.sh_lp.shell_r_fouling = 0.00088
        if self.sh_lp.config.has_radiation:
            self.sh_lp.emissivity_wall.fix(0.7)
        self.sh_lp.fcorrection_htc.fix(1.1)
        self.sh_lp.fcorrection_dp_tube.fix(1.0)
        self.sh_lp.fcorrection_dp_shell.fix(1.0)

        ######### IP Section ###########
        self.econ_ip1.tube_di.fix(0.038)
        self.econ_ip1.tube_thickness.fix(0.003)
        self.econ_ip1.pitch_x.fix(0.09)
        self.econ_ip1.pitch_y.fix(0.09)
        self.econ_ip1.tube_length.fix(7 * 14)
        self.econ_ip1.tube_nrow.fix(6)
        self.econ_ip1.tube_ncol.fix(12)
        self.econ_ip1.nrow_inlet.fix(3)
        self.econ_ip1.delta_elevation.fix(1.0)
        self.econ_ip1.tube_r_fouling = 0.000176
        self.econ_ip1.shell_r_fouling = 0.00088
        if self.econ_ip1.config.has_radiation:
            self.econ_ip1.emissivity_wall.fix(0.7)
        self.econ_ip1.fcorrection_htc.fix(1.9)
        self.econ_ip1.fcorrection_dp_tube.fix(1.0)
        self.econ_ip1.fcorrection_dp_shell.fix(1.0)

        self.splitter_ip1.split_fraction[0, "toNGPH"].fix(0.4592)

        self.econ_ip2.tube_di.fix(0.051)
        self.econ_ip2.tube_thickness.fix(0.003)
        self.econ_ip2.pitch_x.fix(0.09)
        self.econ_ip2.pitch_y.fix(0.09)
        self.econ_ip2.tube_length.fix(7)
        self.econ_ip2.tube_nrow.fix(15)
        self.econ_ip2.tube_ncol.fix(71)
        self.econ_ip2.nrow_inlet.fix(2)
        self.econ_ip2.delta_elevation.fix(12)
        self.econ_ip2.tube_r_fouling = 0.000176
        self.econ_ip2.shell_r_fouling = 0.00088
        if self.econ_ip2.config.has_radiation:
            self.econ_ip2.emissivity_wall.fix(0.7)
        self.econ_ip2.fcorrection_htc.fix(1.0)
        self.econ_ip2.fcorrection_dp_tube.fix(1.0)
        self.econ_ip2.fcorrection_dp_shell.fix(1.0)

        self.evap_ip.area.fix(8368.6)
        self.evap_ip.overall_heat_transfer_coefficient.fix(150)

        self.sh_ip1.tube_di.fix(0.051)
        self.sh_ip1.tube_thickness.fix(0.003)
        self.sh_ip1.pitch_x.fix(0.11)
        self.sh_ip1.pitch_y.fix(0.11)
        self.sh_ip1.tube_length.fix(7 * 4)
        self.sh_ip1.tube_nrow.fix(10)
        self.sh_ip1.tube_ncol.fix(8)
        self.sh_ip1.nrow_inlet.fix(5)
        self.sh_ip1.delta_elevation.fix(1)
        self.sh_ip1.tube_r_fouling = 0.000176
        self.sh_ip1.shell_r_fouling = 0.00088
        if self.sh_ip1.config.has_radiation:
            self.sh_ip1.emissivity_wall.fix(0.7)
        self.sh_ip1.fcorrection_htc.fix(0.85*1.45)
        self.sh_ip1.fcorrection_dp_tube.fix(1.0)
        self.sh_ip1.fcorrection_dp_shell.fix(1.0)

        self.sh_ip2.tube_di.fix(0.0635)
        self.sh_ip2.tube_thickness.fix(0.003)
        self.sh_ip2.pitch_x.fix(0.11)
        self.sh_ip2.pitch_y.fix(0.11)
        self.sh_ip2.tube_length.fix(7 * 4)
        self.sh_ip2.tube_nrow.fix(16)
        self.sh_ip2.tube_ncol.fix(37)
        self.sh_ip2.nrow_inlet.fix(4)
        self.sh_ip2.delta_elevation.fix(1)
        self.sh_ip2.tube_r_fouling = 0.000176
        self.sh_ip2.shell_r_fouling = 0.00088
        if self.sh_ip2.config.has_radiation:
            self.sh_ip2.emissivity_wall.fix(0.7)
        self.sh_ip2.fcorrection_htc.fix(0.85*1.45)
        self.sh_ip2.fcorrection_dp_tube.fix(1.0)
        self.sh_ip2.fcorrection_dp_shell.fix(1.0)

        self.sh_ip3.tube_di.fix(0.0635)
        self.sh_ip3.tube_thickness.fix(0.003)
        self.sh_ip3.pitch_x.fix(0.11)
        self.sh_ip3.pitch_y.fix(0.11)
        self.sh_ip3.tube_length.fix(7 * 10)
        self.sh_ip3.tube_nrow.fix(16)
        self.sh_ip3.tube_ncol.fix(37)
        self.sh_ip3.nrow_inlet.fix(4)
        self.sh_ip3.delta_elevation.fix(1)
        self.sh_ip3.tube_r_fouling = 0.000176
        self.sh_ip3.shell_r_fouling = 0.00088
        if self.sh_ip3.config.has_radiation:
            self.sh_ip3.emissivity_wall.fix(0.7)
        self.sh_ip3.fcorrection_htc.fix(0.95*1.45)
        self.sh_ip3.fcorrection_dp_tube.fix(1.0)
        self.sh_ip3.fcorrection_dp_shell.fix(1.0)

        ######### HP Section ###########
        self.econ_hp1.tube_di.fix(0.051)
        self.econ_hp1.tube_thickness.fix(0.003)
        self.econ_hp1.pitch_x.fix(0.1)
        self.econ_hp1.pitch_y.fix(0.1)
        self.econ_hp1.tube_length.fix(7 * 10)
        self.econ_hp1.tube_nrow.fix(8)
        self.econ_hp1.tube_ncol.fix(61)
        self.econ_hp1.nrow_inlet.fix(4)
        self.econ_hp1.delta_elevation.fix(1.0)
        self.econ_hp1.tube_r_fouling = 0.000176
        self.econ_hp1.shell_r_fouling = 0.00088
        if self.econ_hp1.config.has_radiation:
            self.econ_hp1.emissivity_wall.fix(0.7)
        self.econ_hp1.fcorrection_htc.fix(1.0)
        self.econ_hp1.fcorrection_dp_tube.fix(0.05)
        self.econ_hp1.fcorrection_dp_shell.fix(0.05)

        self.econ_hp2.tube_di.fix(0.051)
        self.econ_hp2.tube_thickness.fix(0.003)
        self.econ_hp2.pitch_x.fix(0.1)
        self.econ_hp2.pitch_y.fix(0.1)
        self.econ_hp2.tube_length.fix(7 * 10)
        self.econ_hp2.tube_nrow.fix(8)
        self.econ_hp2.tube_ncol.fix(61)
        self.econ_hp2.nrow_inlet.fix(4)
        self.econ_hp2.delta_elevation.fix(1.0)
        self.econ_hp2.tube_r_fouling = 0.000176
        self.econ_hp2.shell_r_fouling = 0.00088
        if self.econ_hp2.config.has_radiation:
            self.econ_hp2.emissivity_wall.fix(0.7)
        self.econ_hp2.fcorrection_htc.fix(1.0)
        self.econ_hp2.fcorrection_dp_tube.fix(0.05)
        self.econ_hp2.fcorrection_dp_shell.fix(0.05)

        self.econ_hp3.tube_di.fix(0.051)
        self.econ_hp3.tube_thickness.fix(0.003)
        self.econ_hp3.pitch_x.fix(0.1)
        self.econ_hp3.pitch_y.fix(0.1)
        self.econ_hp3.tube_length.fix(7 * 12)
        self.econ_hp3.tube_nrow.fix(8)
        self.econ_hp3.tube_ncol.fix(61)
        self.econ_hp3.nrow_inlet.fix(4)
        self.econ_hp3.delta_elevation.fix(1.0)
        self.econ_hp3.tube_r_fouling = 0.000176
        self.econ_hp3.shell_r_fouling = 0.00088
        if self.econ_hp3.config.has_radiation:
            self.econ_hp3.emissivity_wall.fix(0.7)
        self.econ_hp3.fcorrection_htc.fix(1.2)
        self.econ_hp3.fcorrection_dp_tube.fix(0.05)
        self.econ_hp3.fcorrection_dp_shell.fix(0.05)

        self.econ_hp4.tube_di.fix(0.051)
        self.econ_hp4.tube_thickness.fix(0.003)
        self.econ_hp4.pitch_x.fix(0.1)
        self.econ_hp4.pitch_y.fix(0.1)
        self.econ_hp4.tube_length.fix(7 * 10)
        self.econ_hp4.tube_nrow.fix(8)
        self.econ_hp4.tube_ncol.fix(61)
        self.econ_hp4.nrow_inlet.fix(4)
        self.econ_hp4.delta_elevation.fix(1.0)
        self.econ_hp4.tube_r_fouling = 0.000176
        self.econ_hp4.shell_r_fouling = 0.00088
        if self.econ_hp4.config.has_radiation:
            self.econ_hp4.emissivity_wall.fix(0.7)
        self.econ_hp4.fcorrection_htc.fix(1.0)
        self.econ_hp4.fcorrection_dp_tube.fix(0.05)
        self.econ_hp4.fcorrection_dp_shell.fix(0.05)

        self.econ_hp5.tube_di.fix(0.051)
        self.econ_hp5.tube_thickness.fix(0.003)
        self.econ_hp5.pitch_x.fix(0.1)
        self.econ_hp5.pitch_y.fix(0.1)
        self.econ_hp5.tube_length.fix(7 * 10)
        self.econ_hp5.tube_nrow.fix(8)
        self.econ_hp5.tube_ncol.fix(61)
        self.econ_hp5.nrow_inlet.fix(4)
        self.econ_hp5.delta_elevation.fix(1.0)
        self.econ_hp5.tube_r_fouling = 0.000176
        self.econ_hp5.shell_r_fouling = 0.00088
        if self.econ_hp5.config.has_radiation:
            self.econ_hp5.emissivity_wall.fix(0.7)
        self.econ_hp5.fcorrection_htc.fix(1.0)
        self.econ_hp5.fcorrection_dp_tube.fix(0.05)
        self.econ_hp5.fcorrection_dp_shell.fix(0.05)

        self.sh_hp1.tube_di.fix(0.0635)
        self.sh_hp1.tube_thickness.fix(0.003)
        self.sh_hp1.pitch_x.fix(0.11)
        self.sh_hp1.pitch_y.fix(0.11)
        self.sh_hp1.tube_length.fix(7 * 10)
        self.sh_hp1.tube_nrow.fix(16)
        self.sh_hp1.tube_ncol.fix(37)
        self.sh_hp1.nrow_inlet.fix(4)
        self.sh_hp1.delta_elevation.fix(1)
        self.sh_hp1.tube_r_fouling = 0.000176
        self.sh_hp1.shell_r_fouling = 0.00088
        if self.sh_hp1.config.has_radiation:
            self.sh_hp1.emissivity_wall.fix(0.7)
        self.sh_hp1.fcorrection_htc.fix(1.2*0.75)
        self.sh_hp1.fcorrection_dp_tube.fix(1.0)
        self.sh_hp1.fcorrection_dp_shell.fix(1.0)

        self.sh_hp2.tube_di.fix(0.0635)
        self.sh_hp2.tube_thickness.fix(0.003)
        self.sh_hp2.pitch_x.fix(0.11)
        self.sh_hp2.pitch_y.fix(0.11)
        self.sh_hp2.tube_length.fix(7 * 10)
        self.sh_hp2.tube_nrow.fix(16)
        self.sh_hp2.tube_ncol.fix(37)
        self.sh_hp2.nrow_inlet.fix(4)
        self.sh_hp2.delta_elevation.fix(1)
        self.sh_hp2.tube_r_fouling = 0.000176
        self.sh_hp2.shell_r_fouling = 0.00088
        if self.sh_hp2.config.has_radiation:
            self.sh_hp2.emissivity_wall.fix(0.7)
        self.sh_hp2.fcorrection_htc.fix(0.95*0.75)
        self.sh_hp2.fcorrection_dp_tube.fix(1.0)
        self.sh_hp2.fcorrection_dp_shell.fix(1.0)

        self.sh_hp3.tube_di.fix(0.0635)
        self.sh_hp3.tube_thickness.fix(0.003)
        self.sh_hp3.pitch_x.fix(0.11)
        self.sh_hp3.pitch_y.fix(0.11)
        self.sh_hp3.tube_length.fix(7 * 10)
        self.sh_hp3.tube_nrow.fix(16)
        self.sh_hp3.tube_ncol.fix(37)
        self.sh_hp3.nrow_inlet.fix(4)
        self.sh_hp3.delta_elevation.fix(1)
        self.sh_hp3.tube_r_fouling = 0.000176
        self.sh_hp3.shell_r_fouling = 0.00088
        if self.sh_hp3.config.has_radiation:
            self.sh_hp3.emissivity_wall.fix(0.7)
        self.sh_hp3.fcorrection_htc.fix(0.95*0.75)
        self.sh_hp3.fcorrection_dp_tube.fix(1.0)
        self.sh_hp3.fcorrection_dp_shell.fix(1.0)

        self.sh_hp4.tube_di.fix(0.0635)
        self.sh_hp4.tube_thickness.fix(0.003)
        self.sh_hp4.pitch_x.fix(0.11)
        self.sh_hp4.pitch_y.fix(0.11)
        self.sh_hp4.tube_length.fix(7 * 10)
        self.sh_hp4.tube_nrow.fix(16)
        self.sh_hp4.tube_ncol.fix(37)
        self.sh_hp4.nrow_inlet.fix(4)
        self.sh_hp4.delta_elevation.fix(1)
        self.sh_hp4.tube_r_fouling = 0.000176
        self.sh_hp4.shell_r_fouling = 0.00088
        if self.sh_hp4.config.has_radiation:
            self.sh_hp4.emissivity_wall.fix(0.7)
        self.sh_hp4.fcorrection_htc.fix(0.95*0.75)
        self.sh_hp4.fcorrection_dp_tube.fix(1.0)
        self.sh_hp4.fcorrection_dp_shell.fix(1.0)

    def _set_scaling_factors(self):
        hx_list = [
            self.econ_lp,
            self.sh_lp,
            self.econ_ip1,
            self.econ_ip2,
            self.sh_ip1,
            self.sh_ip2,
            self.sh_ip3,
            self.econ_hp1,
            self.econ_hp2,
            self.econ_hp3,
            self.econ_hp4,
            self.econ_hp5,
            self.sh_hp1,
            self.sh_hp2,
            self.sh_hp3,
            self.sh_hp4,
        ]
        for unit in hx_list:
            iscale.set_scaling_factor(unit.side_1.heat, 1e-6)
            iscale.set_scaling_factor(unit.side_2.heat, 1e-6)
            iscale.set_scaling_factor(unit.overall_heat_transfer_coefficient, 1e-2)
            iscale.set_scaling_factor(unit.area, 1e-5)
            iscale.set_scaling_factor(unit.area, 1e-5)
            iscale.set_scaling_factor(unit.N_Re_tube, 1e-5)
            iscale.set_scaling_factor(unit.N_Nu_tube, 1e-2)
            iscale.set_scaling_factor(unit.N_Re_shell, 1e-4)
            iscale.set_scaling_factor(unit.N_Nu_shell, 1e-2)
            iscale.set_scaling_factor(unit.hconv_tube, 1e-3)
            iscale.set_scaling_factor(unit.hconv_shell_conv, 1e-2)
            iscale.set_scaling_factor(unit.hconv_shell_total, 1e-2)
            iscale.set_scaling_factor(unit.rcond_wall, 1e5)
            if hasattr(unit, "deltaP_tube_friction"):
                iscale.set_scaling_factor(unit.deltaP_tube_friction, 1e-4)
                iscale.set_scaling_factor(unit.friction_factor_tube, 1e3)
                iscale.set_scaling_factor(unit.friction_factor_shell, 1e3)
            if hasattr(unit, "deltaP_tube_uturn"):
                iscale.set_scaling_factor(unit.deltaP_tube_uturn, 1e-2)

        iscale.set_scaling_factor(self.evap_lp.shell.heat, 1e-8)
        iscale.set_scaling_factor(self.evap_lp.tube.heat, 1e-8)
        iscale.set_scaling_factor(self.evap_lp.tube.heat, 1e-8)
        iscale.set_scaling_factor(self.evap_lp.area, 1e-4)
        iscale.set_scaling_factor(self.evap_lp.overall_heat_transfer_coefficient, 1e-2)

        iscale.set_scaling_factor(self.evap_ip.shell.heat, 1e-8)
        iscale.set_scaling_factor(self.evap_ip.tube.heat, 1e-8)
        iscale.set_scaling_factor(self.evap_ip.tube.heat, 1e-8)
        iscale.set_scaling_factor(self.evap_ip.area, 1e-4)
        iscale.set_scaling_factor(self.evap_ip.overall_heat_transfer_coefficient, 1e-2)

        iscale.set_scaling_factor(self.pump_ip.control_volume.work, 1e-6)
        iscale.set_scaling_factor(self.pump_hp.control_volume.work, 1e-7)

        iscale.set_scaling_factor(self.evap_hp.shell.heat, 1e-8)
        iscale.set_scaling_factor(self.evap_hp.tube.heat, 1e-8)
        iscale.set_scaling_factor(self.evap_hp.tube.heat, 1e-8)
        iscale.set_scaling_factor(self.evap_hp.area, 1e-4)
        iscale.set_scaling_factor(self.evap_hp.overall_heat_transfer_coefficient, 1e-2)

        for t in self.config.time:
            iscale.constraint_scaling_transform(self.mixer1.mixer1_pressure_eqn[t], 1e-6)
            iscale.constraint_scaling_transform(self.mixer_soec.mixer_soec_pressure_eqn[t], 1e-6)
            iscale.constraint_scaling_transform(self.mixer_ip1.mixer_ip1_pressure_eqn[t], 1e-6)

    def initialize(
        self,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
        load_from="hrsg_init.json.gz",
        save_to="hrsg_init.json.gz",
    ):
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="flowsheet")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="flowsheet")

        solver_obj = iutil.get_solver(solver, optarg)

        init_log.info_high("HRSG Initialization Starting")

        if load_from is not None:
            if os.path.exists(load_from):
                init_log.info_high(f"HRSG load initial from {load_from}")
                # here suffix=False avoids loading scaling factors
                iutil.from_json(
                    self, fname=load_from, wts=iutil.StoreSpec(suffix=False)
                )
                return

        ######### LP Section ###########
        # FLUE GAS Inlet to econ_lp (F = 138406 kgmol/hr, T = 15 C, P=0.1 MPa abs)
        fg_rate = 38446.11  # mol/s  from Baseline report table 5-22

        self.econ_lp.side_1_inlet.flow_mol[0].fix(9790.55)
        self.econ_lp.side_1_inlet.enth_mol[0].fix(
            iapws95.htpx(T=357.03 * pyo.units.K, P=599844 * pyo.units.Pa)
        )
        self.econ_lp.side_1_inlet.pressure[0].fix(655000)
        self.econ_lp.side_2_inlet.flow_mol_comp[0, "H2O"].fix(fg_rate * 0.0875)
        self.econ_lp.side_2_inlet.flow_mol_comp[0, "CO2"].fix(fg_rate * 0.0408)
        self.econ_lp.side_2_inlet.flow_mol_comp[0, "N2"].fix(fg_rate * 0.75)
        self.econ_lp.side_2_inlet.flow_mol_comp[0, "O2"].fix(fg_rate * 0.1217)
        self.econ_lp.side_2_inlet.temperature[0].fix(433)
        self.econ_lp.side_2_inlet.pressure[0].fix(103421)
        self.econ_lp.initialize(solver=solver, outlvl=outlvl, optarg=optarg)

        propagate_state(self.mixer1.econ_lp, self.econ_lp.side_1_outlet)
        self.mixer1.Preheater.flow_mol[0].fix(833)
        self.mixer1.Preheater.enth_mol[0].fix(
            iapws95.htpx(T=333.15 * pyo.units.K, P=42 * pyo.units.bar)
        )
        self.mixer1.Preheater.pressure[0].fix(3.509e6)
        self.mixer1.initialize(solver=solver, outlvl=outlvl, optarg=optarg)

        propagate_state(self.evap_lp.tube_inlet, self.mixer1.outlet)
        self.evap_lp.shell_inlet.flow_mol_comp[0, "H2O"].fix(fg_rate * 0.0875)
        self.evap_lp.shell_inlet.flow_mol_comp[0, "CO2"].fix(fg_rate * 0.0408)
        self.evap_lp.shell_inlet.flow_mol_comp[0, "N2"].fix(fg_rate * 0.75)
        self.evap_lp.shell_inlet.flow_mol_comp[0, "O2"].fix(fg_rate * 0.1217)
        self.evap_lp.shell_inlet.temperature[0].fix(543.31)
        self.evap_lp.shell_inlet.pressure[0].fix(103421)
        self.evap_lp.overall_heat_transfer_coefficient.fix(212)
        self.evap_lp.lp_vap_frac_eqn.deactivate()
        self.evap_lp.tube_inlet.flow_mol[0].fix()
        self.evap_lp.tube_inlet.enth_mol[0].fix()
        self.evap_lp.tube_inlet.pressure[0].fix()
        self.evap_lp.initialize(solver=solver, outlvl=outlvl, optarg=optarg)
        self.evap_lp.lp_vap_frac_eqn.activate()
        self.evap_lp.heat_transfer_equation.deactivate()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(self.evap_lp, tee=slc.tee)
        if not pyo.check_optimal_termination(res):
            raise InitializationError(f"evap_lp failed to initialize successfully.")

        propagate_state(self.lp13)
        self.mixer_soec.initialize

        propagate_state(self.drum_lp.inlet, self.evap_lp.tube_outlet)
        self.drum_lp.initialize()

        self.drum_lp.inlet.flow_mol[0].fix()
        self.drum_lp.inlet.enth_mol[0].fix()
        self.drum_lp.inlet.pressure[0].fix()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(self.drum_lp, tee=slc.tee)
        if not pyo.check_optimal_termination(res):
            raise InitializationError(f"drum_lp failed to initialize successfully.")

        propagate_state(self.splitter1.inlet, self.drum_lp.liq_outlet)
        self.splitter1.initialize(solver=solver, outlvl=outlvl, optarg=optarg)

        propagate_state(self.pump_ip.inlet, self.splitter1.toIP)
        self.pump_ip.outlet.pressure[0].fix(4.385e6)
        self.pump_ip.initialize(solver=solver, outlvl=outlvl, optarg=optarg)

        propagate_state(self.pump_hp.inlet, self.splitter1.toHP)
        self.pump_hp.outlet.pressure[0].fix(2.25e7)
        self.pump_hp.initialize(solver=solver, outlvl=outlvl, optarg=optarg)

        self.split_fg_lp.inlet.flow_mol_comp[0, "H2O"].fix(fg_rate * 0.0875)
        self.split_fg_lp.inlet.flow_mol_comp[0, "CO2"].fix(fg_rate * 0.0408)
        self.split_fg_lp.inlet.flow_mol_comp[0, "N2"].fix(fg_rate * 0.75)
        self.split_fg_lp.inlet.flow_mol_comp[0, "O2"].fix(fg_rate * 0.1217)
        self.split_fg_lp.inlet.temperature[0].fix(640.15)
        self.split_fg_lp.inlet.pressure[0].fix(103421)
        self.split_fg_lp.split_fraction[0, "toLP_SH"].fix(0.55)
        self.split_fg_lp.initialize(solver=solver, outlvl=outlvl, optarg=optarg)

        propagate_state(self.sh_lp.side_1_inlet, self.drum_lp.vap_outlet)
        propagate_state(self.sh_lp.side_2_inlet, self.split_fg_lp.toLP_SH)

        self.sh_lp.initialize(solver=solver, outlvl=outlvl, optarg=optarg)

        self.sh_lp.side_1_inlet.flow_mol[0].fix()
        self.sh_lp.side_1_inlet.enth_mol[0].fix()
        self.sh_lp.side_1_inlet.pressure[0].fix()

        init_log.info("Low pressure system initialization - Completed")

        propagate_state(self.mixer_lp2.fromLP_SH, self.sh_lp.side_2_inlet)
        propagate_state(self.mixer_lp2.bypass, self.split_fg_lp.toMixer)
        self.mixer_lp2.initialize(solver=solver, outlvl=outlvl, optarg=optarg)

        propagate_state(self.econ_ip1.side_1_inlet, self.pump_ip.outlet)
        self.econ_ip1.side_2_inlet.flow_mol_comp[0, "H2O"].fix(fg_rate * 0.0875)
        self.econ_ip1.side_2_inlet.flow_mol_comp[0, "CO2"].fix(fg_rate * 0.0408)
        self.econ_ip1.side_2_inlet.flow_mol_comp[0, "N2"].fix(fg_rate * 0.75)
        self.econ_ip1.side_2_inlet.flow_mol_comp[0, "O2"].fix(fg_rate * 0.1217)
        self.econ_ip1.side_2_inlet.temperature[0].fix(579.46)
        self.econ_ip1.side_2_inlet.pressure[0].fix(103421)
        self.econ_ip1.initialize(solver=solver, outlvl=outlvl, optarg=optarg)
        self.econ_ip1.deltaP_tube_eqn.deactivate()
        self.econ_ip1.deltaP_tube.fix(-1.79e05)
        self.econ_ip1.side_1_inlet.flow_mol[0].fix()
        self.econ_ip1.side_1_inlet.enth_mol[0].fix()
        self.econ_ip1.side_1_inlet.pressure[0].fix()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(self.econ_ip1, tee=slc.tee)
        if not pyo.check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} econ_ip1 failed to initialize successfully."
            )

        propagate_state(self.splitter_ip1.inlet, self.econ_ip1.side_1_outlet)
        self.splitter_ip1.initialize(solver=solver, outlvl=outlvl, optarg=optarg)

        propagate_state(self.econ_ip2.side_1_inlet, self.splitter_ip1.toIP_ECON2)
        self.econ_ip2.side_2_inlet.flow_mol_comp[0, "H2O"].fix(fg_rate * 0.0875)
        self.econ_ip2.side_2_inlet.flow_mol_comp[0, "CO2"].fix(fg_rate * 0.0408)
        self.econ_ip2.side_2_inlet.flow_mol_comp[0, "N2"].fix(fg_rate * 0.75)
        self.econ_ip2.side_2_inlet.flow_mol_comp[0, "O2"].fix(fg_rate * 0.1217)
        self.econ_ip2.side_2_inlet.temperature[0].fix(607)
        self.econ_ip2.side_2_inlet.pressure[0].fix(103421)
        self.econ_ip2.initialize(solver=solver, outlvl=outlvl, optarg=optarg)
        self.econ_ip2.deltaP_tube_eqn.deactivate()
        self.econ_ip2.deltaP_tube.fix(-1.63e05)
        self.econ_ip2.side_1_inlet.flow_mol[0].fix()
        self.econ_ip2.side_1_inlet.enth_mol[0].fix()
        self.econ_ip2.side_1_inlet.pressure[0].fix()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(self.econ_ip2, tee=False)
        if not pyo.check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} econ_ip2 failed to initialize successfully."
            )

        propagate_state(self.evap_ip.tube_inlet, self.econ_ip2.side_1_outlet)
        self.evap_ip.shell_inlet.flow_mol_comp[0, "H2O"].fix(fg_rate * 0.0875)
        self.evap_ip.shell_inlet.flow_mol_comp[0, "CO2"].fix(fg_rate * 0.0408)
        self.evap_ip.shell_inlet.flow_mol_comp[0, "N2"].fix(fg_rate * 0.75)
        self.evap_ip.shell_inlet.flow_mol_comp[0, "O2"].fix(fg_rate * 0.1217)
        self.evap_ip.shell_inlet.temperature[0].fix(618.0)
        self.evap_ip.shell_inlet.pressure[0].fix(103421)
        self.evap_ip.ip_sat_vap_eqn.deactivate()
        self.evap_ip.heat_transfer_equation.activate()
        self.evap_ip.tube.deltaP.fix(-1.62e05)
        self.evap_ip.initialize(solver=solver, outlvl=outlvl, optarg=optarg)
        self.evap_ip.tube_inlet.flow_mol[0].fix()
        self.evap_ip.tube_inlet.pressure[0].fix()
        self.evap_ip.tube_inlet.enth_mol[0].fix()
        self.evap_ip.overall_heat_transfer_coefficient.fix(150)
        self.evap_ip.ip_sat_vap_eqn.activate()
        self.evap_ip.heat_transfer_equation.deactivate()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(self.evap_ip, tee=False)
        if not pyo.check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} evap_ip failed to initialize successfully."
            )

        propagate_state(self.sh_ip1.side_1_inlet, self.evap_ip.tube_outlet)
        self.sh_ip1.side_2_inlet.flow_mol_comp[0, "H2O"].fix(fg_rate * 0.0875)
        self.sh_ip1.side_2_inlet.flow_mol_comp[0, "CO2"].fix(fg_rate * 0.0408)
        self.sh_ip1.side_2_inlet.flow_mol_comp[0, "N2"].fix(fg_rate * 0.75)
        self.sh_ip1.side_2_inlet.flow_mol_comp[0, "O2"].fix(fg_rate * 0.1217)
        self.sh_ip1.side_2_inlet.temperature[0].fix(668.15)
        self.sh_ip1.side_2_inlet.pressure[0].fix(103421)
        self.sh_ip1.initialize(solver=solver, outlvl=outlvl, optarg=optarg)
        self.sh_ip1.deltaP_tube_eqn.deactivate()
        self.sh_ip1.deltaP_tube.fix(-1.55e05)
        self.sh_ip1.side_1_inlet.flow_mol[0].fix()
        self.sh_ip1.side_1_inlet.enth_mol[0].fix()
        self.sh_ip1.side_1_inlet.pressure[0].fix()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(self.sh_ip1, tee=False)
        if not pyo.check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} sh_ip1 failed to initialize successfully."
            )

        self.splitter_ip2.inlet.flow_mol[0].fix(7503.7337)
        self.splitter_ip2.inlet.enth_mol[0].fix(
            iapws95.htpx(T=628.03 * pyo.units.K, P=3.737e6 * pyo.units.Pa)
        )
        self.splitter_ip2.inlet.pressure[0].fix(3.737e6)
        self.splitter_ip2.split_fraction[0, "Cold_reheat"].fix(0.9941)
        self.splitter_ip2.split_fraction[0, "toEjector"].fix(0.00074)
        self.splitter_ip2.split_fraction[0, "toDryer"].fix(0.000274)
        self.splitter_ip2.initialize(solver=solver, outlvl=outlvl, optarg=optarg)
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(self.splitter_ip2, tee=False)
        if not pyo.check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} splitter_ip2 failed to initialize successfully."
            )

        propagate_state(self.mixer_ip1.sh_ip1, self.sh_ip1.side_1_outlet)
        propagate_state(self.mixer_ip1.Cold_reheat, self.splitter_ip2.Cold_reheat)
        self.mixer_ip1.sh_ip1.flow_mol[0].fix()
        self.mixer_ip1.sh_ip1.enth_mol[0].fix()
        self.mixer_ip1.sh_ip1.pressure[0].fix()
        self.mixer_ip1.Cold_reheat.flow_mol[0].fix()
        self.mixer_ip1.Cold_reheat.enth_mol[0].fix()
        self.mixer_ip1.Cold_reheat.pressure[0].fix()
        self.mixer_ip1.initialize(solver=solver, outlvl=outlvl, optarg=optarg)

        propagate_state(self.sh_ip2.side_1_inlet, self.mixer_ip1.outlet)
        self.sh_ip2.side_2_inlet.flow_mol_comp[0, "H2O"].fix(fg_rate * 0.0875)
        self.sh_ip2.side_2_inlet.flow_mol_comp[0, "CO2"].fix(fg_rate * 0.0408)
        self.sh_ip2.side_2_inlet.flow_mol_comp[0, "N2"].fix(fg_rate * 0.75)
        self.sh_ip2.side_2_inlet.flow_mol_comp[0, "O2"].fix(fg_rate * 0.1217)
        self.sh_ip2.side_2_inlet.temperature[0].fix(802.35)  # 878.15)  # K
        self.sh_ip2.side_2_inlet.pressure[0].fix(103421)  # Pa
        self.sh_ip2.initialize(solver=solver, outlvl=outlvl, optarg=optarg)
        self.sh_ip2.deltaP_tube_eqn.deactivate()
        self.sh_ip2.deltaP_tube.fix(-7.45e04)
        self.sh_ip2.side_1_inlet.flow_mol[0].fix()
        self.sh_ip2.side_1_inlet.enth_mol[0].fix()
        self.sh_ip2.side_1_inlet.pressure[0].fix()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(self.sh_ip2, tee=False)
        if not pyo.check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} sh_ip2 failed to initialize successfully."
            )

        propagate_state(self.sh_ip3.side_1_inlet, self.sh_ip2.side_1_outlet)
        self.sh_ip3.side_2_inlet.flow_mol_comp[0, "H2O"].fix(fg_rate * 0.0875)
        self.sh_ip3.side_2_inlet.flow_mol_comp[0, "CO2"].fix(fg_rate * 0.0408)
        self.sh_ip3.side_2_inlet.flow_mol_comp[0, "N2"].fix(fg_rate * 0.75)
        self.sh_ip3.side_2_inlet.flow_mol_comp[0, "O2"].fix(fg_rate * 0.1217)
        self.sh_ip3.side_2_inlet.temperature[0].fix(868.35)
        self.sh_ip3.side_2_inlet.pressure[0].fix(103421)
        self.sh_ip3.initialize(solver=solver, outlvl=outlvl, optarg=optarg)
        self.sh_ip3.deltaP_tube_eqn.deactivate()
        self.sh_ip3.deltaP_tube.fix(-7.31e04)
        self.sh_ip3.side_1_inlet.flow_mol[0].fix()
        self.sh_ip3.side_1_inlet.enth_mol[0].fix()
        self.sh_ip3.side_1_inlet.pressure[0].fix()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(self.sh_ip3, tee=False)
        if not pyo.check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} sh_ip3 failed to initialize successfully."
            )
        init_log.info("Intermediate pressure system initialization - Completed")

        propagate_state(self.econ_hp1.side_1_inlet, self.pump_hp.outlet)
        self.econ_hp1.side_2_inlet.flow_mol_comp[0, "H2O"].fix(fg_rate * 0.0875)
        self.econ_hp1.side_2_inlet.flow_mol_comp[0, "CO2"].fix(fg_rate * 0.0408)
        self.econ_hp1.side_2_inlet.flow_mol_comp[0, "N2"].fix(fg_rate * 0.75)
        self.econ_hp1.side_2_inlet.flow_mol_comp[0, "O2"].fix(fg_rate * 0.1217)
        self.econ_hp1.side_2_inlet.temperature[0].fix(566.35)
        self.econ_hp1.side_2_inlet.pressure[0].fix(103421)
        self.econ_hp1.initialize(solver=solver, outlvl=outlvl, optarg=optarg)
        self.econ_hp1.deltaP_tube_eqn.deactivate()
        self.econ_hp1.deltaP_tube.fix(-9.79e05)
        self.econ_hp1.side_1_inlet.flow_mol[0].fix()
        self.econ_hp1.side_1_inlet.enth_mol[0].fix()
        self.econ_hp1.side_1_inlet.pressure[0].fix()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(self.econ_hp1, tee=False)
        if not pyo.check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} econ_hp1 failed to initialize successfully."
            )

        propagate_state(self.econ_hp2.side_1_inlet, self.econ_hp1.side_1_outlet)
        self.econ_hp2.side_2_inlet.flow_mol_comp[0, "H2O"].fix(fg_rate * 0.0875)
        self.econ_hp2.side_2_inlet.flow_mol_comp[0, "CO2"].fix(fg_rate * 0.0408)
        self.econ_hp2.side_2_inlet.flow_mol_comp[0, "N2"].fix(fg_rate * 0.75)
        self.econ_hp2.side_2_inlet.flow_mol_comp[0, "O2"].fix(fg_rate * 0.1217)
        self.econ_hp2.side_2_inlet.temperature[0].fix(596.35)
        self.econ_hp2.side_2_inlet.pressure[0].fix(103421)
        self.econ_hp2.initialize(solver=solver, outlvl=outlvl, optarg=optarg)
        self.econ_hp2.deltaP_tube_eqn.deactivate()
        self.econ_hp2.deltaP_tube.fix(-9.38e05)
        self.econ_hp2.side_1_inlet.flow_mol[0].fix()
        self.econ_hp2.side_1_inlet.enth_mol[0].fix()
        self.econ_hp2.side_1_inlet.pressure[0].fix()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(self.econ_hp2, tee=False)
        if not pyo.check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} econ_hp2 failed to initialize successfully."
            )

        propagate_state(self.econ_hp3.side_1_inlet, self.econ_hp2.side_1_outlet)
        self.econ_hp3.side_2_inlet.flow_mol_comp[0, "H2O"].fix(fg_rate * 0.0875)
        self.econ_hp3.side_2_inlet.flow_mol_comp[0, "CO2"].fix(fg_rate * 0.0408)
        self.econ_hp3.side_2_inlet.flow_mol_comp[0, "N2"].fix(fg_rate * 0.75)
        self.econ_hp3.side_2_inlet.flow_mol_comp[0, "O2"].fix(fg_rate * 0.1217)
        self.econ_hp3.side_2_inlet.temperature[0].fix(667.35)
        self.econ_hp3.side_2_inlet.pressure[0].fix(103421)
        self.econ_hp3.initialize(solver=solver, outlvl=outlvl, optarg=optarg)
        self.econ_hp3.deltaP_tube_eqn.deactivate()
        self.econ_hp3.deltaP_tube.fix(-8.96e05)
        self.econ_hp3.side_1_inlet.flow_mol[0].fix()
        self.econ_hp3.side_1_inlet.enth_mol[0].fix()
        self.econ_hp3.side_1_inlet.pressure[0].fix()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(self.econ_hp3, tee=False)
        if not pyo.check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} econ_hp3 failed to initialize successfully."
            )

        propagate_state(self.econ_hp4.side_1_inlet, self.econ_hp3.side_1_outlet)
        self.econ_hp4.side_2_inlet.flow_mol_comp[0, "H2O"].fix(fg_rate * 0.0875)
        self.econ_hp4.side_2_inlet.flow_mol_comp[0, "CO2"].fix(fg_rate * 0.0408)
        self.econ_hp4.side_2_inlet.flow_mol_comp[0, "N2"].fix(fg_rate * 0.75)
        self.econ_hp4.side_2_inlet.flow_mol_comp[0, "O2"].fix(fg_rate * 0.1217)
        self.econ_hp4.side_2_inlet.temperature[0].fix(624.35)
        self.econ_hp4.side_2_inlet.pressure[0].fix(103421)
        self.econ_hp4.initialize(solver=solver, outlvl=outlvl, optarg=optarg)
        self.econ_hp4.deltaP_tube_eqn.deactivate()
        self.econ_hp4.deltaP_tube.fix(-8.62e05)
        self.econ_hp4.side_1_inlet.flow_mol[0].fix()
        self.econ_hp4.side_1_inlet.enth_mol[0].fix()
        self.econ_hp4.side_1_inlet.pressure[0].fix()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(self.econ_hp4, tee=False)
        if not pyo.check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} econ_hp4 failed to initialize successfully."
            )

        propagate_state(self.econ_hp5.side_1_inlet, self.econ_hp4.side_1_outlet)
        self.econ_hp5.side_2_inlet.flow_mol_comp[0, "H2O"].fix(fg_rate * 0.0875)
        self.econ_hp5.side_2_inlet.flow_mol_comp[0, "CO2"].fix(fg_rate * 0.0408)
        self.econ_hp5.side_2_inlet.flow_mol_comp[0, "N2"].fix(fg_rate * 0.75)
        self.econ_hp5.side_2_inlet.flow_mol_comp[0, "O2"].fix(fg_rate * 0.1217)
        self.econ_hp5.side_2_inlet.temperature[0].fix(631.15)
        self.econ_hp5.side_2_inlet.pressure[0].fix(103421)
        self.econ_hp5.initialize(solver=solver, outlvl=outlvl, optarg=optarg)
        self.econ_hp5.deltaP_tube_eqn.deactivate()
        self.econ_hp5.deltaP_tube.fix(-8.27e05)
        self.econ_hp5.side_1_inlet.flow_mol[0].fix()
        self.econ_hp5.side_1_inlet.enth_mol[0].fix()
        self.econ_hp5.side_1_inlet.pressure[0].fix()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(self.econ_hp5, tee=False)
        if not pyo.check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} econ_hp5 failed to initialize successfully."
            )

        self.evap_hp.area.fix(8368.6)
        self.evap_hp.overall_heat_transfer_coefficient.fix(150)

        propagate_state(self.evap_hp.tube_inlet, self.econ_hp5.side_1_outlet)
        self.evap_hp.shell_inlet.flow_mol_comp[0, "H2O"].fix(fg_rate * 0.0875)
        self.evap_hp.shell_inlet.flow_mol_comp[0, "CO2"].fix(fg_rate * 0.0408)
        self.evap_hp.shell_inlet.flow_mol_comp[0, "N2"].fix(fg_rate * 0.75)
        self.evap_hp.shell_inlet.flow_mol_comp[0, "O2"].fix(fg_rate * 0.1217)
        self.evap_hp.shell_inlet.temperature[0].fix(729.15)
        self.evap_hp.shell_inlet.pressure[0].fix(103421)
        self.evap_hp.hp_sat_vap_eqn.deactivate()
        self.evap_hp.heat_transfer_equation.activate()
        self.evap_hp.initialize(solver=solver, outlvl=outlvl, optarg=optarg)
        self.evap_hp.hp_sat_vap_eqn.activate()
        self.evap_hp.heat_transfer_equation.deactivate()
        self.evap_hp.tube_inlet.flow_mol[0].fix()
        self.evap_hp.tube_inlet.enth_mol[0].fix()
        self.evap_hp.tube_inlet.pressure[0].fix()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(self.evap_hp, tee=False)
        if not pyo.check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} evap_hp failed to initialize successfully."
            )

        propagate_state(self.sh_hp1.side_1_inlet, self.evap_hp.tube_outlet)
        self.sh_hp1.side_2_inlet.flow_mol_comp[0, "H2O"].fix(fg_rate * 0.0875)
        self.sh_hp1.side_2_inlet.flow_mol_comp[0, "CO2"].fix(fg_rate * 0.0408)
        self.sh_hp1.side_2_inlet.flow_mol_comp[0, "N2"].fix(fg_rate * 0.75)
        self.sh_hp1.side_2_inlet.flow_mol_comp[0, "O2"].fix(fg_rate * 0.1217)
        self.sh_hp1.side_2_inlet.temperature[0].fix(760.35)
        self.sh_hp1.side_2_inlet.pressure[0].fix(103421)
        self.sh_hp1.initialize(solver=solver, outlvl=outlvl, optarg=optarg)
        self.sh_hp1.deltaP_tube_eqn.deactivate()
        self.sh_hp1.deltaP_tube.fix(-8.27e05)
        self.sh_hp1.side_1_inlet.flow_mol[0].fix()
        self.sh_hp1.side_1_inlet.enth_mol[0].fix()
        self.sh_hp1.side_1_inlet.pressure[0].fix()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(self.sh_hp1, tee=False)
        if not pyo.check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} sh_hp1 failed to initialize successfully."
            )

        propagate_state(self.sh_hp2.side_1_inlet, self.sh_hp1.side_1_outlet)
        self.sh_hp2.side_2_inlet.flow_mol_comp[0, "H2O"].fix(fg_rate * 0.0875)
        self.sh_hp2.side_2_inlet.flow_mol_comp[0, "CO2"].fix(fg_rate * 0.0408)
        self.sh_hp2.side_2_inlet.flow_mol_comp[0, "N2"].fix(fg_rate * 0.75)
        self.sh_hp2.side_2_inlet.flow_mol_comp[0, "O2"].fix(fg_rate * 0.1217)
        self.sh_hp2.side_2_inlet.temperature[0].fix(806.35)
        self.sh_hp2.side_2_inlet.pressure[0].fix(103421)
        self.sh_hp2.initialize(solver=solver, outlvl=outlvl, optarg=optarg)
        self.sh_hp2.deltaP_tube_eqn.deactivate()
        self.sh_hp2.deltaP_tube.fix(-7.93e05)
        self.sh_hp2.side_1_inlet.flow_mol[0].fix()
        self.sh_hp2.side_1_inlet.enth_mol[0].fix()
        self.sh_hp2.side_1_inlet.pressure[0].fix()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(self.sh_hp2, tee=False)
        if not pyo.check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} sh_hp2 failed to initialize successfully."
            )

        propagate_state(self.sh_hp3.side_1_inlet, self.sh_hp2.side_1_outlet)
        self.sh_hp3.side_2_inlet.flow_mol_comp[0, "H2O"].fix(fg_rate * 0.0875)
        self.sh_hp3.side_2_inlet.flow_mol_comp[0, "CO2"].fix(fg_rate * 0.0408)
        self.sh_hp3.side_2_inlet.flow_mol_comp[0, "N2"].fix(fg_rate * 0.75)
        self.sh_hp3.side_2_inlet.flow_mol_comp[0, "O2"].fix(fg_rate * 0.1217)
        self.sh_hp3.side_2_inlet.temperature[0].fix(835.35)
        self.sh_hp3.side_2_inlet.pressure[0].fix(103421)
        self.sh_hp3.initialize(solver=solver, outlvl=outlvl, optarg=optarg)
        self.sh_hp3.deltaP_tube_eqn.deactivate()
        self.sh_hp3.deltaP_tube.fix(-5.52e05)
        self.sh_hp3.side_1_inlet.flow_mol[0].fix()
        self.sh_hp3.side_1_inlet.enth_mol[0].fix()
        self.sh_hp3.side_1_inlet.pressure[0].fix()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(self.sh_hp3, tee=False)
        if not pyo.check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} sh_hp3 failed to initialize successfully."
            )

        propagate_state(self.sh_hp4.side_1_inlet, self.sh_hp3.side_1_outlet)
        self.sh_hp4.side_2_inlet.flow_mol_comp[0, "H2O"].fix(fg_rate * 0.0875)
        self.sh_hp4.side_2_inlet.flow_mol_comp[0, "CO2"].fix(fg_rate * 0.0408)
        self.sh_hp4.side_2_inlet.flow_mol_comp[0, "N2"].fix(fg_rate * 0.75)
        self.sh_hp4.side_2_inlet.flow_mol_comp[0, "O2"].fix(fg_rate * 0.1217)
        self.sh_hp4.side_2_inlet.temperature[0].fix(878.15)
        self.sh_hp4.side_2_inlet.pressure[0].fix(103421)
        self.sh_hp4.initialize(solver=solver, outlvl=outlvl, optarg=optarg)
        self.sh_hp1.deltaP_tube_eqn.deactivate()
        self.sh_hp1.deltaP_tube.fix(-5.45e05)
        self.sh_hp4.side_1_inlet.flow_mol[0].fix()
        self.sh_hp4.side_1_inlet.enth_mol[0].fix()
        self.sh_hp4.side_1_inlet.pressure[0].fix()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(self.sh_hp4, tee=False)
        if not pyo.check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} sh_hp4 failed to initialize successfully."
            )

        # unfix inlets that were fixed for initialization
        self.drum_lp.inlet.flow_mol[0].unfix()
        self.drum_lp.inlet.enth_mol[0].unfix()
        self.drum_lp.inlet.pressure[0].unfix()
        self.econ_ip1.side_1_inlet.flow_mol[0].unfix()
        self.econ_ip1.side_1_inlet.enth_mol[0].unfix()
        self.econ_ip1.side_1_inlet.pressure[0].unfix()
        self.evap_lp.tube_inlet.flow_mol[0].unfix()
        self.evap_lp.tube_inlet.enth_mol[0].unfix()
        self.evap_lp.tube_inlet.pressure[0].unfix()
        self.sh_lp.side_1_inlet.flow_mol[0].unfix()
        self.sh_lp.side_1_inlet.enth_mol[0].unfix()
        self.sh_lp.side_1_inlet.pressure[0].unfix()
        self.econ_ip2.side_1_inlet.flow_mol[0].unfix()
        self.econ_ip2.side_1_inlet.enth_mol[0].unfix()
        self.econ_ip2.side_1_inlet.pressure[0].unfix()
        self.evap_ip.tube_inlet.flow_mol[0].unfix()
        self.evap_ip.tube_inlet.pressure[0].unfix()
        self.evap_ip.tube_inlet.enth_mol[0].unfix()
        self.sh_ip1.side_1_inlet.flow_mol[0].unfix()
        self.sh_ip1.side_1_inlet.enth_mol[0].unfix()
        self.sh_ip1.side_1_inlet.pressure[0].unfix()
        self.mixer_ip1.sh_ip1.flow_mol[0].unfix()
        self.mixer_ip1.sh_ip1.enth_mol[0].unfix()
        self.mixer_ip1.sh_ip1.pressure[0].unfix()
        self.mixer_ip1.Cold_reheat.flow_mol[0].unfix()
        self.mixer_ip1.Cold_reheat.enth_mol[0].unfix()
        self.mixer_ip1.Cold_reheat.pressure[0].unfix()
        self.sh_ip2.side_1_inlet.flow_mol[0].unfix()
        self.sh_ip2.side_1_inlet.enth_mol[0].unfix()
        self.sh_ip2.side_1_inlet.pressure[0].unfix()
        self.sh_ip3.side_1_inlet.flow_mol[0].unfix()
        self.sh_ip3.side_1_inlet.enth_mol[0].unfix()
        self.sh_ip3.side_1_inlet.pressure[0].unfix()
        self.econ_hp1.side_1_inlet.flow_mol[0].unfix()
        self.econ_hp1.side_1_inlet.enth_mol[0].unfix()
        self.econ_hp1.side_1_inlet.pressure[0].unfix()
        self.econ_hp2.side_1_inlet.flow_mol[0].unfix()
        self.econ_hp2.side_1_inlet.enth_mol[0].unfix()
        self.econ_hp2.side_1_inlet.pressure[0].unfix()
        self.econ_hp3.side_1_inlet.flow_mol[0].unfix()
        self.econ_hp3.side_1_inlet.enth_mol[0].unfix()
        self.econ_hp3.side_1_inlet.pressure[0].unfix()
        self.econ_hp4.side_1_inlet.flow_mol[0].unfix()
        self.econ_hp4.side_1_inlet.enth_mol[0].unfix()
        self.econ_hp4.side_1_inlet.pressure[0].unfix()
        self.econ_hp5.side_1_inlet.flow_mol[0].unfix()
        self.econ_hp5.side_1_inlet.enth_mol[0].unfix()
        self.econ_hp5.side_1_inlet.pressure[0].unfix()
        self.evap_hp.tube_inlet.flow_mol[0].unfix()
        self.evap_hp.tube_inlet.enth_mol[0].unfix()
        self.evap_hp.tube_inlet.pressure[0].unfix()
        self.sh_hp1.side_1_inlet.flow_mol[0].unfix()
        self.sh_hp1.side_1_inlet.enth_mol[0].unfix()
        self.sh_hp1.side_1_inlet.pressure[0].unfix()
        self.sh_hp2.side_1_inlet.flow_mol[0].unfix()
        self.sh_hp2.side_1_inlet.enth_mol[0].unfix()
        self.sh_hp2.side_1_inlet.pressure[0].unfix()
        self.sh_hp3.side_1_inlet.flow_mol[0].unfix()
        self.sh_hp3.side_1_inlet.enth_mol[0].unfix()
        self.sh_hp3.side_1_inlet.pressure[0].unfix()
        self.sh_hp4.side_1_inlet.flow_mol[0].unfix()
        self.sh_hp4.side_1_inlet.enth_mol[0].unfix()
        self.sh_hp4.side_1_inlet.pressure[0].unfix()
        self.sh_ip3.side_2_inlet.flow_mol_comp[:, :].unfix()
        self.sh_ip3.side_2_inlet.pressure[0].unfix()
        self.sh_ip3.side_2_inlet.temperature[0].unfix()
        self.sh_hp3.side_2_inlet.flow_mol_comp[:, :].unfix()
        self.sh_hp3.side_2_inlet.pressure[0].unfix()
        self.sh_hp3.side_2_inlet.temperature[0].unfix()
        self.sh_hp2.side_2_inlet.flow_mol_comp[:, :].unfix()
        self.sh_hp2.side_2_inlet.pressure[0].unfix()
        self.sh_hp2.side_2_inlet.temperature[0].unfix()
        self.sh_ip2.side_2_inlet.flow_mol_comp[:, :].unfix()
        self.sh_ip2.side_2_inlet.pressure[0].unfix()
        self.sh_ip2.side_2_inlet.temperature[0].unfix()
        self.sh_hp1.side_2_inlet.flow_mol_comp[:, :].unfix()
        self.sh_hp1.side_2_inlet.pressure[0].unfix()
        self.sh_hp1.side_2_inlet.temperature[0].unfix()
        self.evap_hp.shell_inlet.flow_mol_comp[:, :].unfix()
        self.evap_hp.shell_inlet.pressure[0].unfix()
        self.evap_hp.shell_inlet.temperature[0].unfix()
        self.econ_hp5.side_2_inlet.flow_mol_comp[:, :].unfix()
        self.econ_hp5.side_2_inlet.pressure[0].unfix()
        self.econ_hp5.side_2_inlet.temperature[0].unfix()
        self.sh_ip1.side_2_inlet.flow_mol_comp[:, :].unfix()
        self.sh_ip1.side_2_inlet.pressure[0].unfix()
        self.sh_ip1.side_2_inlet.temperature[0].unfix()
        self.econ_hp4.side_2_inlet.flow_mol_comp[:, :].unfix()
        self.econ_hp4.side_2_inlet.pressure[0].unfix()
        self.econ_hp4.side_2_inlet.temperature[0].unfix()
        self.econ_hp3.side_2_inlet.flow_mol_comp[:, :].unfix()
        self.econ_hp3.side_2_inlet.pressure[0].unfix()
        self.econ_hp3.side_2_inlet.temperature[0].unfix()
        self.split_fg_lp.inlet.flow_mol_comp[:, :].unfix()
        self.split_fg_lp.inlet.pressure[0].unfix()
        self.split_fg_lp.inlet.temperature[0].unfix()
        self.sh_lp.side_2_inlet.flow_mol_comp[:, :].unfix()
        self.sh_lp.side_2_inlet.pressure[0].unfix()
        self.sh_lp.side_2_inlet.temperature[0].unfix()
        self.split_fg_lp.inlet.flow_mol_comp[:, :].unfix()
        self.split_fg_lp.inlet.pressure[0].unfix()
        self.split_fg_lp.inlet.temperature[0].unfix()
        self.evap_ip.shell_inlet.flow_mol_comp[:, :].unfix()
        self.evap_ip.shell_inlet.pressure[0].unfix()
        self.evap_ip.shell_inlet.temperature[0].unfix()
        self.econ_ip2.side_2_inlet.flow_mol_comp[:, :].unfix()
        self.econ_ip2.side_2_inlet.pressure[0].unfix()
        self.econ_ip2.side_2_inlet.temperature[0].unfix()
        self.econ_hp2.side_2_inlet.flow_mol_comp[:, :].unfix()
        self.econ_hp2.side_2_inlet.pressure[0].unfix()
        self.econ_hp2.side_2_inlet.temperature[0].unfix()
        self.econ_ip1.side_2_inlet.flow_mol_comp[:, :].unfix()
        self.econ_ip1.side_2_inlet.pressure[0].unfix()
        self.econ_ip1.side_2_inlet.temperature[0].unfix()
        self.econ_hp1.side_2_inlet.flow_mol_comp[:, :].unfix()
        self.econ_hp1.side_2_inlet.pressure[0].unfix()
        self.econ_hp1.side_2_inlet.temperature[0].unfix()
        self.evap_lp.shell_inlet.flow_mol_comp[:, :].unfix()
        self.evap_lp.shell_inlet.pressure[0].unfix()
        self.evap_lp.shell_inlet.temperature[0].unfix()
        self.econ_lp.side_2_inlet.flow_mol_comp[:, :].unfix()
        self.econ_lp.side_2_inlet.pressure[0].unfix()
        self.econ_lp.side_2_inlet.temperature[0].unfix()

        if save_to is not None:
            iutil.to_json(self, fname=save_to)
            init_log.info_low(f"Initialization saved to {save_to}")
        init_log.info("High pressure system initialization - Completed")

    def _stream_tags(self):
        tag_stm = iutil.ModelTagGroup()
        tag_gas = iutil.ModelTagGroup()
        self.tags_steam_streams = tag_stm
        self.tags_flue_gas_streams = tag_gas
        stream_states = ta.stream_states_dict(
            ta.arcs_to_stream_dict(
                self,
                descend_into=False,
                additional={  # these are streams in or out without an arc
                    "lp01": self.econ_lp.side_1_inlet,
                    "lp03": self.mixer1.Preheater,
                    "lp12": self.mixer_soec.soec_makeup,
                    "ip04": self.splitter_ip1.toNGPH,
                    "lp11": self.sh_lp.side_1_outlet,
                    "hp11": self.sh_hp4.side_1_outlet,
                    "ip10": self.sh_ip3.side_1_outlet,
                    "ip11": self.splitter_ip2.inlet,
                    "ip15": self.splitter_ip2.Cold_reheat,
                    "ip13": self.splitter_ip2.toReclaimer,
                    "ip12": self.splitter_ip2.toEjector,
                    "ip14": self.splitter_ip2.toDryer,
                    "g30": self.econ_lp.side_2_outlet,
                    "g19": self.econ_hp3.side_2_outlet,
                    "g08": self.sh_hp4.side_2_inlet,
                },
            )
        )
        for i, s in stream_states.items():  # create the tags for steam quantities
            if isinstance(s, iapws95.Iapws95StateBlockData):
                tag_group = tag_stm
                is_steam = True
            else:
                tag_group = tag_gas
                is_steam = False
            tag_group[f"{i}_F"] = iutil.ModelTag(
                doc=f"{i}: mass flow",
                expr=s.flow_mass,
                format_string="{:.3f}",
                display_units=pyo.units.kg / pyo.units.s,
            )
            tag_group[f"{i}_Fmol"] = iutil.ModelTag(
                doc=f"{i}: mole flow",
                expr=s.flow_mol,
                format_string="{:.3f}",
                display_units=pyo.units.kmol / pyo.units.s,
            )
            tag_group[f"{i}_Fvol"] = iutil.ModelTag(
                doc=f"{i}: volumetric flow",
                expr=s.flow_vol,
                format_string="{:.3f}",
                display_units=pyo.units.m ** 3 / pyo.units.s,
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
            if is_steam:
                tag_group[f"{i}_X"] = iutil.ModelTag(
                    doc=f"{i}: vapor fraction",
                    expr=s.phase_frac["Vap"],
                    format_string="{:.3f}",
                    display_units=None,
                )
                tag_group[f"{i}_H"] = iutil.ModelTag(
                    doc=f"{i}: molar enthalpy",
                    expr=s.enth_mol,
                    format_string="{:.3f}",
                    display_units=pyo.units.kJ / pyo.units.mol,
                )
            else:
                for c in s.mole_frac_comp:
                    tag_group[f"{i}_y{c}"] = iutil.ModelTag(
                        doc=f"{i}: mole percent {c}",
                        expr=s.mole_frac_comp[c] * 100,
                        format_string="{:.3f}",
                        display_units="%",
                    )

    def write_pfd(self, fname=None):
        """Add model results to the flowsheet template.  If fname is specified,
        this saves the resulting svg to a file.  If fname is not specified, it
        returns the svg string.

        Args:
            fname: Name of file to save svg.  If None, return the svg string

        Returns: (None or Str)
        """
        infilename = os.path.join(this_file_dir(), "hrsg_template.svg")
        with open(infilename, "r") as f:
            s = svg_tag(svg=f, tag_group=self.tags_steam_streams)
        s = svg_tag(svg=s, tag_group=self.tags_flue_gas_streams, outfile=fname)
        if fname is None:
            return s

    @staticmethod
    def _stream_col_gen(tag_group):
        for tag in tag_group.values():
            spltstr = tag.doc.split(":")
            stream = spltstr[0].strip()
            col = f"{spltstr[1].strip()} ({tag.get_unit_str()})"
            yield tag, stream, col

    @staticmethod
    def _stream_table(tag_group):
        rows = set()
        cols = set()
        tags = []
        for tag, stream, col in HrsgFlowsheetData._stream_col_gen(tag_group):
            rows.add(stream)
            cols.add(col)
            tags.append((tag, stream, col))
        df = pd.DataFrame(index=sorted(rows), columns=sorted(cols))
        for tag, stream, col in tags:
            df.at[stream, col] = tag.get_display_value()
        return df

    def steam_streams_dataframe(self):
        return self._stream_table(self.tags_steam_streams)

    def flue_gas_streams_dataframe(self):
        return self._stream_table(self.tags_flue_gas_streams)
