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
import ngcc
import soec_design
import idaes.core.util as iutil
from idaes.core.util.initialization import propagate_state

@declare_process_block_class(
    "NgccSoecDesignFlowsheet",
    doc=(
        "This flowsheet combines an NGCC and SOEC that uses steam from the HRSG."
        "This is a desing point only model."
    ),
)
class NgccFlowsheetData(FlowsheetBlockData):
    def build(self):
        super().build()
        self._add_flowsheets()
        self._add_arcs()
        self._add_constraints()

    def _add_flowsheets(self):
        self.ngcc = ngcc.NgccFlowsheet(
            default={
                "dynamic":self.config.dynamic,
                "time":self.time,
                "time_units": self.config.time_units
            }
        )
        self.soec = soec_design.SoecDesignFlowsheet(
            default={
                "dynamic":self.config.dynamic,
                "time":self.time,
                "time_units": self.config.time_units
            }
        )

    def _add_arcs(self):
        self.c01 = Arc(
            source=self.ngcc.st.steam_turbine_lp_split.soec,
            destination=self.soec.feed_hx01.tube_inlet,
        )
        self.c02 = Arc(
            source=self.soec.makeup_mix.outlet,
            destination=self.ngcc.hrsg.mixer_soec.soec_makeup,
        )
        pyo.TransformationFactory("network.expand_arcs").apply_to(self)

    def _add_constraints(self):
        @self.Constraint(self.config.time)
        def sweep_pressure_eqn(b, t):
            return (
                b.soec.sweep_compressor.control_volume.properties_out[t].pressure ==
                b.soec.feed_hx01.tube.properties_in[t].pressure
            )

    def initialize(
        self,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
        load_from="ngcc_soec_design_init.json.gz",
        save_to="ngcc_soec_design_init.json.gz",
    ):

        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="flowsheet")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="flowsheet")

        if load_from is not None:
            if os.path.exists(load_from):
                init_log.info_high(f"NGCC/SOEC design load initial from {load_from}")
                # here suffix=False avoids loading scaling factors
                iutil.from_json(
                    self, fname=load_from, wts=iutil.StoreSpec(suffix=False)
                )
                return
        solver_obj = iutil.get_solver(solver, optarg)
        self.ngcc.initialize(
            outlvl=outlvl,
            solver=solver,
            optarg=optarg,
            load_from="ngcc_init.json.gz",
            save_to="ngcc_init.json.gz",
        )
        self.soec.initialize(
            outlvl=outlvl,
            solver=solver,
            optarg=optarg,
            load_from="soec_design_init.json.gz",
            save_to="soec_desgin_init.json.gz",
        )
        self.soec.sweep_compressor.control_volume.properties_out[:].pressure.unfix()
        self.c01.expanded_block.deactivate()
        self.c02.expanded_block.deactivate()
        self.ngcc.hrsg.mixer_soec.soec_makeup.fix()
        self.soec.feed_hx01.tube_inlet.fix()
        self.soec.water_pump.inlet.flow_mol.fix()
        solver_obj.solve(self, tee=True)

        self.ngcc.st.steam_turbine_lp_split.split_fraction[0, "soec"].unfix()
        self.ngcc.st.steam_turbine_lp_split.soec.flow_mol.fix(1200)
        propagate_state(self.c02)
        solver_obj.solve(self, tee=True)

        propagate_state(self.c01)
        propagate_state(self.c02)
        solver_obj.solve(self, tee=True)

        self.ngcc.hrsg.mixer_soec.soec_makeup.unfix()
        self.soec.feed_hx01.tube_inlet.unfix()
        self.soec.feed_hx01.tube_inlet.flow_mol.unfix()
        self.soec.water_pump.inlet.flow_mol.fix(1200)
        self.c01.expanded_block.activate()
        self.c02.expanded_block.activate()
        solver_obj.solve(self, tee=True)

        # switch to 97% capture and set GT power to 477 for full load
        self.ngcc.net_power_mw.unfix()
        self.ngcc.gt.gt_power.fix(-477e6)
        self.ngcc.cap_specific_reboiler_duty.fix(2.4e6)
        self.ngcc.cap_fraction.fix(0.97)
        self.soec.cmp01.ratioP.fix(2.07)
        self.soec.cmp02.ratioP.fix(2.07)
        self.soec.cmp03.ratioP.fix(2.07)
        self.soec.cmp04.ratioP.fix(2.07)
        self.soec.cmp05.ratioP.fix(2.07)
        self.soec.cmp06.ratioP.fix(2.07)
        self.soec.cmp06.ratioP.unfix()
        self.soec.cmp06.control_volume.properties_out[:].pressure.fix(320*pyo.units.bar)
        solver_obj.solve(self, tee=True)

        if save_to is not None:
            iutil.to_json(self)
            if save_to is not None:
                iutil.to_json(self, fname=save_to)
                init_log.info_high(f"Initialization saved to {save_to}")
