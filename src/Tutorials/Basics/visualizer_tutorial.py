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
IDAES Visualizer Tutorial script

This runs the same code as the Jupyter notebook by the same name.
See that notebook for details.
"""
# stdlib
import sys
# Pyomo
from pyomo.environ import (Constraint,
                           Var,
                           ConcreteModel,
                           Expression,
                           Objective,
                           SolverFactory,
                           TransformationFactory,
                           value)
from pyomo.network import Arc, SequentialDecomposition
# IDAES
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import (PressureChanger,
                                        Mixer,
                                        Separator as Splitter,
                                        Heater,
                                        StoichiometricReactor)
from idaes.generic_models.unit_models import Flash
from idaes.generic_models.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.core.util.model_statistics import degrees_of_freedom


__author__ = "Dan Gunter"  # GH: dangunter, email: dkgunter@lbl.gov


def create_flowsheet() -> FlowsheetBlock:
    """Create an IDAES flowsheet.

    {description of the flowsheet created}
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    return m.fs


def initialize_flowsheet(fs: FlowsheetBlock):
    pass


def main():
    print("This script doesn't do anything yet!")
    return 0


if __name__ == "__main__":
    sys.exit(main())

