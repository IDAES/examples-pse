##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
Flowsheets for HDA with Flash and HDA with Distillation.
"""

# Import Pyomo libraries
from pyomo.environ import (Constraint,
                           Var,
                           ConcreteModel,
                           SolverFactory,
                           TransformationFactory,
                           units as pyunits)
from pyomo.network import Arc, SequentialDecomposition

# Import IDAES core libraries
from idaes.core import FlowsheetBlock

# Import IDAES generic unit models
from idaes.models.unit_models import (PressureChanger,
                                      Mixer,
                                      Separator as Splitter,
                                      Heater,
                                      StoichiometricReactor,
                                      CSTR,
                                      Flash,
                                      Translator)
from idaes.models.unit_models.pressure_changer import \
    ThermodynamicAssumption

# Import IDAES distillation unit models
from idaes.models_extra.column_models.tray_column import TrayColumn
from idaes.models_extra.column_models.condenser import \
    CondenserType, TemperatureSpec

# Import IDAES core utilities
from idaes.core.util.initialization import propagate_state
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
from pyomo.util.check_units import assert_units_consistent

# Import idaes logger to set output levels
import idaes.logger as idaeslog


def hda_with_flash(tee=True):
    if tee is True:
        outlvl = idaeslog.INFO
    else:
        outlvl = idaeslog.ERROR

    # Import thermodynamic and reaction property packages
    import hda_ideal_VLE as thermo_props
    import hda_reaction as reaction_props

    # build flowsheet
    print('Building flowsheet...')
    print()

    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.thermo_params = thermo_props.HDAParameterBlock()
    m.fs.reaction_params = reaction_props.HDAReactionParameterBlock(
        default={"property_package": m.fs.thermo_params})

    m.fs.M101 = Mixer(default={"property_package": m.fs.thermo_params,
                               "inlet_list": ["toluene_feed", "hydrogen_feed",
                                              "vapor_recycle"]})
    m.fs.H101 = Heater(default={"property_package": m.fs.thermo_params,
                                "has_pressure_change": False,
                                "has_phase_equilibrium": True})
    m.fs.R101 = StoichiometricReactor(
                default={"property_package": m.fs.thermo_params,
                         "reaction_package": m.fs.reaction_params,
                         "has_heat_of_reaction": True,
                         "has_heat_transfer": True,
                         "has_pressure_change": False})
    m.fs.F101 = Flash(default={"property_package": m.fs.thermo_params,
                               "has_heat_transfer": True,
                               "has_pressure_change": True})
    m.fs.S101 = Splitter(default={"property_package": m.fs.thermo_params,
                                  "ideal_separation": False,
                                  "outlet_list": ["purge", "recycle"]})
    m.fs.C101 = PressureChanger(default={
                "property_package": m.fs.thermo_params,
                "compressor": True,
                "thermodynamic_assumption":
                ThermodynamicAssumption.isothermal})
    m.fs.F102 = Flash(default={"property_package": m.fs.thermo_params,
                               "has_heat_transfer": True,
                               "has_pressure_change": True})

    m.fs.s03 = Arc(source=m.fs.M101.outlet, destination=m.fs.H101.inlet)
    m.fs.s04 = Arc(source=m.fs.H101.outlet, destination=m.fs.R101.inlet)
    m.fs.s05 = Arc(source=m.fs.R101.outlet, destination=m.fs.F101.inlet)
    m.fs.s06 = Arc(source=m.fs.F101.vap_outlet, destination=m.fs.S101.inlet)
    m.fs.s08 = Arc(source=m.fs.S101.recycle, destination=m.fs.C101.inlet)
    m.fs.s09 = Arc(source=m.fs.C101.outlet,
                   destination=m.fs.M101.vapor_recycle)
    m.fs.s10 = Arc(source=m.fs.F101.liq_outlet, destination=m.fs.F102.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    # set inputs
    print('Setting inputs...')
    print()

    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "benzene"].fix(1e-5*pyunits.mol/pyunits.s)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "toluene"].fix(1e-5*pyunits.mol/pyunits.s)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "hydrogen"].fix(1e-5*pyunits.mol/pyunits.s)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "methane"].fix(1e-5*pyunits.mol/pyunits.s)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "benzene"].fix(1e-5*pyunits.mol/pyunits.s)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "toluene"].fix(0.30*pyunits.mol/pyunits.s)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "hydrogen"].fix(1e-5*pyunits.mol/pyunits.s)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "methane"].fix(1e-5*pyunits.mol/pyunits.s)
    m.fs.M101.toluene_feed.temperature.fix(303.2*pyunits.K)
    m.fs.M101.toluene_feed.pressure.fix(350000*pyunits.Pa)

    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "benzene"].fix(1e-5*pyunits.mol/pyunits.s)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "toluene"].fix(1e-5*pyunits.mol/pyunits.s)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "hydrogen"].fix(0.30*pyunits.mol/pyunits.s)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "methane"].fix(0.02*pyunits.mol/pyunits.s)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "benzene"].fix(1e-5*pyunits.mol/pyunits.s)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "toluene"].fix(1e-5*pyunits.mol/pyunits.s)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "hydrogen"].fix(1e-5*pyunits.mol/pyunits.s)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "methane"].fix(1e-5*pyunits.mol/pyunits.s)
    m.fs.M101.hydrogen_feed.temperature.fix(303.2*pyunits.K)
    m.fs.M101.hydrogen_feed.pressure.fix(350000*pyunits.Pa)

    m.fs.H101.outlet.temperature.fix(600*pyunits.K)

    m.fs.R101.conversion = Var(initialize=0.75, bounds=(0, 1))

    m.fs.R101.conv_constraint = Constraint(
        expr=m.fs.R101.conversion*m.fs.R101.inlet.
        flow_mol_phase_comp[0, "Vap", "toluene"] ==
        (m.fs.R101.inlet.flow_mol_phase_comp[0, "Vap", "toluene"] -
         m.fs.R101.outlet.flow_mol_phase_comp[0, "Vap", "toluene"]))

    m.fs.R101.conversion.fix(0.75*pyunits.dimensionless)
    m.fs.R101.heat_duty.fix(0*pyunits.W)

    m.fs.F101.vap_outlet.temperature.fix(325.0*pyunits.K)
    m.fs.F101.deltaP.fix(0*pyunits.Pa)

    m.fs.F102.vap_outlet.temperature.fix(375*pyunits.K)
    m.fs.F102.deltaP.fix(-200000*pyunits.Pa)

    m.fs.S101.split_fraction[0, "purge"].fix(0.2)
    m.fs.C101.outlet.pressure.fix(350000*pyunits.Pa)

    # initialize flowsheet
    print('Initializing flowsheet...')
    print()

    seq = SequentialDecomposition()
    seq.options.select_tear_method = "heuristic"
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = 5
    print('Limiting Wegstein tear to 5 iterations to obtain initial solution,'
          ' if not converged IPOPT will pick up and continue.')
    print()

    G = seq.create_graph(m)
    heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
    order = seq.calculation_order(G)
    if tee is True:
        for o in heuristic_tear_set:
            print(o.name)
        for o in order:
            print(o[0].name)

    tear_guesses = {
        "flow_mol_phase_comp": {
                (0, "Vap", "benzene"): 1e-5,
                (0, "Vap", "toluene"): 1e-5,
                (0, "Vap", "hydrogen"): 0.30,
                (0, "Vap", "methane"): 0.02,
                (0, "Liq", "benzene"): 1e-5,
                (0, "Liq", "toluene"): 0.30,
                (0, "Liq", "hydrogen"): 1e-5,
                (0, "Liq", "methane"): 1e-5},
        "temperature": {0: 303},
        "pressure": {0: 350000}}

    seq.set_guesses_for(m.fs.H101.inlet, tear_guesses)

    def function(unit):
        unit.initialize(outlvl=outlvl)

    seq.run(m, function)

    # solve model
    print('Solving flowsheet...')
    print()

    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6, 'max_iter': 5000}
    results = solver.solve(m, tee=tee)
    from pyomo.environ import TerminationCondition
    assert results.solver.termination_condition == TerminationCondition.optimal
    assert_units_consistent(m)

    print('Complete.')
    print()

    return m


def hda_with_distillation(tee=True):
    if tee is True:
        outlvl = idaeslog.INFO
    else:
        outlvl = idaeslog.ERROR

    # Import thermodynamic and reaction property packages
    import hda_reaction as reaction_props
    from idaes.models.properties.activity_coeff_models.\
        BTX_activity_coeff_VLE import BTXParameterBlock
    from hda_ideal_VLE import HDAParameterBlock

    # build flowsheet
    print('Building flowsheet...')
    print()

    n = ConcreteModel()
    n.fs = FlowsheetBlock(default={"dynamic": False})

    n.fs.BTHM_params = HDAParameterBlock()
    n.fs.BT_params = BTXParameterBlock(default={
        "valid_phase": ('Liq', 'Vap'),
        "activity_coeff_model": "Ideal"
        })
    n.fs.reaction_params = reaction_props.HDAReactionParameterBlock(
        default={"property_package": n.fs.BTHM_params})

    n.fs.M101 = Mixer(default={"property_package": n.fs.BTHM_params,
                               "inlet_list": ["toluene_feed", "hydrogen_feed",
                                              "vapor_recycle"]})
    n.fs.H101 = Heater(default={"property_package": n.fs.BTHM_params,
                                "has_phase_equilibrium": True})
    n.fs.R101 = CSTR(
            default={"property_package": n.fs.BTHM_params,
                     "reaction_package": n.fs.reaction_params,
                     "has_heat_of_reaction": True,
                     "has_heat_transfer": True})
    n.fs.F101 = Flash(default={"property_package": n.fs.BTHM_params,
                               "has_heat_transfer": True,
                               "has_pressure_change": True})
    n.fs.S101 = Splitter(default={"property_package": n.fs.BTHM_params,
                                  "outlet_list": ["purge", "recycle"]})
    n.fs.C101 = PressureChanger(default={
                "property_package": n.fs.BTHM_params,
                "compressor": True,
                "thermodynamic_assumption":
                ThermodynamicAssumption.isothermal})
    n.fs.translator = Translator(default={
        "inlet_property_package": n.fs.BTHM_params,
        "outlet_property_package": n.fs.BT_params
        })

    # Add constraint: Total flow = benzene flow + toluene flow (molar)
    n.fs.translator.eq_total_flow = Constraint(
        expr=n.fs.translator.outlet.flow_mol[0] ==
        n.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "benzene"] +
        n.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "toluene"])

    # Add constraint: Outlet temperature = Inlet temperature
    n.fs.translator.eq_temperature = Constraint(
        expr=n.fs.translator.outlet.temperature[0] ==
        n.fs.translator.inlet.temperature[0])

    n.fs.translator.eq_pressure = Constraint(
        expr=n.fs.translator.outlet.pressure[0] ==
        n.fs.translator.inlet.pressure[0])

    # Add constraint: Benzene mole fraction definition
    n.fs.translator.eq_mole_frac_benzene = Constraint(
        expr=n.fs.translator.outlet.mole_frac_comp[0, "benzene"] ==
        n.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "benzene"] /
        (n.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "benzene"] +
         n.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "toluene"]))

    # Add constraint: Toluene mole fraction definition
    n.fs.translator.eq_mole_frac_toluene = Constraint(
        expr=n.fs.translator.outlet.mole_frac_comp[0, "toluene"] ==
        n.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "toluene"] /
        (n.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "benzene"] +
         n.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "toluene"]))

    n.fs.H102 = Heater(default={"property_package": n.fs.BT_params,
                                "has_pressure_change": True,
                                "has_phase_equilibrium": True})

    n.fs.s03 = Arc(source=n.fs.M101.outlet, destination=n.fs.H101.inlet)
    n.fs.s04 = Arc(source=n.fs.H101.outlet, destination=n.fs.R101.inlet)
    n.fs.s05 = Arc(source=n.fs.R101.outlet, destination=n.fs.F101.inlet)
    n.fs.s06 = Arc(source=n.fs.F101.vap_outlet, destination=n.fs.S101.inlet)
    n.fs.s08 = Arc(source=n.fs.S101.recycle, destination=n.fs.C101.inlet)
    n.fs.s09 = Arc(source=n.fs.C101.outlet,
                   destination=n.fs.M101.vapor_recycle)
    n.fs.s10a = Arc(source=n.fs.F101.liq_outlet,
                    destination=n.fs.translator.inlet)
    n.fs.s10b = Arc(source=n.fs.translator.outlet, destination=n.fs.H102.inlet)

    TransformationFactory("network.expand_arcs").apply_to(n)

    # set inputs
    print('Setting inputs...')
    print()

    n.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "benzene"].fix(1e-5*pyunits.mol/pyunits.s)
    n.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "toluene"].fix(1e-5*pyunits.mol/pyunits.s)
    n.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "hydrogen"].fix(1e-5*pyunits.mol/pyunits.s)
    n.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "methane"].fix(1e-5*pyunits.mol/pyunits.s)
    n.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "benzene"].fix(1e-5*pyunits.mol/pyunits.s)
    n.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "toluene"].fix(0.30*pyunits.mol/pyunits.s)
    n.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "hydrogen"].fix(1e-5*pyunits.mol/pyunits.s)
    n.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "methane"].fix(1e-5*pyunits.mol/pyunits.s)
    n.fs.M101.toluene_feed.temperature.fix(303.2*pyunits.K)
    n.fs.M101.toluene_feed.pressure.fix(350000*pyunits.Pa)

    n.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "benzene"].fix(1e-5*pyunits.mol/pyunits.s)
    n.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "toluene"].fix(1e-5*pyunits.mol/pyunits.s)
    n.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "hydrogen"].fix(0.30*pyunits.mol/pyunits.s)
    n.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "methane"].fix(0.02*pyunits.mol/pyunits.s)
    n.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "benzene"].fix(1e-5*pyunits.mol/pyunits.s)
    n.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "toluene"].fix(1e-5*pyunits.mol/pyunits.s)
    n.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "hydrogen"].fix(1e-5*pyunits.mol/pyunits.s)
    n.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "methane"].fix(1e-5*pyunits.mol/pyunits.s)
    n.fs.M101.hydrogen_feed.temperature.fix(303.2*pyunits.K)
    n.fs.M101.hydrogen_feed.pressure.fix(350000*pyunits.Pa)

    n.fs.H101.outlet.temperature.fix(600*pyunits.K)

    n.fs.R101.conversion = Var(initialize=0.75, bounds=(0, 1))

    n.fs.R101.conv_constraint = Constraint(
        expr=n.fs.R101.conversion*n.fs.R101.inlet.
        flow_mol_phase_comp[0, "Vap", "toluene"] ==
        (n.fs.R101.inlet.flow_mol_phase_comp[0, "Vap", "toluene"] -
         n.fs.R101.outlet.flow_mol_phase_comp[0, "Vap", "toluene"]))

    n.fs.R101.conversion.fix(0.75*pyunits.dimensionless)
    n.fs.R101.heat_duty.fix(0*pyunits.W)

    n.fs.F101.vap_outlet.temperature.fix(325.0*pyunits.K)
    n.fs.F101.deltaP.fix(0*pyunits.Pa)

    n.fs.S101.split_fraction[0, "purge"].fix(0.2)
    n.fs.C101.outlet.pressure.fix(350000*pyunits.Pa)

    n.fs.H102.outlet.temperature.fix(375*pyunits.K)
    n.fs.H102.deltaP.fix(-200000*pyunits.Pa)

    # set scaling factors
    # Set scaling factors for heat duty, reaction extent and volume
    iscale.set_scaling_factor(n.fs.H101.control_volume.heat, 1e-2)
    iscale.set_scaling_factor(n.fs.R101.control_volume.heat, 1e-2)
    iscale.set_scaling_factor(n.fs.R101.control_volume.rate_reaction_extent, 1)
    iscale.set_scaling_factor(n.fs.R101.control_volume.volume, 1)
    iscale.set_scaling_factor(n.fs.F101.control_volume.heat, 1e-2)
    iscale.set_scaling_factor(n.fs.H102.control_volume.heat, 1e-2)

    # Set the scaling factors for the remaining variables and all constraints
    iscale.calculate_scaling_factors(n.fs.H101)
    iscale.calculate_scaling_factors(n.fs.R101)
    iscale.calculate_scaling_factors(n.fs.F101)
    iscale.calculate_scaling_factors(n.fs.H102)

    # initialize flowsheet
    print('Initializing flowsheet...')
    print()

    seq = SequentialDecomposition()
    seq.options.select_tear_method = "heuristic"
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = 3
    print('Limiting Wegstein tear to 3 iterations to obtain initial solution,'
          ' if not converged IPOPT will pick up and continue.')
    print()

    G = seq.create_graph(n)
    heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
    order = seq.calculation_order(G)
    if tee is True:
        for o in heuristic_tear_set:
            print(o.name)
        for o in order:
            print(o[0].name)

    tear_guesses = {
        "flow_mol_phase_comp": {
                (0, "Vap", "benzene"): 1e-5,
                (0, "Vap", "toluene"): 1e-5,
                (0, "Vap", "hydrogen"): 0.30,
                (0, "Vap", "methane"): 0.02,
                (0, "Liq", "benzene"): 1e-5,
                (0, "Liq", "toluene"): 0.30,
                (0, "Liq", "hydrogen"): 1e-5,
                (0, "Liq", "methane"): 1e-5},
        "temperature": {0: 303},
        "pressure": {0: 350000}}

    seq.set_guesses_for(n.fs.H101.inlet, tear_guesses)

    def function(unit):
        unit.initialize(outlvl=outlvl)

    seq.run(n, function)

    # solve model
    print('Solving flowsheet...')
    print()

    solver = get_solver()
    results = solver.solve(n, tee=tee)
    from pyomo.environ import TerminationCondition
    assert results.solver.termination_condition == TerminationCondition.optimal

    # add and initialize distilation column, and resolve
    print('Adding distillation column and resolving flowsheet...')
    print()

    n.fs.D101 = TrayColumn(default={
                        "number_of_trays": 10,
                        "feed_tray_location": 5,
                        "condenser_type": CondenserType.totalCondenser,
                        "condenser_temperature_spec":
                        TemperatureSpec.atBubblePoint,
                        "property_package": n.fs.BT_params})

    n.fs.s11 = Arc(source=n.fs.H102.outlet, destination=n.fs.D101.feed)

    TransformationFactory("network.expand_arcs").apply_to(n)

    propagate_state(n.fs.s11)

    n.fs.D101.condenser.reflux_ratio.fix(0.5*pyunits.dimensionless)
    n.fs.D101.reboiler.boilup_ratio.fix(0.5*pyunits.dimensionless)
    n.fs.D101.condenser.condenser_pressure.fix(150000*pyunits.Pa)

    # set scaling factors
    # Set scaling factors for heat duty
    iscale.set_scaling_factor(n.fs.D101.condenser.control_volume.heat, 1e-2)
    iscale.set_scaling_factor(n.fs.D101.reboiler.control_volume.heat, 1e-2)

    # Set the scaling factors for the remaining variables and all constraints
    iscale.calculate_scaling_factors(n.fs.D101)

    n.fs.D101.initialize(outlvl=outlvl)
    results = solver.solve(n, tee=tee)
    from pyomo.environ import TerminationCondition
    assert results.solver.termination_condition == TerminationCondition.optimal
    assert_units_consistent(n)

    print('Complete.')
    print()

    return n
