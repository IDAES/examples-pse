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

__author__ = "Alex Noring"

import pyomo.environ as pyo
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.models_extra.power_generation.costing.power_plant_costing import (
    get_PP_costing,
    get_total_TPC,
    costing_initialization,
    get_fixed_OM_costs,
    get_variable_OM_costs,
    initialize_fixed_OM_costs,
    initialize_variable_OM_costs,
)


def add_total_plant_cost(b, project_contingency, process_contingency):

    # b.costing.total_plant_cost = pyo.Var(bounds=(0, 1e4))
    try:
        b.parent_block().generic_costing_units.append(b)
    except AttributeError:
        b.parent_block().generic_costing_units = []
        b.parent_block().generic_costing_units.append(b)

    @b.costing.Expression()
    def total_plant_cost(c):
        return c.purchase_cost * project_contingency * process_contingency / 1e6


def get_solo_soec_capital_costing(m):
    m.get_costing(year="2018")
    m.costing.total_plant_cost = pyo.Var(
        initialize=130,
        # bounds=(0, 1e4),
        doc="total plant cost in $MM",
    )

    m.soec_module.costing = pyo.Block()
    m.soec_module.costing.total_plant_cost = pyo.Var(
        initialize=130,
        # bounds=(0, 1e4),
        doc="total plant cost in $MM",
    )
    @m.soec_module.costing.Constraint()
    def soec_cost(c):
        return c.total_plant_cost * 1e6 == 60.87 * m.soec_module.number_cells

    # water heaters - U-tube HXs
    # costed with IDAES generic heat exchanger correlation
    m.water_heater01.get_costing(hx_type="U-tube")
    add_total_plant_cost(m.water_heater01, 1.15, 1.15)

    m.water_heater02.get_costing(hx_type="U-tube")
    add_total_plant_cost(m.water_heater02, 1.15, 1.15)

    # sweep hx and the two feed hx's
    # costed with price/ft^2 from NGFC Pathways study
    m.sweep_hx.costing = pyo.Block()
    m.sweep_hx.costing.total_plant_cost = pyo.Var(
        initialize=20,
        bounds=(0, 1e4),
        doc='total plant cost in $MM')
    @m.sweep_hx.costing.Constraint()
    def sweep_hx_cost(c):
        return c.total_plant_cost*1e6 == (
            81.88*pyo.units.convert(m.sweep_hx.area, pyo.units.ft**2))

    m.feed_hx01.costing = pyo.Block()
    m.feed_hx01.costing.total_plant_cost = pyo.Var(
        initialize=20,
        bounds=(0, 1e4),
        doc='total plant cost in $MM')
    @m.feed_hx01.costing.Constraint()
    def feed_hx01_cost(c):
        return c.total_plant_cost*1e6 == (
            81.88*pyo.units.convert(m.feed_hx01.area, pyo.units.ft**2))

    # sweep compressor
    m.sweep_compressor.costing = pyo.Block()
    m.sweep_compressor.costing.total_plant_cost = pyo.Var(
        initialize=2,
        bounds=(0, 1e4),
        doc='total plant cost in $MM')

    # air is coming in at standard conditions
    sweep_compressor_scfm = pyo.units.convert(
        m.sweep_compressor.control_volume.properties_in[0].flow_vol,
        pyo.units.ft**3/pyo.units.min)

    n_sections = 5

    @m.sweep_compressor.costing.Constraint()
    def sweep_compressor_cost(c):  # in 2018 $
        return c.total_plant_cost == n_sections * 10.9995 * pyo.exp(
            0.4692 +
            0.1203*pyo.log(sweep_compressor_scfm/1000/n_sections) +
            0.0931*pyo.log(sweep_compressor_scfm/1000/n_sections)**2) * 1e-3

    # sweep turbine
    # costed with IDAES generic turbine correlation
    m.sweep_turbine.get_costing()
    add_total_plant_cost(m.sweep_turbine, 1.15, 1.15)

    # H2 compressor
    m.cmp01.costing = pyo.Block()
    m.cmp01.costing.total_plant_cost = pyo.Var(
        initialize=1,
        bounds=(0, 1e4),
        doc='total plant cost in $MM')

    h2_comp_process_param = pyo.units.convert(
        m.cmp01.control_volume.properties_out[0].flow_mass,
        pyo.units.lb/pyo.units.hr)

    @m.cmp01.costing.Constraint()
    def h2_comp_cost(c):
        ref_cost = 11.408  # MM$ 2018
        ref_param = 44369*pyo.units.lb/pyo.units.hr
        alpha = 0.7
        return c.total_plant_cost == (
            ref_cost*(h2_comp_process_param/ref_param)**alpha)

    m.cmp02.costing = pyo.Block()
    m.cmp02.costing.total_plant_cost = pyo.Var(
        initialize=1,
        bounds=(0, 1e4),
        doc='total plant cost in $MM')

    h2_comp_process_param = pyo.units.convert(
        m.cmp02.control_volume.properties_out[0].flow_mass,
        pyo.units.lb/pyo.units.hr)

    @m.cmp02.costing.Constraint()
    def h2_comp_cost(c):
        ref_cost = 11.408  # MM$ 2018
        ref_param = 44369*pyo.units.lb/pyo.units.hr
        alpha = 0.7
        return c.total_plant_cost == (
            ref_cost*(h2_comp_process_param/ref_param)**alpha)

    m.cmp03.costing = pyo.Block()
    m.cmp03.costing.total_plant_cost = pyo.Var(
        initialize=1,
        bounds=(0, 1e4),
        doc='total plant cost in $MM')

    h2_comp_process_param = pyo.units.convert(
        m.cmp03.control_volume.properties_out[0].flow_mass,
        pyo.units.lb/pyo.units.hr)

    @m.cmp03.costing.Constraint()
    def h2_comp_cost(c):
        ref_cost = 11.408  # MM$ 2018
        ref_param = 44369*pyo.units.lb/pyo.units.hr
        alpha = 0.7
        return c.total_plant_cost == (
            ref_cost*(h2_comp_process_param/ref_param)**alpha)

    m.cmp04.costing = pyo.Block()
    m.cmp04.costing.total_plant_cost = pyo.Var(
        initialize=1,
        bounds=(0, 1e4),
        doc='total plant cost in $MM')

    h2_comp_process_param = pyo.units.convert(
        m.cmp04.control_volume.properties_out[0].flow_mass,
        pyo.units.lb/pyo.units.hr)

    @m.cmp04.costing.Constraint()
    def h2_comp_cost(c):
        ref_cost = 11.408  # MM$ 2018
        ref_param = 44369*pyo.units.lb/pyo.units.hr
        alpha = 0.7
        return c.total_plant_cost == (
            ref_cost*(h2_comp_process_param/ref_param)**alpha)

    # build constraint summing total plant costs
    get_total_TPC(m)

    # costing initialization
    variables = [
        m.soec_module.costing.total_plant_cost,
        m.cmp01.costing.total_plant_cost,
        m.cmp02.costing.total_plant_cost,
        m.cmp03.costing.total_plant_cost,
        m.cmp04.costing.total_plant_cost,
        m.sweep_hx.costing.total_plant_cost,
        m.feed_hx01.costing.total_plant_cost,
        m.sweep_compressor.costing.total_plant_cost,
    ]

    constraints = [
        m.soec_module.costing.soec_cost,
        m.cmp01.costing.h2_comp_cost,
        m.cmp02.costing.h2_comp_cost,
        m.cmp03.costing.h2_comp_cost,
        m.cmp04.costing.h2_comp_cost,
        m.sweep_hx.costing.sweep_hx_cost,
        m.feed_hx01.costing.feed_hx01_cost,
        m.sweep_compressor.costing.sweep_compressor_cost,
    ]

    for v, c in zip(variables, constraints):
        print(v.name)
        for i in v.keys():
            calculate_variable_from_constraint(v[i], c[i])


    costing_initialization(m)
    calculate_variable_from_constraint(
        m.costing.total_TPC, m.costing.total_TPC_eq
    )


def lock_capital_cost(m):
    for b in m.generic_costing_units:
        for v in b.costing.total_plant_cost.values():
            print(b)
            v.display()
            v.set_value(pyo.value(v))

        b.costing.deactivate()
    m.costing.deactivate()
    m.costing.total_TPC.fix()


def get_soec_OM_costing(m, design_h2_production=2.5 * pyo.units.kg / pyo.units.s):
    # fixed O&M costs
    get_fixed_OM_costs(m.fs, design_h2_production, tech=6)

    @m.Constraint()
    def stack_replacement_cost(fs):
        return fs.costing.other_fixed_costs*1e6 == 4.425*m.soec_module.number_cells

    m.costing.other_fixed_costs.unfix()

    # initialize fixed costs
    calculate_variable_from_constraint(
        m.costing.other_fixed_costs, m.stack_replacement_cost
    )
    initialize_fixed_OM_costs(m.fs)

    # variable O&M costs

    # water use
    @m.Expression(m.time)
    def water_use(fs, t):  # gps
        density = 8.34*pyo.units.lb/pyo.units.gal
        return m.raw_water_withdrawal[t]/density

    # water treatment
    @m.Expression(m.time)
    def water_treatment_use(fs, t):
        density = 8.34*pyo.units.lb/pyo.units.gal
        use_rate = 0.00297886*pyo.units.lb/pyo.units.gal
        return use_rate * m.process_water_discharge[t]/density

    resources = [
        "electricity",
        "water",
        "water treatment chemicals"]
    rates = [
        m.total_electric_power,
        m.water_use,
        m.water_treatment_use]

    prices = {"electricity": 30*pyo.units.USD/pyo.units.MWh}

    get_variable_OM_costs(m.fs, m.h2_mass_production, resources, rates, prices=prices)

    # initialize variable costs
    initialize_variable_OM_costs(m.fs)


def display_soec_costing(m):
    print("Capital cost: ${:.0f}M".format(pyo.value(m.costing.total_TPC)))
    print(
        "Fixed O&M cost: ${:.1f}M/yr".format(
            pyo.value(m.costing.total_fixed_OM_cost)
        )
    )
    print(
        "Electricity cost: ${:.2f}/kg H2".format(
            pyo.value(m.H2_costing.variable_operating_costs[0, "electricity"])
        )
    )
    print(
        "Fuel cost: ${:.2f}/kg H2".format(
            pyo.value(m.H2_costing.variable_operating_costs[0, "natural gas"])
        )
    )
    print(
        "Total variable O&M cost: ${:.2f}/kg H2".format(
            pyo.value(m.H2_costing.total_variable_OM_cost[0])
        )
    )
