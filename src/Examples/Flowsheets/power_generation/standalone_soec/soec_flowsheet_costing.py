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
from pyomo.environ import units as pyunits
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.common.config import ConfigValue
from idaes.core import UnitModelCostingBlock
from idaes.core.util.unit_costing import initialize as initialize_costing_block, pressure_changer_costing
from idaes.core.util.constants import Constants
import idaes.core.util.scaling as iscale
from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption
import idaes.core.base.costing_base as cost_base

from idaes.models_extra.power_generation.costing.power_plant_costing import (
    get_PP_costing,
    costing_initialization,
    get_fixed_OM_costs,
    get_variable_OM_costs,
    initialize_fixed_OM_costs,
    initialize_variable_OM_costs,
)

from idaes.models_extra.power_generation.costing.power_plant_capcost import (
    QGESSCosting,
    QGESSCostingData,
)
from idaes.models.costing.SSLW import(
    SSLWCosting,
    SSLWCostingData,
    HXType)

ssf = iscale.set_scaling_factor
cst = iscale.constraint_scaling_transform


def add_total_plant_cost(b, installation_labor=1.25, eng_fee=1.2,
                         project_contingency=1.15, process_contingency=1.15):

    # b.costing.total_plant_cost = pyo.Var(bounds=(0, 1e4))
    try:
        b.parent_block().generic_costing_units.append(b)
    except AttributeError:
        b.parent_block().generic_costing_units = []
        b.parent_block().generic_costing_units.append(b)

    @b.costing.Expression()
    def total_plant_cost(c):
        return (c.purchase_cost * installation_labor * eng_fee *
                (project_contingency + process_contingency) / 1e6)


def get_total_TPC(b):
    # Have to modify this from power_plant_costing because manual costing blocks renamed ersatz_costing
    # This method accepts either a concrete model or a flowsheet block
    TPC_list = []
    for o in b.component_objects(descend_into=True):
        # look for costing blocks
        costing_block = None
        total_plant_cost = None
        if hasattr(o, "costing"):
            costing_block = o.costing
        elif hasattr(o, "ersatz_costing"):
            costing_block = o.ersatz_costing

        if costing_block is not None and hasattr(costing_block, "total_plant_cost"):
            total_plant_cost = costing_block.total_plant_cost

        if total_plant_cost is not None:
            for key in total_plant_cost.keys():
                TPC_list.append(total_plant_cost[key])

    b.costing.total_TPC = pyo.Var(initialize=0, bounds=(0, 1e4), doc="total TPC in $MM")

    @b.costing.Constraint()
    def total_TPC_eqn(c):
        return c.total_TPC == sum(TPC_list)

def get_solo_soec_capital_costing(fs):
    fs.get_costing(year="2018")
    # fs.costing.total_plant_cost = pyo.Var(
    #     initialize=130,
    #     # bounds=(0, 1e4),
    #     doc="total plant cost in $MM",
    # )

    fs.soec_module.ersatz_costing = pyo.Block()
    fs.soec_module.ersatz_costing.total_plant_cost = pyo.Var(
        initialize=130,
        # bounds=(0, 1e4),
        doc="total plant cost in $MM",
    )
    @fs.soec_module.ersatz_costing.Constraint()
    def total_plant_cost_eqn(b):
        extra_installed_area = 1.10  # include 10% extra stack area to account for degradation over SOEC lifetime
        return b.total_plant_cost * 1e6 == 60.87 * fs.soec_module.number_cells * extra_installed_area
    ssf(fs.soec_module.number_cells, 1e-5)
    cst(fs.soec_module.ersatz_costing.total_plant_cost_eqn, 1e-5)



    # water heaters - U-tube HXs
    # costed with IDAES generic heat exchanger correlation
    water_heaters = [
        fs.water_evaporator01,
        fs.water_evaporator02,
        fs.water_evaporator03,
        fs.water_evaporator04,
        fs.water_evaporator05,
        fs.water_preheater,
    ]
    for hx in water_heaters:
        hx.get_costing(hx_type="U-tube")
        add_total_plant_cost(hx)

    # sweep hx and two feed hx's
    # costed with price/ft^2 from NGFC Pathways study

    cross_flow_exchangers = [
        fs.feed_hot_exchanger,
        fs.sweep_hot_exchanger,
        fs.sweep_medium_exchanger,
    ]
    for hx in cross_flow_exchangers:
        hx.ersatz_costing = pyo.Block()
        hx.ersatz_costing.total_plant_cost = pyo.Var(
            initialize=20,
            bounds=(0, 1e4),
            doc='total plant cost in $MM')
        @hx.ersatz_costing.Constraint()
        def total_plant_cost_eqn(b):
            return b.total_plant_cost*1e6 == (
                81.88*pyo.units.convert(hx.area, pyo.units.ft**2))

        ssf(hx.area, 1e-3)
        cst(hx.ersatz_costing.total_plant_cost_eqn, 1e-5)

    for heater in [fs.feed_heater, fs.sweep_heater]:
        heater.ersatz_costing = pyo.Block()
        heater.ersatz_costing.total_plant_cost = pyo.Var(
            initialize=20,
            bounds=(0, 1e4),
            doc='total plant cost in $MM')
        @heater.ersatz_costing.Constraint()
        def total_plant_cost_eqn(b):
            U = 100 * pyo.units.W / pyo.units.m ** 2 / pyo.units.K
            DT = 50 * pyo.units.K
            area = heater.heat_duty[0] / U / DT
            # Add factor of two to pathways cost to account for corrosion-resistant materials for trim heaters
            return b.total_plant_cost*1e6 == (
                2*81.88*pyo.units.convert(area, pyo.units.ft**2))

        #ssf(hx.area, 1e-3)
        cst(heater.ersatz_costing.total_plant_cost_eqn, 1e-5)


    fs.sweep_blower.get_costing(mover_type="blower", Mat_factor="carbon_steel")
    add_total_plant_cost(fs.sweep_blower)

    flash_vessels = [fs.product_flash01,
                     #fs.product_flash02,
                     fs.product_flash03, fs.product_flash04]
    for flash in flash_vessels:
        flash.get_costing(Mat_factor="low_alloy_steel")
        # Heuristics for flash vessel taken from Analysis, Synthesis, and Design of Chemical Processes by Turton et al.
        # Fourth edition, copyright 2012. Page 344
        @flash.Constraint()
        def length_diameter_heuristic(b):
            return b.length == 3 * b.diameter
        @flash.Constraint()
        def capacity_heuristic(b):
            # Vessel should take 5-10 minutes to fill halfway based on nominal flow rate.
            # Molar volume of water = 1.8068e-5 m^3/mol
            return (Constants.pi / 4 * b.diameter ** 2 * b.length / 2
                    == (1.8068e-5 * pyo.units.m ** 3 / pyo.units.mol) * b.liq_outlet.flow_mol[0] * (10 * 60 *pyo.units.s))
        add_total_plant_cost(flash)
        flash.diameter.bounds = (0, None)
        flash.length.bounds = (0, None)

    # H2 compressors
    reciprocating_piston_compressors = [fs.cmp01, fs.cmp02, fs.cmp03]
    for cmp in reciprocating_piston_compressors:
        cmp.ersatz_costing = pyo.Block()
        cmp.ersatz_costing.total_plant_cost = pyo.Var(
            initialize=1,
            bounds=(0, 1e4),
            doc='total plant cost in $MM')

        @cmp.ersatz_costing.Constraint()
        def total_plant_cost_eqn(b):
            ref_cost = 0.5*11.408  # MM$ 2018 Original number for two stages
            ref_param = 44369*pyo.units.lb/pyo.units.hr
            alpha = 0.7
            h2_comp_process_param = pyo.units.convert(
                cmp.control_volume.properties_out[0].flow_mass,
                pyo.units.lb / pyo.units.hr)
            return b.total_plant_cost == (
                ref_cost*(h2_comp_process_param/ref_param)**alpha)

    # Have to call method because Helmholtz model doesn't automatically call it
    fs.water_compressor.costing = pyo.Block()
    # Despite me specifying the mover type, it's still going to check the ThermodynamicAssumption
    fs.water_compressor.config.declare(
        "thermodynamic_assumption",
        ConfigValue(
            default=ThermodynamicAssumption.isentropic
        ),
    )
    fs.water_compressor.config.declare(
        "compressor",
        ConfigValue(
            default=True
        ),
    )
    pressure_changer_costing(fs.water_compressor.costing, mover_type="compressor", compressor_type="screw")
    add_total_plant_cost(fs.water_compressor)

    fs.cmp04.ersatz_costing = pyo.Block()
    fs.cmp04.ersatz_costing.total_plant_cost = pyo.Var(
            initialize=1,
            bounds=(0, 1e4),
            doc='total plant cost in $MM'
    )
    fs.cmp04.ersatz_costing.number_units = pyo.Param(
        initialize=2,
        domain=pyo.Integers,
        mutable=True
    )
    @fs.cmp04.ersatz_costing.Constraint()
    def total_plant_cost_eqn(b):
        ref_cost = (
            b.number_units
            * 583.1/556.8  # Correct from ChemE index from 2015 to 2018
            * 4.05  # Cost in $MM, increased to take into account sixth intercooler
        )
        ref_param = b.number_units * 240e3 * pyo.units.kg/pyo.units.day  # ~2.7 kg/s
        # We cannot use this exponent to scale up because the original design is constrained by material limitations,
        # but we should be able to scale down using standard centrifugal compressor exponent
        alpha = 0.8
        process_param = pyo.units.convert(fs.cmp04.inlet.flow_mol[0] * 0.002016 * (pyo.units.kg / pyo.units.mol),
                                          pyo.units.kg/pyo.units.day)
        return b.total_plant_cost == ref_cost * (process_param/ref_param) ** alpha

    fs.heat_pump.ersatz_costing = pyo.Block()
    fs.heat_pump.ersatz_costing.total_plant_cost = pyo.Var(
            initialize=1,
            bounds=(0, 1e4),
            doc='total plant cost in $MM'
    )

    @fs.heat_pump.ersatz_costing.Constraint()
    def total_plant_cost_eqn(b):
        equipment_cost = (583.1/550 # Correct for the ChemE index from 2008-2009 to that of 2018
                             * 1.25 # $MM
                             * 2) # Safety factor
        ref_cost = equipment_cost * 1.25 * 1.2 * 1.3 # 25% installation labor, 20% engineering fee, 15% + 15% process and product contingencies
        ref_param = 10275e3 * pyo.units.BTU/pyo.units.hr
        process_param = pyo.units.convert(-fs.heat_pump.heat_out[0], pyo.units.BTU/pyo.units.hr)
        alpha = 1
        return b.total_plant_cost == ref_cost*(process_param/ref_param)**alpha

    ##########################################################################
    # Power Plant Costing Library

    # accounts with feedwater flow to the SOEC as the process parameter
    fw_accounts = ["3.1", "3.3"]

    fs.fw_system = pyo.Block()

    @fs.Expression(fs.time)
    def feedwater_flow(b, t):
        return pyo.units.convert(
            fs.water_preheater.tube.properties_in[t].flow_mass,
            pyo.units.lb/pyo.units.hr)

    get_PP_costing(
        fs.fw_system, fw_accounts, fs.feedwater_flow[0], "lb/hr", 6)

    # accounts with raw water withdrawal as the process parameter
    water_withdrawl_accounts = ["3.2", "3.4", "9.5", "14.6"]

    fs.water_withdrawl_system = pyo.Block()

    # need to do some water accounting here
    @fs.Expression(fs.time)
    def water_demand(b, t):
        return (
            fs.water_preheater.tube.properties_in[t].flow_mol
            + fs.heat_source.control_volume.properties_in[t].flow_mol
        )

    @fs.Expression(fs.time)
    def water_recycle(b, t):
        return (
            fs.product_flash01.liq_outlet.flow_mol[t]
            + fs.product_flash03.liq_outlet.flow_mol[t]
            + fs.product_flash04.liq_outlet.flow_mol[t]
            + fs.product_flash05.liq_outlet.flow_mol[t]
        )

    @fs.Expression(fs.time)
    def raw_water_withdrawal(b, t):
        # Molar volume of water = 1.8068e-5 m^3/mol
        #return pyo.units.convert(
        # 1.8068e-5 * pyo.units.m**3/pyo.units.mol
        # * (fs.water_demand[t] - fs.water_recycle[t]),
        # pyo.units.gal/pyo.units.min
        return 3000 * pyo.units.gal/pyo.units.min

    get_PP_costing(
        fs.water_withdrawl_system, water_withdrawl_accounts,
        fs.raw_water_withdrawal[0], "gpm", 6)

    # accounts with process water discharge as the process parameter
    water_discharge_accounts = ["3.7"]

    fs.discharge_system = pyo.Block()

    @fs.Expression(fs.time)
    def process_water_discharge(b, t):
        # return pyo.units.convert(
        #     fs.heat_source.control_volume.properties_out[0].flow_vol,
        #     pyo.units.gal/pyo.units.min)
        return 3000 * pyo.units.gal/pyo.units.min

    get_PP_costing(
        fs.discharge_system, water_discharge_accounts,
        fs.process_water_discharge[0], "gpm", 6)

    # accounts with net power as the process parameter
    # includes some buildings, site improvements, and electrical equipment
    net_power_accounts = ["11.1", "11.7", "11.9",
                          "13.1", "13.2", "13.3",
                          "14.4", "14.7", "14.8", "14.9", "14.10"]

    fs.net_power_costs = pyo.Block()

    get_PP_costing(
        fs.net_power_costs, net_power_accounts,
        pyo.units.convert(fs.total_electric_power[0], pyo.units.kW), "kW", 6)

    # costing for transformers
    transformer_accounts = ["11.8"]

    fs.transformers = pyo.Block()

    get_PP_costing(
        fs.transformers, transformer_accounts,
        pyo.units.convert(fs.total_electric_power[0], pyo.units.kW), "kW", 6)

    # accounts with auxilliary load (excluding the SOEC) as the process parameter
    aux_load_accounts = ["11.2", "11.3", "11.4", "11.5", "11.6",
                         "12.4", "12.5", "12.6",
                         "12.7", "12.8", "12.9"]

    fs.aux_load_costs = pyo.Block()

    @fs.Expression(fs.time)
    def aux_load(b, t):
        return fs.total_electric_power[t] - fs.soec_module.electrical_work[t]

    get_PP_costing(fs.aux_load_costs, aux_load_accounts,
                   pyo.units.convert(fs.aux_load[0], pyo.units.kW), "kW", 6)

    # adding expressions to group Fossil Energy Baseline costs
    @fs.costing.Expression()
    def water_systems_cost(b):
        return (fs.fw_system.costing.total_plant_cost["3.1"] +
                fs.water_withdrawl_system.costing.total_plant_cost["3.2"] +
                fs.fw_system.costing.total_plant_cost["3.3"] +
                fs.water_withdrawl_system.costing.total_plant_cost["3.4"] +
                fs.discharge_system.costing.total_plant_cost["3.7"] +
                fs.water_withdrawl_system.costing.total_plant_cost["9.5"])

    @fs.costing.Expression()
    def accessory_electric_plant_cost(b):
        return (fs.net_power_costs.costing.total_plant_cost["11.1"] +
                fs.aux_load_costs.costing.total_plant_cost["11.2"] +
                fs.aux_load_costs.costing.total_plant_cost["11.3"] +
                fs.aux_load_costs.costing.total_plant_cost["11.4"] +
                fs.aux_load_costs.costing.total_plant_cost["11.5"] +
                fs.aux_load_costs.costing.total_plant_cost["11.6"] +
                fs.net_power_costs.costing.total_plant_cost["11.7"] +
                fs.transformers.costing.total_plant_cost["11.8"] +
                fs.net_power_costs.costing.total_plant_cost["11.9"])

    @fs.costing.Expression()
    def instrumentation_and_control_cost(b):
        return (fs.aux_load_costs.costing.total_plant_cost["12.4"] +
                fs.aux_load_costs.costing.total_plant_cost["12.5"] +
                fs.aux_load_costs.costing.total_plant_cost["12.6"] +
                fs.aux_load_costs.costing.total_plant_cost["12.7"] +
                fs.aux_load_costs.costing.total_plant_cost["12.8"] +
                fs.aux_load_costs.costing.total_plant_cost["12.9"])

    @fs.costing.Expression()
    def improvements_to_site_cost(b):
        return (fs.net_power_costs.costing.total_plant_cost["13.1"] +
                fs.net_power_costs.costing.total_plant_cost["13.2"] +
                fs.net_power_costs.costing.total_plant_cost["13.3"])

    @fs.costing.Expression()
    def buildings_and_structures_cost(b):
        return (fs.net_power_costs.costing.total_plant_cost["14.4"] +
                fs.water_withdrawl_system.costing.total_plant_cost["14.6"] +
                fs.net_power_costs.costing.total_plant_cost["14.7"] +
                fs.net_power_costs.costing.total_plant_cost["14.8"] +
                fs.net_power_costs.costing.total_plant_cost["14.9"] +
                fs.net_power_costs.costing.total_plant_cost["14.10"])

    # build constraint summing total plant costs
    # Get total total plant cost.
    get_total_TPC(fs)

    @fs.costing.Expression()
    def owners_cost(b):
        return 0.21 * b.total_TPC

    @fs.costing.Expression()
    def total_overnight_cost(b):
        return b.total_TPC + b.owners_cost

    @fs.costing.Expression()
    def total_as_spent_cost(b):
        return 1.093 * b.total_overnight_cost

    @fs.costing.Expression()
    def total_annualized_cost(b):
        return 0.0707 * b.total_as_spent_cost

@declare_process_block_class("SOFCPathwaysCosting")
class SOFCPathwaysCostingData(FlowsheetCostingBlockData):

    # Register currency and conversion rates based on CE Index
    register_idaes_currency_units()

    def build_global_params(self):
        """
        This is where we can declare any global parameters we need, such as
        Lang factors, or coefficients for costing methods that should be
        shared across the process.

        You can do what you want here, so you could have e.g. sub-Blocks
        for each costing method to separate the parameters for each method.
        """
        # Set the base year for all costs
        self.base_currency = pyo.units.USD_2018
        # Set a base period for all operating costs
        self.base_period = pyo.units.year

    def build_process_costs(self):
        """
        This is where you do all your process wide costing.
        This is completely up to you, but you will have access to the
        following aggregate costs:

        1. self.aggregate_capital_cost
        2. self.aggregate_fixed_operating_cost
        3. self.aggregate_variable_operating_cost
        4. self.aggregate_flow_costs (indexed by flow type)
        """
        # TODO: Do we have any process level methods to add here?
        pass

    @staticmethod
    def initialize_build(self):
        """
        Here we can add initialization steps for the things we built in
        build_process_costs.

        Note that the aggregate costs will be initialized by the framework.
        """
        # TODO: For now,  no additional process level costs to initialize
        pass

    def cost_solid_oxide_cell(blk):
        """
        SOFC\SOEC costing method based on number of cells, assuming 400 mA/cm^2 per cell in SOFC mode
        and 1000 mA/cm^2 in SOEC mode.
        """

        #TODO Is replacement cost necessary here?
        # @fs.soec_module.ersatz_costing.Expression()
        # def annual_soec_replacement_cost(b):
        #     # In MM$
        #     return 4.2765e-6 * fs.soec_module.number_cells
        @blk.Expression()
        def total_plant_cost(c):
            extra_installed_area = 1.10  # accounts for cell degradation
            return pyunits.convert(60.87 * blk.unit_model.number_cells *
                                   extra_installed_area * pyunits.USD_2018,
                                   to_units=blk.base_currency)

    def cost_cross_flow_heat_exchanger(blk):
        """
        Cross flow heat exchangers
        """
        @blk.Expression()
        def total_plant_cost(c):
            return pyunits.convert(81.88 * pyo.value(
                pyo.units.convert(blk.unit_model.area, pyo.units.ft ** 2)) *
                                   pyunits.USD_2018, to_units=blk.base_currency)

    def cost_trim_heater(blk):
        """
        Trim heaters
        """
        @blk.costing.Expression()
        def total_plant_cost(b):
            U = 56 * pyo.units.W / pyo.units.m ** 2 / pyo.units.K
            DT = 20 * pyo.units.K
            area = blk.heat_duty[0] / U / DT
            # Add factor of two to pathways cost to account for corrosion-resistant materials for trim heaters
            return pyunits.convert(
                2*81.88 * pyo.value(
                    pyo.units.convert(area, pyo.units.ft**2)
                    ) *
                pyunits.USD_2018,
                to_units=blk.base_currency
            )
    def cost_centrifugal_hydrogen_compressor(blk):
        """
        Costing method for the six-stage centrifugal compressor
        designed for hydrogen pipeline service.
        Di Bella, Francis A. Development of a centrifugal hydrogen pipeline gas
        compressor. No. TM-1785. Concepts ETI, Inc. dba Concepts NREC, 2015.

        The compressor is rated to compress 240,000 kg/day of H2, ~2.78 kg/s from
        350 psig to 1250 psig. Intercoolers are used to cool hydrogen to 105 F between
        stages. A cooler after the final stage is not included in the original study,
        but the cost has adjusted upwards here to account for the addition of one.
        """
        blk.number_units = pyo.Param(
            initialize=2,
            domain=pyo.Integers,
            mutable=True
        )

        @blk.Expression()
        def total_plant_cost(b):
            ref_cost = (
                b.number_units
                * 583.1/556.8  # Correct from ChemE index from 2015 to 2018
                * 4.05  # Cost in $MM, increased to take into account sixth intercooler
            )
            ref_param = b.number_units * 240e3 * pyo.units.kg/pyo.units.day  # ~2.78 kg/s
            # We cannot use this exponent to scale up because the original design is constrained by material limitations,
            # but we should be able to scale down using standard centrifugal compressor exponent
            alpha = 0.8
            process_param = pyo.units.convert(blk.unit_model.inlet.flow_mol[0] * 0.002016 * (pyo.units.kg / pyo.units.mol),
                                              pyo.units.kg/pyo.units.day)
            return b.total_plant_cost == ref_cost * (process_param/ref_param) ** alpha

    def cost_heat_pump(blk):
        """
        Heat pump costing based on information found in an Emerson white paper
        located at https://climate.emerson.com/documents/vilter-heat-pump-white-paper-en-us-5411194.pdf
        Based off of the example providing 190 F water for dairy HTST, with a safety factor
        of two because it will be a new system rather than one linked into an existing refrigeration
        system. We're not accounting for makeup ammonia for the heat pump, but we're almost certainly
        overestimating how much we need to pay for water so hopefully things will balance out
        """
        @blk.Expression()
        def total_plant_cost(b):
            equipment_cost = (583.1 / 550  # Correct for the ChemE index from 2008-2009 to that of 2018
                              * 1.25  # $MM
                              * 2)  # Safety factor
            ref_cost = equipment_cost * 1.25 * 1.2 * 1.3  # 25% installation labor, 20% engineering fee, 15% + 15% process and product contingencies
            ref_param = 10275e3 * pyo.units.BTU / pyo.units.hr
            process_param = pyo.units.convert(-b.unit_model.heat_out[0], pyo.units.BTU / pyo.units.hr)
            alpha = 1
            return b.total_plant_cost == ref_cost * (process_param / ref_param) ** alpha

def get_soec_capital_cost(m, CE_index_year):
    CE_index_units = getattr(pyunits, "MUSD_" + CE_index_year)  # millions of USD in base year

    ###########################################################################
    # SOEC CAPITAL COSTING

    # soec modules
    m.fs.soec_module.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=SSLWCostingData.cost_none,
    )

    @m.fs.soec.soec_module.costing.Expression()
    def total_plant_cost(c):
        extra_installed_area = 1.10  # accounts for cell degradation
        return pyunits.convert(60.87 * m.fs.soec_module.number_cells *
                               extra_installed_area * pyunits.USD_2018,
                               CE_index_units)

    # high temperature soec trim heaters
    for heater in [m.fs.feed_heater, m.fs.sweep_heater]:
        # create dummy area variables
        heater.area = pyo.Var(initialize=0, units=pyunits.m ** 2)
        heater.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=SSLWCostingData.cost_none,
        )

        @heater.costing.Expression()
        def total_plant_cost(b):
            #TODO revisit numbers
            U = 56 * pyo.units.W / pyo.units.m ** 2 / pyo.units.K
            DT = 20 * pyo.units.K
            area = heater.heat_duty[0] / U / DT
            # Add factor of two to pathways cost to account for corrosion-resistant materials for trim heaters
            return pyunits.convert(
                2 * 81.88 * pyo.value(
                    pyo.units.convert(area, pyo.units.ft ** 2)
                ) *
                pyunits.USD_2018,
                CE_index_units)

        # sweep and feed hx
        # costed with price/ft^2 from NGFC Pathways study
        for hx in [m.fs.sweep_hot_exchanger, m.fs.sweep_medium_exchanger, m.fs.feed_hot_exchanger]:
            hx.costing = UnitModelCostingBlock(
                flowsheet_costing_block=m.fs.costing,
                costing_method=SSLWCostingData.cost_none,
            )

            @hx.costing.Expression()
            def total_plant_cost(c):
                return pyunits.convert(81.88 * pyo.value(
                    pyo.units.convert(hx.area, pyo.units.ft ** 2)) *
                                       pyunits.USD_2018, CE_index_units)



def initialize_flowsheet_costing(fs):
    flash_vessels = [fs.product_flash01,
                     #fs.product_flash02,
                        fs.product_flash03, fs.product_flash04]
    for flash in flash_vessels:
        flash.diameter.value = pyo.value(6*1.81e-3*8*flash.liq_outlet.flow_mol[0]/(3*Constants.pi))**(1/3)
        flash.length.value = 3 * flash.diameter.value
    for o in fs.component_objects(descend_into=True):
        # look for costing blocks
        if hasattr(o, "costing"):
            if hasattr(o.costing, "purchase_cost"):
                initialize_costing_block(o.costing)
            else:
                for key in o.costing.bare_erected_cost.keys():
                    calculate_variable_from_constraint(
                        o.costing.bare_erected_cost[key],
                        o.costing.bare_erected_cost_eq[key],
                    )
                    calculate_variable_from_constraint(
                        o.costing.total_plant_cost[key],
                        o.costing.total_plant_cost_eq[key],
                    )
        elif hasattr(o, "ersatz_costing"):
            calculate_variable_from_constraint(
                o.ersatz_costing.total_plant_cost,
                o.ersatz_costing.total_plant_cost_eqn
            )
    calculate_variable_from_constraint(
        fs.costing.total_TPC,
        fs.costing.total_TPC_eqn
    )

def scale_flowsheet_costing(fs):
    ssf = iscale.set_scaling_factor
    cst = iscale.constraint_scaling_transform

    ssf(fs.sweep_blower.costing.base_cost_per_unit, 1e-5)
    ssf(fs.sweep_blower.costing.purchase_cost, 1e-5)

    ssf(fs.water_compressor.costing.base_cost_per_unit, 1e-6)
    ssf(fs.water_compressor.costing.purchase_cost, 1e-6)
    # Since I'm hacking this costing method into the Helmholtz compressor, I need to scale this equations manually
    cst(fs.water_compressor.costing.base_cost_per_unit_eq, 1e-6)
    cst(fs.water_compressor.costing.cp_cost_eq, 1e-6)

    flash_vessels = [fs.product_flash01,
                     #fs.product_flash02,
                     fs.product_flash03, fs.product_flash04]
    for flash in flash_vessels:
        ssf(flash.costing.base_cost_per_unit, 3e-5)
        ssf(flash.costing.purchase_cost, 3e-5)
        ssf(flash.costing.vessel_purchase_cost, 1e-4)
        ssf(flash.costing.weight, 5e-4)
        ssf(flash.costing.base_cost_platf_ladders, 3e-4)
        cst(flash.costing.cv_cost_eq, 3e-5)
        cst(flash.costing.CPL_eq, 3e-4)
        cst(flash.costing.weight_eq, 5e-4)
        cst(flash.costing.cp_vessel_eq, 1e-4)
        cst(flash.costing.cp_cost_eq, 3e-5)
    water_heaters = [
        fs.water_evaporator01,
        fs.water_evaporator02,
        fs.water_evaporator03,
        fs.water_evaporator04,
        fs.water_evaporator05,
        fs.water_preheater,
    ]
    for hx in water_heaters:
        ssf(hx.costing.base_cost_per_unit, 1e-5)
        ssf(hx.costing.purchase_cost, 1e-5)


def lock_capital_cost(fs):
    for b in fs.generic_costing_units:
        for v in b.costing.total_plant_cost.values():
            print(b)
            v.display()
            v.set_value(pyo.value(v))

        b.costing.deactivate()
    fs.costing.deactivate()
    fs.costing.total_TPC.fix()


def get_soec_OM_costing(fs):
    # fixed O&M costs
    @fs.costing.Expression()
    def annual_operating_labor_cost(b):
        return (
                6.3  # Five operators for 3 shifts + redundancy.
                * 38.50  # Hourly rate in 2018 $/hr
                * 24 * 365.2425  # Average hours in Gregorian year
                * 1.3  # 30% Labor burden
                * 1e-6  # To put in $MM
        )

    @fs.costing.Expression()
    def maintenance_labor_cost(b):
        return b.total_TPC * 0.4 * 0.019

    @fs.costing.Expression()
    def maintenance_material_cost(b):
        return b.total_TPC * 0.6 * 0.019

    @fs.costing.Expression()
    def admin_and_support_labor_cost(b):
        return 0.25 * (b.annual_operating_labor_cost + b.maintenance_labor_cost)

    @fs.costing.Expression()
    def property_tax_and_insurance_cost(b):
        return 0.02 * b.total_TPC

    @fs.costing.Expression()
    def annual_fixed_operations_and_maintenance_cost(b):
        return (
                b.annual_operating_labor_cost
                + b.maintenance_labor_cost
                + b.maintenance_material_cost
                + b.admin_and_support_labor_cost
                + b.property_tax_and_insurance_cost
                + fs.soec_module.ersatz_costing.annual_soec_replacement_cost
        )

    # variable O&M costs

    cost_base.register_idaes_currency_units()
    fs.costing.water_price = pyo.Var(
        initialize=9.352e-4,
        units=pyo.units.USD_2018 / pyo.units.kg,
        doc="Price of water including treatment chemicals",
    )
    fs.costing.water_price.fix()
    fs.costing.electricity_price = pyo.Var(
        initialize=71.70,
        units=pyo.units.USD_2018 / pyo.units.MWh,
        doc="Price of electricity",
    )
    fs.costing.electricity_price.fix()
    fs.costing.air_price = pyo.Var(
        initialize=(
                0.35  # Price of dry air at 3.3 barg in dollars per 100 standard (1 ATM and 15 C) cubic meters
                / 100  # 100 standard meters cubed to 1 standard meter cubed
                * 0.02364  # Molar volume of dry air at standard conditions in m^3/mol
                * 671.1 / 584.6  # Correct from 2012 to 2018 using cost index numbers from unit_costing.py
                * 2/3  # Pressurizing air to 3.3 at 71.70 $/MWh isothermally makes up 72.5% of this cost
                       # We probably can't get away with that much of a reduction, but multiply price
                       # by 2/3rds to attempt to account for different needs and on-site generation
                ),
        units=pyo.units.USD_2018/pyo.units.mol
    )
    fs.costing.air_price.fix()

    fs.costing.plant_uptime = pyo.Var(
        initialize=3.1556952e7,
        units=pyo.units.s,
        doc="Amount of time plant is run in year",
    )
    fs.costing.plant_uptime.fix()

    @fs.costing.Expression()
    def annual_electricity_cost(b):
        return (
                fs.total_electric_power[0]
                * b.plant_uptime
                * pyo.units.convert(b.electricity_price, pyo.units.USD_2018 / pyo.units.J)
                * 1e-6 / pyo.units.USD_2018  # To put in MM$
        )

    @fs.costing.Expression()
    def annual_water_cost(b):
        return (
                b.water_price
                * b.plant_uptime
                * 0.01802 * pyo.units.kg / pyo.units.mol
                * (fs.water_demand[0] - fs.water_recycle[0])
                * 1e-6 / pyo.units.USD_2018  # To put in MM$
        )

    @fs.costing.Expression()
    def annual_air_cost(b):
        return (
            b.air_price
            * b.plant_uptime
            * fs.sweep_blower.inlet.flow_mol[0]
            * 1e-6 / pyo.units.USD_2018  # To put in MM$
        )


def display_soec_costing(fs):
    print("Capital cost: ${:.0f}M".format(pyo.value(fs.costing.total_TPC)))
    print(
        "Fixed O&M cost: ${:.1f}M/yr".format(
            pyo.value(fs.costing.annual_fixed_operations_and_maintenance_cost)
        )
    )
    print(
        "Electricity cost: ${:.2f}/kg H2".format(
            pyo.value(1e6 * fs.costing.annual_electricity_cost / (fs.h2_mass_production[0] * fs.costing.plant_uptime))
        )
    )
    print(
        "Water cost: ${:.2f}/kg H2".format(
            pyo.value(1e6 * fs.costing.annual_water_cost / (fs.h2_mass_production[0] * fs.costing.plant_uptime))
        )
    )
    print(
        "Air cost: ${:.2f}/kg H2".format(
            pyo.value(1e6 * fs.costing.annual_air_cost / (fs.h2_mass_production[0] * fs.costing.plant_uptime))
        )
    )