
Overview
--------

This notebook is based on a workshop tutorial originally given at FOCAPD
in 2019. It demonstrated the following capabilities:

-  Construct a steady-state flowsheet using the IDAES unit model library
-  Connecting unit models in a flowsheet using Arcs
-  Using the SequentialDecomposition tool to initialize a flowsheet with
   recycle
-  Fomulate and solve an optimization problem

   -  Defining an objective function
   -  Setting variable bounds
   -  Adding additional constraints

Changes were made to use the notebook to demonstrate use of the IDAES
Data Management Framework (DMF): - The underlying Python modules were
adapted to use the DMF for loading property data - The optimization
problem was removed to shorten the example

Problem Statement
-----------------

Hydrodealkylation is a chemical reaction that often involves reacting an
aromatic hydrocarbon in the presence of hydrogen gas to form a simpler
aromatic hydrocarbon devoid of functional groups,. In this example,
toluene will be reacted with hydrogen gas at high temperatures to form
benzene via the following reaction:

**C6H5CH3 + H2 → C6H6 + CH4**

This reaction is often accompanied by an equilibrium side reaction which
forms diphenyl, which we will neglect for this example.

This example is based on the 1967 AIChE Student Contest problem as
present by Douglas, J.M., Chemical Design of Chemical Processes, 1988,
McGraw-Hill.

The flowsheet that we will be using for this module is shown below with
the stream conditions. We will be processing toluene and hydrogen to
produce at least 370 TPY of benzene. As shown in the flowsheet, there
are two flash tanks, F101 to separate out the non-condensibles and F102
to further separate the benzene-toluene mixture to improve the benzene
purity. Note that typically a distillation column is required to obtain
high purity benzene but that is beyond the scope of this workshop. The
non-condensibles separated out in F101 will be partially recycled back
to M101 and the rest will be either purged or combusted for power
generation.We will assume ideal gas for this flowsheet. The properties
required for this module are available in the same directory:

-  hda_ideal_VLE.py
-  hda_reaction.py

The state variables chosen for the property package are **flows of
component by phase, temperature and pressure**. The components
considered are: **toluene, hydrogen, benzene and methane**. Therefore,
every stream has 8 flow variables, 1 temperature and 1 pressure
variable.

Data Management Framework
~~~~~~~~~~~~~~~~~~~~~~~~~

In this version of the example, the DMF “workspace” is used to store
these state variables. This involves the use of additional files:

-  config.yaml - DMF Configuration
-  resourcedb.json - The DMF’s file-based datastore, where the values
   and their provenance are stored.

|image0|

.. |image0| image:: module_2_flowsheet.png

Importing required pyomo and idaes components
---------------------------------------------

To construct a flowsheet, we will need several components from the pyomo
and idaes package. Let us first import the following components from
Pyomo: - Constraint (to write constraints) - Var (to declare variables)
- ConcreteModel (to create the concrete model object) - Expression (to
evaluate values as a function of variables defined in the model) -
Objective (to define an objective function for optimization) -
SolverFactory (to solve the problem) - TransformationFactory (to apply
certain transformations) - Arc (to connect two unit models) -
SequentialDecomposition (to initialize the flowsheet in a sequential
mode)

For further details on these components, please refer to the pyomo
documentation: https://pyomo.readthedocs.io/en/latest/

.. code:: ipython3

    from pyomo.environ import (Constraint,
                               Var,
                               ConcreteModel,
                               Expression,
                               Objective,
                               SolverFactory,
                               TransformationFactory,
                               value)
    from pyomo.network import Arc, SequentialDecomposition

From idaes, we will be needing the FlowsheetBlock and the following unit
models: - Mixer - Heater - StoichiometricReactor - Flash - Separator
(splitter) - PressureChanger

.. code:: ipython3

    from idaes.core import FlowsheetBlock

.. code:: ipython3

    from idaes.unit_models import (PressureChanger,
                                   Mixer,
                                   Separator as Splitter,
                                   Heater,
                                   StoichiometricReactor)

.. raw:: html

   <div class="alert alert-block alert-info">

Inline Exercise: Now, import the remaining unit models highlighted in
blue above and run the cell using ``Shift+Enter`` after typing in the
code.

.. raw:: html

   </div>

.. code:: ipython3

    from idaes.unit_models import Flash

We will also be needing some utility tools to put together the flowsheet
and calculate the degrees of freedom.

.. code:: ipython3

    from idaes.unit_models.pressure_changer import ThermodynamicAssumption
    from idaes.core.util.model_statistics import degrees_of_freedom

Importing required thermo and reaction package
----------------------------------------------

The final set of imports are to import the thermo and reaction package
for the HDA process. We have created a custom thermo package that
assumes Ideal Gas with support for VLE.

The reaction package here is very simple as we will be using only a
StochiometricReactor and the reaction package consists of the
stochiometric coefficients for the reaction and the parameter for the
heat of reaction.

Let us import the following modules and they are in the same directory
as this jupyter notebook:

.. raw:: html

   <ul>

.. raw:: html

   <li>

hda_ideal_VLE as thermo_props

.. raw:: html

   </li>

.. raw:: html

   <li>

hda_reaction as reaction_props

.. raw:: html

   </li>

.. raw:: html

   </ul>

.. raw:: html

   </div>

.. code:: ipython3

    import hda_ideal_VLE as thermo_props
    import hda_reaction as reaction_props

Constructing the Flowsheet
--------------------------

We have now imported all the components, unit models, and property
modules we need to construct a flowsheet. Let us create a ConcreteModel
and add the flowsheet block as we did in module 1.

.. code:: ipython3

    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

We now need to add the property packages to the flowsheet. Unlike Module
1, where we only had a thermo property package, for this flowsheet we
will also need to add a reaction property package.

.. code:: ipython3

    m.fs.thermo_params = thermo_props.HDAParameterBlock()
    m.fs.reaction_params = reaction_props.HDAReactionParameterBlock(
            default={"property_package": m.fs.thermo_params})

Viewing imported properties with the DMF
----------------------------------------

Since the property package above used the DMF to load the constant
values used for calculating property values, we can look at the
provenance and data programmatically.

.. code:: ipython3

    from idaes.dmf.resource import TidyUnitData as TDU
    prop_name = 'dh_vap'
    rsrc = m.fs.thermo_params.get_property_resource(prop_name)
    print(f"Property {prop_name} source: {rsrc.formatted_source()}")
    print(f"Property {prop_name} data:")
    TDU(rsrc.data).table


.. parsed-literal::

    Property dh_vap source: Reid, R. C., Prausnitz, J. M., & Poling, B. E. (1987). The properties of gases and liquids. ISBN: 9780070517998 Date: 1987
    Property dh_vap data:




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>compound</th>
          <th>value</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>benzene</td>
          <td>33870.0</td>
        </tr>
        <tr>
          <th>1</th>
          <td>toluene</td>
          <td>38262.0</td>
        </tr>
        <tr>
          <th>2</th>
          <td>hydrogen</td>
          <td>0.0</td>
        </tr>
        <tr>
          <th>3</th>
          <td>methane</td>
          <td>0.0</td>
        </tr>
        <tr>
          <th>4</th>
          <td>diphenyl</td>
          <td>62710.0</td>
        </tr>
      </tbody>
    </table>
    </div>



Adding Unit Models
------------------

Let us start adding the unit models we have imported to the flowsheet.
Here, we are adding the Mixer (assigned a name M101) and a Heater
(assigned a name H101). Note that, all unit models need to be given a
property package argument. In addition to that, there are several
arguments depending on the unit model, please refer to the documentation
for more details
(https://idaes-pse.readthedocs.io/en/latest/models/index.html). For
example, the Mixer unit model here is given a ``list`` consisting of
names to the three inlets.

.. code:: ipython3

    m.fs.M101 = Mixer(default={"property_package": m.fs.thermo_params,
                               "inlet_list": ["toluene_feed", "hydrogen_feed", "vapor_recycle"]})
    
    m.fs.H101 = Heater(default={"property_package": m.fs.thermo_params,
                                "has_pressure_change": False,
                                "has_phase_equilibrium": True})

.. raw:: html

   <div class="alert alert-block alert-info">

Inline Exercise: Let us now add the StoichiometricReactor(assign the
name R101) and pass the following arguments:

.. raw:: html

   <ul>

.. raw:: html

   <li>

“property_package”: m.fs.thermo_params

.. raw:: html

   </li>

.. raw:: html

   <li>

“reaction_package”: m.fs.reaction_params

.. raw:: html

   </li>

.. raw:: html

   <li>

“has_heat_of_reaction”: True

.. raw:: html

   </li>

.. raw:: html

   <li>

“has_heat_transfer”: True

.. raw:: html

   </li>

.. raw:: html

   <li>

“has_pressure_change”: False

.. raw:: html

   </li>

.. raw:: html

   </ul>

.. raw:: html

   </div>

.. code:: ipython3

    m.fs.R101 = StoichiometricReactor(
                default={"property_package": m.fs.thermo_params,
                         "reaction_package": m.fs.reaction_params,
                         "has_heat_of_reaction": True,
                         "has_heat_transfer": True,
                         "has_pressure_change": False})

Let us now add the Flash(assign the name F101) and pass the following
arguments:

.. raw:: html

   <ul>

.. raw:: html

   <li>

“property_package”: m.fs.thermo_params

.. raw:: html

   </li>

.. raw:: html

   <li>

“has_heat_transfer”: True

.. raw:: html

   </li>

.. raw:: html

   <li>

“has_pressure_change”: False

.. raw:: html

   </li>

.. raw:: html

   </ul>

.. code:: ipython3

    m.fs.F101 = Flash(default={"property_package": m.fs.thermo_params,
                                   "has_heat_transfer": True,
                                   "has_pressure_change": True})

Let us now add the Splitter(S101), PressureChanger(C101) and the second
Flash(F102).

.. code:: ipython3

    m.fs.S101 = Splitter(default={"property_package": m.fs.thermo_params,
                                   "ideal_separation": False,
                                   "outlet_list": ["purge", "recycle"]})
        
    
    m.fs.C101 = PressureChanger(default={
                "property_package": m.fs.thermo_params,
                "compressor": True,
                "thermodynamic_assumption": ThermodynamicAssumption.isothermal})
        
    m.fs.F102 = Flash(default={"property_package": m.fs.thermo_params,
                               "has_heat_transfer": True,
                               "has_pressure_change": True})

Connecting Unit Models using Arcs
---------------------------------

We have now added all the unit models we need to the flowsheet. However,
we have not yet specifed how the units are to be connected. To do this,
we will be using the ``Arc`` which is a pyomo component that takes in
two arguments: ``source`` and ``destination``. Let us connect the outlet
of the mixer(M101) to the inlet of the heater(H101).

.. code:: ipython3

    m.fs.s03 = Arc(source=m.fs.M101.outlet, destination=m.fs.H101.inlet)

|image0|

.. raw:: html

   <div class="alert alert-block alert-info">

Inline Exercise: Now, connect the H101 outlet to the R101 inlet using
the cell above as a guide.

.. raw:: html

   </div>

.. |image0| image:: module_2_flowsheet.png

.. code:: ipython3

    m.fs.s04 = Arc(source=m.fs.H101.outlet, destination=m.fs.R101.inlet)

We will now be connecting the rest of the flowsheet as shown below.
Notice how the outlet names are different for the flash tanks F101 and
F102 as they have a vapor and a liquid outlet.

.. code:: ipython3

    m.fs.s05 = Arc(source=m.fs.R101.outlet, destination=m.fs.F101.inlet)
    m.fs.s06 = Arc(source=m.fs.F101.vap_outlet, destination=m.fs.S101.inlet)
    m.fs.s08 = Arc(source=m.fs.S101.recycle, destination=m.fs.C101.inlet)
    m.fs.s09 = Arc(source=m.fs.C101.outlet,
                   destination=m.fs.M101.vapor_recycle)
    m.fs.s10 = Arc(source=m.fs.F101.liq_outlet, destination=m.fs.F102.inlet)

We have now connected the unit model block using the arcs. However, each
of these arcs link to ports on the two unit models that are connected.
In this case, the ports consist of the state variables that need to be
linked between the unit models. Pyomo provides a convenient method to
write these equality constraints for us between two ports and this is
done as follows:

.. code:: ipython3

    TransformationFactory("network.expand_arcs").apply_to(m)

Adding expressions to compute purity and operating costs
--------------------------------------------------------

In this section, we will add a few Expressions that allows us to
evaluate the performance. Expressions provide a convenient way of
calculating certain values that are a function of the variables defined
in the model. For more details on Expressions, please refer to:
https://pyomo.readthedocs.io/en/latest/pyomo_modeling_components/Expressions.html

For this flowsheet, we are interested in computing the purity of the
product Benzene stream (i.e. the mole fraction) and the operating cost
which is a sum of the cooling and heating cost.

Let us first add an Expression to compute the mole fraction of benzene
in the ``vap_outlet`` of F102 which is our product stream. Please note
that the var flow_mol_phase_comp has the index - [time, phase,
component]. As this is a steady-state flowsheet, the time index by
default is 0. The valid phases are [“Liq”, “Vap”]. Similarly the valid
component list is [“benzene”, “toluene”, “hydrogen”, “methane”].

.. code:: ipython3

    m.fs.purity = Expression(
            expr=m.fs.F102.vap_outlet.flow_mol_phase_comp[0, "Vap", "benzene"] /
            (m.fs.F102.vap_outlet.flow_mol_phase_comp[0, "Vap", "benzene"]
             + m.fs.F102.vap_outlet.flow_mol_phase_comp[0, "Vap", "toluene"]))

Now, let us add an expression to compute the cooling cost assuming a
cost of 0.212E-4 $/kW. Note that cooling utility is required for the
reactor (R101) and the first flash (F101).

.. code:: ipython3

    m.fs.cooling_cost = Expression(expr=0.212e-7 * (-m.fs.F101.heat_duty[0]) +
                                       0.212e-7 * (-m.fs.R101.heat_duty[0]))

Now, let us add an expression to compute the heating cost assuming the
utility cost as follows:

.. raw:: html

   <ul>

.. raw:: html

   <li>

2.2E-4 dollars/kW for H101

.. raw:: html

   </li>

.. raw:: html

   <li>

1.9E-4 dollars/kW for F102

.. raw:: html

   </li>

.. raw:: html

   </ul>

Note that the heat duty is in units of watt (J/s).

.. code:: ipython3

    m.fs.heating_cost = Expression(expr=2.2e-7 * m.fs.H101.heat_duty[0] +
                                       1.9e-7 * m.fs.F102.heat_duty[0])

Let us now add an expression to compute the total operating cost per
year which is basically the sum of the cooling and heating cost we
defined above.

.. code:: ipython3

    m.fs.operating_cost = Expression(expr=(3600 * 24 * 365 *
                                               (m.fs.heating_cost +
                                                m.fs.cooling_cost)))

Fixing feed conditions
----------------------

Let us first check how many degrees of freedom exist for this flowsheet
using the ``degrees_of_freedom`` tool we imported earlier.

.. code:: ipython3

    print(degrees_of_freedom(m))


.. parsed-literal::

    49


We will now be fixing the toluene feed stream to the conditions shown in
the flowsheet above. Please note that though this is a pure toluene
feed, the remaining components are still assigned a very small non-zero
value to help with convergence and initializing.

.. code:: ipython3

    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "benzene"].fix(1e-5)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "toluene"].fix(1e-5)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "hydrogen"].fix(1e-5)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "methane"].fix(1e-5)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "benzene"].fix(1e-5)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "toluene"].fix(0.30)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "hydrogen"].fix(1e-5)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "methane"].fix(1e-5)
    m.fs.M101.toluene_feed.temperature.fix(303.2)
    m.fs.M101.toluene_feed.pressure.fix(350000)

Similarly, let us fix the hydrogen feed to the following conditions in
the next cell:

.. raw:: html

   <ul>

.. raw:: html

   <li>

FH2 = 0.30 mol/s

.. raw:: html

   </li>

.. raw:: html

   <li>

FCH4 = 0.02 mol/s

.. raw:: html

   </li>

.. raw:: html

   <li>

Remaining components = 1e-5 mol/s

.. raw:: html

   </li>

.. raw:: html

   <li>

T = 303.2 K

.. raw:: html

   </li>

.. raw:: html

   <li>

P = 350000 Pa

.. raw:: html

   </li>

.. raw:: html

   </ul>

.. code:: ipython3

    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "benzene"].fix(1e-5)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "toluene"].fix(1e-5)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "hydrogen"].fix(0.30)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "methane"].fix(0.02)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "benzene"].fix(1e-5)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "toluene"].fix(1e-5)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "hydrogen"].fix(1e-5)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "methane"].fix(1e-5)
    m.fs.M101.hydrogen_feed.temperature.fix(303.2)
    m.fs.M101.hydrogen_feed.pressure.fix(350000)

Fixing unit model specifications
--------------------------------

Now that we have fixed our inlet feed conditions, we will now be fixing
the operating conditions for the unit models in the flowsheet. Let us
set set the H101 outlet temperature to 600 K.

.. code:: ipython3

    m.fs.H101.outlet.temperature.fix(600)

For the StoichiometricReactor, we have to define the conversion in terms
of toluene. This requires us to create a new variable for specifying the
conversion and adding a Constraint that defines the conversion with
respect to toluene. The second degree of freedom for the reactor is to
define the heat duty. In this case, let us assume the reactor to be
adiabatic i.e. Q = 0.

.. code:: ipython3

    m.fs.R101.conversion = Var(initialize=0.75, bounds=(0, 1))
    
    m.fs.R101.conv_constraint = Constraint(
        expr=m.fs.R101.conversion*m.fs.R101.inlet.
        flow_mol_phase_comp[0, "Vap", "toluene"] ==
        (m.fs.R101.inlet.flow_mol_phase_comp[0, "Vap", "toluene"] -
         m.fs.R101.outlet.flow_mol_phase_comp[0, "Vap", "toluene"]))
    
    m.fs.R101.conversion.fix(0.75)
    m.fs.R101.heat_duty.fix(0)

The Flash conditions for F101 can be set as follows.

.. code:: ipython3

    m.fs.F101.vap_outlet.temperature.fix(325.0)
    m.fs.F101.deltaP.fix(0)

.. raw:: html

   <div class="alert alert-block alert-info">

Inline Exercise: Set the conditions for Flash F102 to the following
conditions:

.. raw:: html

   <ul>

.. raw:: html

   <li>

T = 375 K

.. raw:: html

   </li>

.. raw:: html

   <li>

deltaP = -200000

.. raw:: html

   </li>

.. raw:: html

   </ul>

Use Shift+Enter to run the cell once you have typed in your code.

.. raw:: html

   </div>

.. code:: ipython3

    m.fs.F102.vap_outlet.temperature.fix(375)
    m.fs.F102.deltaP.fix(-200000)

Let us fix the purge split fraction to 20% and the outlet pressure of
the compressor is set to 350000 Pa.

.. code:: ipython3

    m.fs.S101.split_fraction[0, "purge"].fix(0.2)
    m.fs.C101.outlet.pressure.fix(350000)

.. raw:: html

   <div class="alert alert-block alert-info">

Inline Exercise: We have now defined all the feed conditions and the
inputs required for the unit models. The system should now have 0
degrees of freedom i.e. should be a square problem. Please check that
the degrees of freedom is 0.

Use Shift+Enter to run the cell once you have typed in your code.

.. raw:: html

   </div>

.. code:: ipython3

    print(degrees_of_freedom(m))


.. parsed-literal::

    20


Initialization
--------------

This section will demonstrate how to use the built-in sequential
decomposition tool to initialize our flowsheet.

|image0|

.. |image0| image:: module_2_flowsheet.png

Let us first create an object for the SequentialDecomposition and
specify our options for this.

.. code:: ipython3

    seq = SequentialDecomposition()
    seq.options.select_tear_method = "heuristic"
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = 5
    
    # Using the SD tool
    G = seq.create_graph(m)
    heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
    order = seq.calculation_order(G)

Which is the tear stream? Display tear set and order

.. code:: ipython3

    for o in heuristic_tear_set:
        print(o.name)


.. parsed-literal::

    fs.s03


What sequence did the SD tool determine to solve this flowsheet with the
least number of tears?

.. code:: ipython3

    for o in order:
        print(o[0].name)


.. parsed-literal::

    fs.H101
    fs.R101
    fs.F101
    fs.S101
    fs.C101
    fs.M101


|image0|

The SequentialDecomposition tool has determined that the tear stream is
the mixer outlet. We will need to provide a reasonable guess for this.

.. |image0| image:: module_2_tear_stream.png

.. code:: ipython3

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
    
    # Pass the tear_guess to the SD tool
    seq.set_guesses_for(m.fs.H101.inlet, tear_guesses)

Next, we need to tell the tool how to initialize a particular unit. We
will be writing a python function which takes in a “unit” and calls the
initialize method on that unit.

.. code:: ipython3

    def function(unit):
            unit.initialize(outlvl=1)

We are now ready to initialize our flowsheet in a sequential mode. Note
that we specifically set the iteration limit to be 5 as we are trying to
use this tool only to get a good set of initial values such that IPOPT
can then take over and solve this flowsheet for us.

.. code:: ipython3

    seq.run(m, function)


.. parsed-literal::

    2019-11-08 04:02:45 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 1 Complete.
    2019-11-08 04:02:45 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 2 Complete.
    2019-11-08 04:02:45 - INFO - idaes.core.unit_model - fs.H101 Initialisation Complete.
    2019-11-08 04:02:45 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 1 Complete.
    2019-11-08 04:02:46 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 2 Complete.
    2019-11-08 04:02:46 - INFO - idaes.core.unit_model - fs.R101 Initialisation Complete.
    2019-11-08 04:02:46 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 1 Complete.
    2019-11-08 04:02:46 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 2 Complete.
    2019-11-08 04:02:46 - INFO - idaes.core.unit_model - fs.F101 Initialisation Complete.
    2019-11-08 04:02:46 - INFO - idaes.unit_models.separator - fs.S101 Initialisation Complete.
    2019-11-08 04:02:46 - INFO - idaes.core.unit_model - fs.F102 Initialisation Step 1 Complete.
    2019-11-08 04:02:46 - INFO - idaes.core.unit_model - fs.F102 Initialisation Step 2 Complete.
    2019-11-08 04:02:46 - INFO - idaes.core.unit_model - fs.F102 Initialisation Complete.
    2019-11-08 04:02:46 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 1 Complete.
    2019-11-08 04:02:46 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 2 Complete.
    2019-11-08 04:02:46 - INFO - idaes.core.unit_model - fs.C101 Initialisation Complete.
    2019-11-08 04:02:46 - INFO - idaes.unit_models.mixer - fs.M101 Initialisation Complete.
    2019-11-08 04:02:46 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 1 Complete.
    2019-11-08 04:02:46 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 2 Complete.
    2019-11-08 04:02:46 - INFO - idaes.core.unit_model - fs.H101 Initialisation Complete.
    2019-11-08 04:02:46 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 1 Complete.
    2019-11-08 04:02:46 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 2 Complete.
    2019-11-08 04:02:46 - INFO - idaes.core.unit_model - fs.R101 Initialisation Complete.
    2019-11-08 04:02:46 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 1 Complete.
    2019-11-08 04:02:46 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 2 Complete.
    2019-11-08 04:02:46 - INFO - idaes.core.unit_model - fs.F101 Initialisation Complete.
    2019-11-08 04:02:46 - INFO - idaes.unit_models.separator - fs.S101 Initialisation Complete.
    2019-11-08 04:02:46 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 1 Complete.
    2019-11-08 04:02:46 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 2 Complete.
    2019-11-08 04:02:46 - INFO - idaes.core.unit_model - fs.C101 Initialisation Complete.
    2019-11-08 04:02:47 - INFO - idaes.unit_models.mixer - fs.M101 Initialisation Complete.
    2019-11-08 04:02:47 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 1 Complete.
    2019-11-08 04:02:47 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 2 Complete.
    2019-11-08 04:02:47 - INFO - idaes.core.unit_model - fs.H101 Initialisation Complete.
    2019-11-08 04:02:47 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 1 Complete.
    2019-11-08 04:02:47 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 2 Complete.
    2019-11-08 04:02:47 - INFO - idaes.core.unit_model - fs.R101 Initialisation Complete.
    2019-11-08 04:02:47 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 1 Complete.
    2019-11-08 04:02:47 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 2 Complete.
    2019-11-08 04:02:47 - INFO - idaes.core.unit_model - fs.F101 Initialisation Complete.
    2019-11-08 04:02:47 - INFO - idaes.unit_models.separator - fs.S101 Initialisation Complete.
    2019-11-08 04:02:47 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 1 Complete.
    2019-11-08 04:02:47 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 2 Complete.
    2019-11-08 04:02:47 - INFO - idaes.core.unit_model - fs.C101 Initialisation Complete.
    2019-11-08 04:02:47 - INFO - idaes.unit_models.mixer - fs.M101 Initialisation Complete.
    2019-11-08 04:02:47 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 1 Complete.
    2019-11-08 04:02:47 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 2 Complete.
    2019-11-08 04:02:47 - INFO - idaes.core.unit_model - fs.H101 Initialisation Complete.
    2019-11-08 04:02:47 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 1 Complete.
    2019-11-08 04:02:47 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 2 Complete.
    2019-11-08 04:02:47 - INFO - idaes.core.unit_model - fs.R101 Initialisation Complete.
    2019-11-08 04:02:47 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 1 Complete.
    2019-11-08 04:02:47 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 2 Complete.
    2019-11-08 04:02:47 - INFO - idaes.core.unit_model - fs.F101 Initialisation Complete.
    2019-11-08 04:02:47 - INFO - idaes.unit_models.separator - fs.S101 Initialisation Complete.
    2019-11-08 04:02:47 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 1 Complete.
    2019-11-08 04:02:47 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 2 Complete.
    2019-11-08 04:02:47 - INFO - idaes.core.unit_model - fs.C101 Initialisation Complete.
    2019-11-08 04:02:48 - INFO - idaes.unit_models.mixer - fs.M101 Initialisation Complete.
    2019-11-08 04:02:48 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 1 Complete.
    2019-11-08 04:02:48 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 2 Complete.
    2019-11-08 04:02:48 - INFO - idaes.core.unit_model - fs.H101 Initialisation Complete.
    2019-11-08 04:02:48 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 1 Complete.
    2019-11-08 04:02:48 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 2 Complete.
    2019-11-08 04:02:48 - INFO - idaes.core.unit_model - fs.R101 Initialisation Complete.
    2019-11-08 04:02:48 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 1 Complete.
    2019-11-08 04:02:48 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 2 Complete.
    2019-11-08 04:02:48 - INFO - idaes.core.unit_model - fs.F101 Initialisation Complete.
    2019-11-08 04:02:48 - INFO - idaes.unit_models.separator - fs.S101 Initialisation Complete.
    2019-11-08 04:02:48 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 1 Complete.
    2019-11-08 04:02:48 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 2 Complete.
    2019-11-08 04:02:48 - INFO - idaes.core.unit_model - fs.C101 Initialisation Complete.
    2019-11-08 04:02:48 - INFO - idaes.unit_models.mixer - fs.M101 Initialisation Complete.
    2019-11-08 04:02:48 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 1 Complete.
    2019-11-08 04:02:48 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 2 Complete.
    2019-11-08 04:02:48 - INFO - idaes.core.unit_model - fs.H101 Initialisation Complete.
    2019-11-08 04:02:48 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 1 Complete.
    2019-11-08 04:02:48 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 2 Complete.
    2019-11-08 04:02:48 - INFO - idaes.core.unit_model - fs.R101 Initialisation Complete.
    2019-11-08 04:02:48 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 1 Complete.
    2019-11-08 04:02:48 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 2 Complete.
    2019-11-08 04:02:48 - INFO - idaes.core.unit_model - fs.F101 Initialisation Complete.
    2019-11-08 04:02:48 - INFO - idaes.unit_models.separator - fs.S101 Initialisation Complete.
    2019-11-08 04:02:49 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 1 Complete.
    2019-11-08 04:02:49 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 2 Complete.
    2019-11-08 04:02:49 - INFO - idaes.core.unit_model - fs.C101 Initialisation Complete.
    2019-11-08 04:02:49 - INFO - idaes.unit_models.mixer - fs.M101 Initialisation Complete.
    2019-11-08 04:02:49 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 1 Complete.
    2019-11-08 04:02:49 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 2 Complete.
    2019-11-08 04:02:49 - INFO - idaes.core.unit_model - fs.H101 Initialisation Complete.
    2019-11-08 04:02:49 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 1 Complete.
    2019-11-08 04:02:49 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 2 Complete.
    2019-11-08 04:02:49 - INFO - idaes.core.unit_model - fs.R101 Initialisation Complete.
    2019-11-08 04:02:49 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 1 Complete.
    2019-11-08 04:02:49 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 2 Complete.
    2019-11-08 04:02:49 - INFO - idaes.core.unit_model - fs.F101 Initialisation Complete.
    2019-11-08 04:02:49 - INFO - idaes.unit_models.separator - fs.S101 Initialisation Complete.
    2019-11-08 04:02:49 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 1 Complete.
    2019-11-08 04:02:49 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 2 Complete.
    2019-11-08 04:02:49 - INFO - idaes.core.unit_model - fs.C101 Initialisation Complete.
    2019-11-08 04:02:49 - INFO - idaes.unit_models.mixer - fs.M101 Initialisation Complete.
    WARNING: Wegstein failed to converge in 5 iterations
    2019-11-08 04:02:49 - INFO - idaes.core.unit_model - fs.F102 Initialisation Step 1 Complete.
    2019-11-08 04:02:49 - INFO - idaes.core.unit_model - fs.F102 Initialisation Step 2 Complete.
    2019-11-08 04:02:49 - INFO - idaes.core.unit_model - fs.F102 Initialisation Complete.


.. raw:: html

   <div class="alert alert-block alert-info">

Inline Exercise: We have now initialized the flowsheet. Let us run the
flowsheet in a simulation mode to look at the results. To do this,
complete the last line of code where we pass the model to the solver.
You will need to type the following:

results = solver.solve(m, tee=True)

Use Shift+Enter to run the cell once you have typed in your code.

.. raw:: html

   </div>

.. code:: ipython3

    # Create the solver object
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6, 'max_iter': 5000}
    
    # Solve the model
    results = solver.solve(m, tee=True)


.. parsed-literal::

    Ipopt 3.12.8: tol=1e-06
    max_iter=5000
    
    
    ******************************************************************************
    This program contains Ipopt, a library for large-scale nonlinear optimization.
     Ipopt is released as open source code under the Eclipse Public License (EPL).
             For more information visit http://projects.coin-or.org/Ipopt
    ******************************************************************************
    
    This is Ipopt version 3.12.8, running with linear solver mumps.
    NOTE: Other linear solvers might be more efficient (see Ipopt documentation).
    
    Number of nonzeros in equality constraint Jacobian...:      994
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:      898
    
    Total number of variables............................:      326
                         variables with only lower bounds:        0
                    variables with lower and upper bounds:      146
                         variables with only upper bounds:        0
    Total number of equality constraints.................:      306
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  0.0000000e+00 4.21e+04 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  0.0000000e+00 3.05e+04 1.21e+02  -1.0 2.19e+04    -  9.15e-01 1.86e-01h  1
       2  0.0000000e+00 1.41e+04 5.73e+02  -1.0 1.60e+04    -  9.89e-01 4.77e-01h  1
       3  0.0000000e+00 1.38e+04 9.90e+04  -1.0 1.57e+04    -  9.90e-01 1.55e-02h  6
       4  0.0000000e+00 1.37e+04 3.00e+05  -1.0 1.32e+04    -  1.00e+00 7.75e-03h  7
       5  0.0000000e+00 1.36e+04 7.03e+05  -1.0 1.16e+04    -  1.00e+00 7.76e-03h  7
       6  0.0000000e+00 1.36e+04 1.51e+06  -1.0 1.27e+04    -  1.00e+00 3.88e-03h  8
       7  0.0000000e+00 1.35e+04 3.11e+06  -1.0 1.46e+04    -  1.00e+00 3.89e-03h  8
       8  0.0000000e+00 1.35e+04 6.28e+06  -1.0 1.77e+04    -  1.00e+00 1.95e-03h  9
       9  0.0000000e+00 1.35e+04 1.25e+07  -1.0 2.88e+04    -  1.00e+00 1.97e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  0.0000000e+00 1.35e+04 2.44e+07  -1.0 6.01e+04    -  1.00e+00 9.94e-04h 10
      11  0.0000000e+00 1.34e+04 4.69e+07  -1.0 9.93e+04    -  1.00e+00 4.83e-04h 11
      12  0.0000000e+00 1.34e+04 8.93e+07  -1.0 1.24e+05    -  1.00e+00 4.07e-04h 11
      13  0.0000000e+00 7.48e+05 1.01e+08  -1.0 1.28e+05    -  1.00e+00 4.07e-01w  1
      14  0.0000000e+00 3.93e+05 1.91e+08  -1.0 2.68e+04    -  6.78e-01 4.91e-01w  1
      15  0.0000000e+00 2.54e+05 1.96e+08  -1.0 6.63e+04    -  4.89e-01 4.37e-01w  1
      16  0.0000000e+00 1.34e+04 1.70e+08  -1.0 1.77e+04    -  1.00e+00 3.97e-04h 10
      17  0.0000000e+00 1.34e+04 3.23e+08  -1.0 1.23e+05    -  1.00e+00 4.08e-04h 11
      18  0.0000000e+00 1.34e+04 6.15e+08  -1.0 1.16e+05    -  1.00e+00 4.29e-04h 11
      19  0.0000000e+00 1.34e+04 1.12e+09  -1.0 1.08e+05    -  9.02e-01 4.51e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      20  0.0000000e+00 1.34e+04 2.15e+09  -1.0 1.03e+05    -  1.00e+00 4.68e-04h 11
      21  0.0000000e+00 1.34e+04 2.92e+09  -1.0 9.90e+04    -  3.93e-01 4.81e-04h 11
      22  0.0000000e+00 1.34e+04 5.60e+09  -1.0 9.75e+04    -  1.00e+00 4.86e-04h 11
      23  0.0000000e+00 1.34e+04 7.69e+09  -1.0 9.56e+04    -  4.07e-01 4.93e-04h 11
      24  0.0000000e+00 1.34e+04 1.48e+10  -1.0 9.48e+04    -  1.00e+00 4.96e-04h 11
      25  0.0000000e+00 1.34e+04 2.04e+10  -1.0 9.39e+04    -  4.13e-01 4.99e-04h 11
      26  0.0000000e+00 8.25e+05 1.91e+10  -1.0 9.34e+04    -  1.00e+00 5.12e-01w  1
      27  0.0000000e+00 4.03e+05 4.76e+10  -1.0 2.26e+04    -  6.55e-01 5.20e-01w  1
      28  0.0000000e+00 1.09e+05 1.49e+12  -1.0 1.13e+04    -  1.97e-01 7.36e-01w  1
      29  0.0000000e+00 1.34e+04 3.92e+10  -1.0 3.98e+03    -  1.00e+00 5.00e-04h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      30  0.0000000e+00 1.34e+04 5.42e+10  -1.0 9.30e+04    -  4.16e-01 5.01e-04h 11
      31  0.0000000e+00 1.34e+04 1.04e+11  -1.0 9.27e+04    -  1.00e+00 5.02e-04h 11
      32  0.0000000e+00 1.34e+04 1.44e+11  -1.0 9.25e+04    -  4.18e-01 5.02e-04h 11
      33  0.0000000e+00 1.33e+04 2.77e+11  -1.0 9.22e+04    -  1.00e+00 5.03e-04h 11
      34  0.0000000e+00 1.33e+04 3.84e+11  -1.0 9.21e+04    -  4.19e-01 5.03e-04h 11
      35  0.0000000e+00 1.33e+04 7.39e+11  -1.0 9.19e+04    -  1.00e+00 5.03e-04h 11
      36  0.0000000e+00 1.33e+04 1.02e+12  -1.0 9.19e+04    -  4.20e-01 5.03e-04h 11
      37  0.0000000e+00 1.33e+04 1.97e+12  -1.0 9.16e+04    -  1.00e+00 5.03e-04h 11
      38  0.0000000e+00 1.33e+04 2.74e+12  -1.0 9.16e+04    -  4.21e-01 5.03e-04h 11
      39  0.0000000e+00 8.13e+05 2.55e+12  -1.0 9.14e+04    -  1.00e+00 5.15e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      40  0.0000000e+00 3.78e+05 6.75e+12  -1.0 2.25e+04    -  6.54e-01 5.43e-01w  1
      41  0.0000000e+00 2.66e+05 2.45e+13  -1.0 7.75e+03    -  1.20e-01 2.98e-01w  1
      42  0.0000000e+00 1.33e+04 5.26e+12  -1.0 7.80e+03    -  1.00e+00 5.03e-04h 10
      43  0.0000000e+00 1.33e+04 7.30e+12  -1.0 9.14e+04    -  4.22e-01 5.03e-04h 11
      44  0.0000000e+00 1.33e+04 1.40e+13  -1.0 9.12e+04    -  1.00e+00 5.03e-04h 11
      45  0.0000000e+00 1.33e+04 1.95e+13  -1.0 9.12e+04    -  4.22e-01 5.03e-04h 11
      46  0.0000000e+00 1.33e+04 3.75e+13  -1.0 9.10e+04    -  1.00e+00 5.03e-04h 11
      47  0.0000000e+00 1.33e+04 5.22e+13  -1.0 9.10e+04    -  4.23e-01 5.03e-04h 11
      48  0.0000000e+00 1.33e+04 1.00e+14  -1.0 9.08e+04    -  1.00e+00 5.03e-04h 11
      49  0.0000000e+00 1.33e+04 1.40e+14  -1.0 9.08e+04    -  4.24e-01 5.03e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      50  0.0000000e+00 1.33e+04 2.68e+14  -1.0 9.06e+04    -  1.00e+00 5.03e-04h 11
      51  0.0000000e+00 1.32e+04 3.73e+14  -1.0 9.06e+04    -  4.24e-01 5.03e-04h 11
      52  0.0000000e+00 7.98e+05 3.49e+14  -1.0 9.04e+04    -  1.00e+00 5.15e-01w  1
      53  0.0000000e+00 3.71e+05 9.23e+14  -1.0 2.24e+04    -  6.55e-01 5.43e-01w  1
      54  0.0000000e+00 2.62e+05 3.28e+15  -1.0 7.65e+03    -  1.19e-01 2.95e-01w  1
      55  0.0000000e+00 1.32e+04 7.18e+14  -1.0 7.78e+03    -  1.00e+00 5.03e-04h 10
      56  0.0000000e+00 1.32e+04 1.00e+15  -1.0 9.03e+04    -  4.25e-01 5.03e-04h 11
      57  0.0000000e+00 1.32e+04 1.13e+15  -1.0 9.01e+04    -  1.00e+00 5.03e-04h 11
      58  0.0000000e+00 1.32e+04 1.13e+15  -1.0 3.96e+04    -  8.00e-01 1.97e-03h  9
      59  0.0000000e+00 1.32e+04 1.13e+15  -1.0 3.18e+04    -  1.00e+00 1.97e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      60  0.0000000e+00 1.32e+04 1.13e+15  -1.0 2.98e+04    -  9.52e-01 1.97e-03h  9
      61  0.0000000e+00 1.31e+04 1.13e+15  -1.0 2.97e+04    -  1.00e+00 1.97e-03h  9
      62  0.0000000e+00 1.31e+04 1.13e+15  -1.0 2.96e+04    -  9.56e-01 1.97e-03h  9
      63  0.0000000e+00 1.31e+04 1.13e+15  -1.0 2.95e+04    -  1.00e+00 1.97e-03h  9
      64  0.0000000e+00 1.30e+04 1.13e+15  -1.0 2.95e+04    -  9.60e-01 1.97e-03h  9
      65  0.0000000e+00 2.85e+05 1.11e+15  -1.0 2.94e+04    -  1.00e+00 5.04e-01w  1
      66  0.0000000e+00 1.32e+05 2.90e+15  -1.0 3.75e+04    -  6.63e-01 5.44e-01w  1
      67  0.0000000e+00 4.04e+05 4.78e+16  -1.0 6.01e+04    -  1.19e-01 5.55e-01w  1
      68  0.0000000e+00 1.30e+04 1.13e+15  -1.0 4.29e+03    -  1.00e+00 1.97e-03h  8
      69  0.0000000e+00 1.30e+04 1.13e+15  -1.0 2.93e+04    -  9.64e-01 1.97e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      70  0.0000000e+00 1.30e+04 1.13e+15  -1.0 2.92e+04    -  1.00e+00 1.97e-03h  9
      71  0.0000000e+00 1.29e+04 1.13e+15  -1.0 2.91e+04    -  9.68e-01 1.97e-03h  9
      72  0.0000000e+00 1.29e+04 1.13e+15  -1.0 2.91e+04    -  1.00e+00 1.97e-03h  9
      73  0.0000000e+00 1.29e+04 1.13e+15  -1.0 2.90e+04    -  9.72e-01 1.97e-03h  9
      74  0.0000000e+00 1.29e+04 1.13e+15  -1.0 2.89e+04    -  1.00e+00 1.97e-03h  9
      75  0.0000000e+00 1.28e+04 1.13e+15  -1.0 2.88e+04    -  9.76e-01 1.97e-03h  9
      76  0.0000000e+00 1.28e+04 1.13e+15  -1.0 2.87e+04    -  1.00e+00 1.97e-03h  9
      77  0.0000000e+00 1.28e+04 1.13e+15  -1.0 2.87e+04    -  9.80e-01 1.97e-03h  9
      78  0.0000000e+00 2.81e+05 1.14e+15  -1.0 2.86e+04    -  1.00e+00 5.04e-01w  1
      79  0.0000000e+00 1.30e+05 3.02e+15  -1.0 3.70e+04    -  6.59e-01 5.45e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      80  0.0000000e+00 3.85e+05 4.73e+16  -1.0 6.10e+04    -  1.17e-01 5.42e-01w  1
      81  0.0000000e+00 1.28e+04 1.13e+15  -1.0 4.29e+03    -  1.00e+00 1.97e-03h  8
      82  0.0000000e+00 1.27e+04 1.13e+15  -1.0 2.85e+04    -  9.84e-01 1.97e-03h  9
      83  0.0000000e+00 1.27e+04 1.13e+15  -1.0 2.84e+04    -  1.00e+00 1.97e-03h  9
      84  0.0000000e+00 1.27e+04 1.13e+15  -1.0 2.84e+04    -  9.88e-01 1.97e-03h  9
      85  0.0000000e+00 1.27e+04 1.13e+15  -1.0 2.83e+04    -  1.00e+00 1.97e-03h  9
      86  0.0000000e+00 1.26e+04 1.13e+15  -1.0 2.82e+04    -  9.92e-01 1.97e-03h  9
      87  0.0000000e+00 1.26e+04 1.13e+15  -1.0 2.81e+04    -  1.00e+00 1.97e-03h  9
      88  0.0000000e+00 1.26e+04 1.13e+15  -1.0 2.80e+04    -  9.95e-01 1.97e-03h  9
      89  0.0000000e+00 1.26e+04 1.13e+15  -1.0 2.80e+04    -  1.00e+00 1.97e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      90  0.0000000e+00 1.25e+04 1.13e+15  -1.0 2.79e+04    -  9.99e-01 1.97e-03h  9
      91  0.0000000e+00 2.69e+05 1.16e+15  -1.0 2.78e+04    -  1.00e+00 5.04e-01w  1
      92  0.0000000e+00 1.24e+05 3.14e+15  -1.0 3.66e+04    -  6.55e-01 5.47e-01w  1
      93  0.0000000e+00 3.72e+05 4.69e+16  -1.0 6.18e+04    -  1.14e-01 5.29e-01w  1
      94  0.0000000e+00 1.25e+04 1.13e+15  -1.0 4.29e+03    -  1.00e+00 1.97e-03h  8
      95  0.0000000e+00 1.25e+04 1.13e+15  -1.0 2.77e+04    -  1.00e+00 1.97e-03h  9
      96  0.0000000e+00 1.25e+04 1.13e+15  -1.0 2.77e+04    -  1.00e+00 1.97e-03h  9
      97  0.0000000e+00 1.24e+04 1.13e+15  -1.0 2.76e+04    -  1.00e+00 1.97e-03h  9
      98  0.0000000e+00 1.24e+04 1.13e+15  -1.0 2.75e+04    -  1.00e+00 1.97e-03h  9
      99  0.0000000e+00 1.24e+04 1.13e+15  -1.0 2.74e+04    -  1.00e+00 1.97e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     100  0.0000000e+00 1.24e+04 1.13e+15  -1.0 2.74e+04    -  1.00e+00 1.97e-03h  9
     101  0.0000000e+00 1.24e+04 1.13e+15  -1.0 2.73e+04    -  1.00e+00 1.97e-03h  9
     102  0.0000000e+00 1.23e+04 1.13e+15  -1.0 2.72e+04    -  1.00e+00 1.97e-03h  9
     103  0.0000000e+00 1.23e+04 1.13e+15  -1.0 2.71e+04    -  1.00e+00 1.97e-03h  9
     104  0.0000000e+00 2.58e+05 1.19e+15  -1.0 2.71e+04    -  1.00e+00 5.03e-01w  1
     105  0.0000000e+00 1.19e+05 3.27e+15  -1.0 3.62e+04    -  6.51e-01 5.48e-01w  1
     106  0.0000000e+00 3.59e+05 4.72e+16  -1.0 6.23e+04    -  1.12e-01 5.20e-01w  1
     107  0.0000000e+00 1.23e+04 1.13e+15  -1.0 4.27e+03    -  1.00e+00 1.97e-03h  8
     108  0.0000000e+00 1.23e+04 1.13e+15  -1.0 2.70e+04    -  1.00e+00 1.97e-03h  9
     109  0.0000000e+00 1.22e+04 1.13e+15  -1.0 2.69e+04    -  1.00e+00 1.97e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     110  0.0000000e+00 1.22e+04 1.13e+15  -1.0 2.68e+04    -  1.00e+00 1.97e-03h  9
     111  0.0000000e+00 1.22e+04 1.13e+15  -1.0 2.68e+04    -  1.00e+00 1.97e-03h  9
     112  0.0000000e+00 1.22e+04 1.13e+15  -1.0 2.67e+04    -  1.00e+00 1.97e-03h  9
     113  0.0000000e+00 1.21e+04 1.13e+15  -1.0 2.66e+04    -  1.00e+00 1.97e-03h  9
     114  0.0000000e+00 1.21e+04 1.14e+15  -1.0 2.66e+04    -  1.00e+00 1.97e-03h  9
     115  0.0000000e+00 1.21e+04 1.14e+15  -1.0 2.65e+04    -  1.00e+00 1.97e-03h  9
     116  0.0000000e+00 1.21e+04 1.14e+15  -1.0 2.64e+04    -  1.00e+00 1.97e-03h  9
     117  0.0000000e+00 2.48e+05 1.21e+15  -1.0 2.63e+04    -  1.00e+00 5.03e-01w  1
     118  0.0000000e+00 1.14e+05 3.41e+15  -1.0 3.58e+04    -  6.47e-01 5.50e-01w  1
     119  0.0000000e+00 3.46e+05 5.02e+16  -1.0 6.13e+04    -  1.10e-01 5.21e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     120  0.0000000e+00 1.20e+04 1.14e+15  -1.0 4.22e+03    -  1.00e+00 1.97e-03h  8
     121  0.0000000e+00 1.20e+04 1.14e+15  -1.0 2.63e+04    -  1.00e+00 1.97e-03h  9
     122  0.0000000e+00 1.20e+04 1.14e+15  -1.0 2.62e+04    -  1.00e+00 1.97e-03h  9
     123  0.0000000e+00 1.20e+04 1.14e+15  -1.0 2.61e+04    -  1.00e+00 1.97e-03h  9
     124  0.0000000e+00 1.19e+04 1.14e+15  -1.0 2.61e+04    -  1.00e+00 1.97e-03h  9
     125  0.0000000e+00 1.19e+04 1.14e+15  -1.0 2.60e+04    -  1.00e+00 1.97e-03h  9
     126  0.0000000e+00 1.19e+04 1.14e+15  -1.0 2.59e+04    -  1.00e+00 1.97e-03h  9
     127  0.0000000e+00 1.19e+04 1.14e+15  -1.0 2.58e+04    -  1.00e+00 1.97e-03h  9
     128  0.0000000e+00 1.19e+04 1.14e+15  -1.0 2.58e+04    -  1.00e+00 1.97e-03h  9
     129  0.0000000e+00 1.18e+04 1.14e+15  -1.0 2.57e+04    -  1.00e+00 1.97e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     130  0.0000000e+00 2.38e+05 1.24e+15  -1.0 2.56e+04    -  1.00e+00 5.03e-01w  1
     131  0.0000000e+00 1.09e+05 3.56e+15  -1.0 3.54e+04    -  6.43e-01 5.51e-01w  1
     132  0.0000000e+00 3.31e+05 6.54e+16  -1.0 5.54e+04    -  1.08e-01 5.63e-01w  1
     133  0.0000000e+00 1.18e+04 1.14e+15  -1.0 4.06e+03    -  1.00e+00 1.97e-03h  8
     134  0.0000000e+00 1.18e+04 1.14e+15  -1.0 2.55e+04    -  1.00e+00 1.97e-03h  9
     135  0.0000000e+00 1.18e+04 1.14e+15  -1.0 2.55e+04    -  1.00e+00 1.97e-03h  9
     136  0.0000000e+00 1.17e+04 1.14e+15  -1.0 2.54e+04    -  1.00e+00 1.97e-03h  9
     137  0.0000000e+00 1.17e+04 1.14e+15  -1.0 2.53e+04    -  1.00e+00 1.97e-03h  9
     138  0.0000000e+00 1.17e+04 1.14e+15  -1.0 2.53e+04    -  1.00e+00 1.97e-03h  9
     139  0.0000000e+00 1.17e+04 1.15e+15  -1.0 2.52e+04    -  1.00e+00 1.97e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     140  0.0000000e+00 1.16e+04 1.15e+15  -1.0 2.51e+04    -  1.00e+00 1.97e-03h  9
     141  0.0000000e+00 1.16e+04 1.15e+15  -1.0 2.50e+04    -  1.00e+00 1.97e-03h  9
     142  0.0000000e+00 1.16e+04 1.15e+15  -1.0 2.50e+04    -  1.00e+00 1.97e-03h  9
     143  0.0000000e+00 2.28e+05 1.27e+15  -1.0 2.49e+04    -  1.00e+00 5.03e-01w  1
     144  0.0000000e+00 1.04e+05 3.72e+15  -1.0 3.52e+04    -  6.37e-01 5.53e-01w  1
     145  0.0000000e+00 3.08e+05 1.85e+17  -1.0 3.86e+04    -  1.07e-01 7.55e-01w  1
     146  0.0000000e+00 1.16e+04 1.15e+15  -1.0 4.51e+03    -  1.00e+00 1.97e-03h  8
     147  0.0000000e+00 1.16e+04 1.15e+15  -1.0 2.48e+04    -  1.00e+00 1.97e-03h  9
     148  0.0000000e+00 1.15e+04 1.15e+15  -1.0 2.47e+04    -  1.00e+00 1.97e-03h  9
     149  0.0000000e+00 1.15e+04 1.15e+15  -1.0 2.46e+04    -  1.00e+00 1.97e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     150  0.0000000e+00 1.15e+04 1.15e+15  -1.0 2.46e+04    -  1.00e+00 1.97e-03h  9
     151  0.0000000e+00 1.15e+04 1.15e+15  -1.0 2.45e+04    -  1.00e+00 1.97e-03h  9
     152  0.0000000e+00 1.14e+04 1.15e+15  -1.0 2.44e+04    -  1.00e+00 1.97e-03h  9
     153  0.0000000e+00 1.14e+04 1.15e+15  -1.0 2.43e+04    -  1.00e+00 1.97e-03h  9
     154  0.0000000e+00 1.14e+04 1.15e+15  -1.0 2.42e+04    -  1.00e+00 1.96e-03h  9
     155  0.0000000e+00 1.14e+04 1.15e+15  -1.0 2.41e+04    -  1.00e+00 1.96e-03h  9
     156  0.0000000e+00 2.17e+05 1.30e+15  -1.0 2.40e+04    -  1.00e+00 5.03e-01w  1
     157  0.0000000e+00 9.99e+04 3.90e+15  -1.0 3.56e+04    -  6.23e-01 5.54e-01w  1
     158  0.0000000e+00 1.31e+05 7.87e+17  -1.0 1.57e+04    -  1.05e-01 9.31e-01w  1
     159  0.0000000e+00 1.14e+04 1.15e+15  -1.0 6.98e+05    -  1.00e+00 1.96e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     160  0.0000000e+00 1.13e+04 1.16e+15  -1.0 2.39e+04    -  1.00e+00 1.96e-03h  9
     161  0.0000000e+00 1.13e+04 1.16e+15  -1.0 2.37e+04    -  1.00e+00 1.96e-03h  9
     162  0.0000000e+00 1.13e+04 1.16e+15  -1.0 2.36e+04    -  1.00e+00 1.96e-03h  9
     163  0.0000000e+00 1.13e+04 1.16e+15  -1.0 2.35e+04    -  1.00e+00 1.96e-03h  9
     164  0.0000000e+00 1.12e+04 1.16e+15  -1.0 2.33e+04    -  1.00e+00 1.96e-03h  9
     165  0.0000000e+00 1.12e+04 1.16e+15  -1.0 2.32e+04    -  1.00e+00 1.96e-03h  9
     166  0.0000000e+00 1.12e+04 1.16e+15  -1.0 2.30e+04    -  1.00e+00 1.96e-03h  9
     167  0.0000000e+00 1.12e+04 1.16e+15  -1.0 2.28e+04    -  1.00e+00 1.96e-03h  9
     168  0.0000000e+00 1.12e+04 1.16e+15  -1.0 2.26e+04    -  1.00e+00 1.96e-03h  9
     169  0.0000000e+00 2.02e+05 1.34e+15  -1.0 2.23e+04    -  1.00e+00 5.03e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     170  0.0000000e+00 9.58e+04 4.16e+15  -1.0 3.79e+04    -  5.87e-01 5.55e-01w  1
     171  0.0000000e+00 8.97e+04 2.18e+18  -1.0 4.09e+03    -  1.04e-01 9.81e-01w  1
     172  0.0000000e+00 1.11e+04 1.16e+15  -1.0 4.36e+06    -  1.00e+00 1.96e-03h  8
     173  0.0000000e+00 1.11e+04 1.16e+15  -1.0 2.21e+04    -  1.00e+00 1.96e-03h  9
     174  0.0000000e+00 1.11e+04 1.16e+15  -1.0 2.17e+04    -  1.00e+00 1.96e-03h  9
     175  0.0000000e+00 1.11e+04 1.17e+15  -1.0 2.14e+04    -  1.00e+00 1.96e-03h  9
     176  0.0000000e+00 1.10e+04 1.17e+15  -1.0 2.09e+04    -  1.00e+00 1.96e-03h  9
     177  0.0000000e+00 1.10e+04 1.17e+15  -1.0 2.04e+04    -  1.00e+00 1.96e-03h  9
     178  0.0000000e+00 1.10e+04 1.17e+15  -1.0 1.98e+04    -  1.00e+00 1.96e-03h  9
     179  0.0000000e+00 1.10e+04 1.17e+15  -1.0 1.91e+04    -  1.00e+00 1.96e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     180  0.0000000e+00 1.10e+04 1.17e+15  -1.0 1.82e+04    -  1.00e+00 1.96e-03h  9
     181  0.0000000e+00 1.09e+04 1.17e+15  -1.0 1.72e+04    -  1.00e+00 3.92e-03h  8
     182  0.0000000e+00 1.59e+05 1.38e+15  -1.0 1.67e+04    -  1.00e+00 5.02e-01w  1
     183  0.0000000e+00 9.33e+04 4.59e+15  -1.0 4.74e+04    -  5.08e-01 5.55e-01w  1
     184  0.0000000e+00 9.75e+04 2.25e+18  -1.0 2.67e+03    -  1.02e-01 9.81e-01w  1
     185  0.0000000e+00 1.09e+04 1.17e+15  -1.0 3.80e+06    -  1.00e+00 3.92e-03h  7
     186  0.0000000e+00 1.08e+04 1.18e+15  -1.0 1.61e+04    -  1.00e+00 3.92e-03h  8
     187  0.0000000e+00 1.08e+04 1.18e+15  -1.0 1.54e+04    -  1.00e+00 3.91e-03h  8
     188  0.0000000e+00 2.03e+04 8.92e+15  -1.0 1.43e+04    -  7.56e-01 7.98e-01H  1
     189  0.0000000e+00 1.74e+04 6.11e+15  -1.0 1.09e+04    -  9.32e-01 1.24e-01h  3
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     190  0.0000000e+00 2.18e+03 3.24e+15  -1.0 7.66e+03    -  9.66e-01 5.03e-01h  1
     191  0.0000000e+00 2.18e+03 1.06e+16  -1.0 2.51e+04    -  2.14e-01 9.54e-04h 11
     192  0.0000000e+00 2.18e+03 1.06e+16  -1.0 2.57e+04    -  2.15e-01 9.54e-04h 11
     193  0.0000000e+00 2.17e+03 1.06e+16  -1.0 2.56e+04    -  2.16e-01 9.54e-04h 11
     194  0.0000000e+00 2.17e+03 1.06e+16  -1.0 2.56e+04    -  2.32e-01 9.54e-04h 11
     195  0.0000000e+00 2.17e+03 1.06e+16  -1.0 2.55e+04    -  1.00e+00 1.91e-03h 10
     196  0.0000000e+00 2.16e+03 1.06e+16  -1.0 2.53e+04    -  2.95e-01 1.91e-03h 10
     197  0.0000000e+00 2.16e+03 1.06e+16  -1.0 2.52e+04    -  1.00e+00 1.91e-03h 10
     198  0.0000000e+00 2.15e+03 1.05e+16  -1.0 2.50e+04    -  2.94e-01 1.91e-03h 10
     199  0.0000000e+00 2.15e+03 1.05e+16  -1.0 2.49e+04    -  1.00e+00 1.91e-03h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     200  0.0000000e+00 2.15e+03 1.05e+16  -1.0 2.47e+04    -  3.02e-01 1.91e-03h 10
     201  0.0000000e+00 1.99e+05 2.31e+17  -1.0 2.46e+04    -  1.00e+00 9.77e-01w  1
     202  0.0000000e+00 1.96e+05 2.27e+17  -1.0 8.73e+04    -  5.63e-02 1.85e-02w  1
     203  0.0000000e+00 5.11e+05 2.22e+17  -1.0 2.37e+05    -  3.03e-02 1.70e-01w  1
     204  0.0000000e+00 2.14e+03 1.05e+16  -1.0 2.44e+05    -  1.00e+00 1.91e-03h  9
     205  0.0000000e+00 2.14e+03 1.05e+16  -1.0 2.44e+04    -  3.10e-01 1.91e-03h 10
     206  0.0000000e+00 2.13e+03 1.04e+16  -1.0 2.43e+04    -  1.00e+00 1.91e-03h 10
     207  0.0000000e+00 2.13e+03 1.04e+16  -1.0 2.41e+04    -  3.18e-01 1.91e-03h 10
     208  0.0000000e+00 2.13e+03 1.04e+16  -1.0 2.40e+04    -  1.00e+00 1.91e-03h 10
     209  0.0000000e+00 2.12e+03 1.04e+16  -1.0 2.38e+04    -  3.26e-01 1.91e-03h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     210  0.0000000e+00 2.12e+03 1.04e+16  -1.0 2.37e+04    -  1.00e+00 1.91e-03h 10
     211  0.0000000e+00 2.11e+03 1.03e+16  -1.0 2.35e+04    -  3.35e-01 1.91e-03h 10
     212  0.0000000e+00 2.11e+03 1.03e+16  -1.0 2.34e+04    -  1.00e+00 1.91e-03h 10
     213  0.0000000e+00 2.11e+03 1.03e+16  -1.0 2.32e+04    -  3.43e-01 1.91e-03h 10
     214  0.0000000e+00 1.77e+05 2.35e+17  -1.0 2.31e+04    -  1.00e+00 9.78e-01w  1
     215  0.0000000e+00 1.75e+05 2.32e+17  -1.0 9.48e+04    -  4.85e-02 1.67e-02w  1
     216  0.0000000e+00 4.91e+05 2.27e+17  -1.0 2.41e+05    -  2.84e-02 1.65e-01w  1
     217  0.0000000e+00 2.10e+03 1.03e+16  -1.0 2.49e+05    -  1.00e+00 1.91e-03h  9
     218  0.0000000e+00 2.10e+03 1.03e+16  -1.0 2.29e+04    -  3.52e-01 1.91e-03h 10
     219  0.0000000e+00 2.09e+03 1.02e+16  -1.0 2.28e+04    -  1.00e+00 1.91e-03h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     220  0.0000000e+00 2.09e+03 1.02e+16  -1.0 2.26e+04    -  3.60e-01 1.91e-03h 10
     221  0.0000000e+00 2.09e+03 1.02e+16  -1.0 2.25e+04    -  1.00e+00 1.91e-03h 10
     222  0.0000000e+00 2.08e+03 1.02e+16  -1.0 2.23e+04    -  3.69e-01 1.91e-03h 10
     223  0.0000000e+00 2.08e+03 1.02e+16  -1.0 2.22e+04    -  1.00e+00 1.91e-03h 10
     224  0.0000000e+00 2.07e+03 1.01e+16  -1.0 2.20e+04    -  3.78e-01 1.91e-03h 10
     225  0.0000000e+00 2.07e+03 1.01e+16  -1.0 2.19e+04    -  1.00e+00 1.91e-03h 10
     226  0.0000000e+00 2.07e+03 1.01e+16  -1.0 2.17e+04    -  3.88e-01 1.91e-03h 10
     227  0.0000000e+00 1.58e+05 2.40e+17  -1.0 2.16e+04    -  1.00e+00 9.79e-01w  1
     228  0.0000000e+00 1.57e+05 2.38e+17  -1.0 1.03e+05    -  4.20e-02 1.50e-02w  1
     229  0.0000000e+00 4.73e+05 2.32e+17  -1.0 2.45e+05    -  2.68e-02 1.61e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     230  0.0000000e+00 2.06e+03 1.01e+16  -1.0 2.53e+05    -  1.00e+00 1.91e-03h  9
     231  0.0000000e+00 2.06e+03 1.01e+16  -1.0 2.15e+04    -  3.97e-01 1.91e-03h 10
     232  0.0000000e+00 2.05e+03 1.00e+16  -1.0 2.14e+04    -  1.00e+00 1.91e-03h 10
     233  0.0000000e+00 2.05e+03 1.00e+16  -1.0 2.12e+04    -  4.07e-01 1.91e-03h 10
     234  0.0000000e+00 2.05e+03 1.00e+16  -1.0 2.11e+04    -  1.00e+00 1.91e-03h 10
     235  0.0000000e+00 2.04e+03 9.99e+15  -1.0 2.09e+04    -  4.16e-01 1.91e-03h 10
     236  0.0000000e+00 2.04e+03 9.97e+15  -1.0 2.08e+04    -  1.00e+00 1.91e-03h 10
     237  0.0000000e+00 2.03e+03 9.95e+15  -1.0 2.07e+04    -  4.26e-01 1.91e-03h 10
     238  0.0000000e+00 2.03e+03 9.93e+15  -1.0 2.06e+04    -  1.00e+00 1.91e-03h 10
     239  0.0000000e+00 2.03e+03 9.91e+15  -1.0 2.04e+04    -  4.36e-01 1.91e-03h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     240  0.0000000e+00 1.41e+05 2.45e+17  -1.0 2.03e+04    -  1.00e+00 9.79e-01w  1
     241  0.0000000e+00 1.40e+05 2.43e+17  -1.0 1.12e+05    -  3.66e-02 1.35e-02w  1
     242  0.0000000e+00 4.56e+05 2.38e+17  -1.0 2.49e+05    -  2.53e-02 1.56e-01w  1
     243  0.0000000e+00 2.02e+03 9.90e+15  -1.0 2.59e+05    -  1.00e+00 1.91e-03h  9
     244  0.0000000e+00 2.02e+03 9.88e+15  -1.0 2.01e+04    -  4.47e-01 1.91e-03h 10
     245  0.0000000e+00 2.01e+03 9.86e+15  -1.0 2.00e+04    -  1.00e+00 1.91e-03h 10
     246  0.0000000e+00 2.01e+03 9.84e+15  -1.0 1.99e+04    -  4.57e-01 1.91e-03h 10
     247  0.0000000e+00 2.01e+03 9.82e+15  -1.0 1.98e+04    -  1.00e+00 1.91e-03h 10
     248  0.0000000e+00 2.00e+03 9.80e+15  -1.0 1.96e+04    -  4.68e-01 1.91e-03h 10
     249  0.0000000e+00 2.00e+03 9.78e+15  -1.0 1.95e+04    -  1.00e+00 1.91e-03h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     250  0.0000000e+00 2.00e+03 9.76e+15  -1.0 1.94e+04    -  4.79e-01 1.91e-03h 10
     251  0.0000000e+00 1.99e+03 9.74e+15  -1.0 1.93e+04    -  1.00e+00 1.91e-03h 10
     252  0.0000000e+00 1.99e+03 9.73e+15  -1.0 1.91e+04    -  4.90e-01 1.91e-03h 10
     253  0.0000000e+00 1.26e+05 2.50e+17  -1.0 1.91e+04    -  1.00e+00 9.80e-01w  1
     254  0.0000000e+00 1.25e+05 2.47e+17  -1.0 1.22e+05    -  3.20e-02 1.22e-02w  1
     255  0.0000000e+00 4.41e+05 2.42e+17  -1.0 2.53e+05    -  2.40e-02 1.52e-01w  1
     256  0.0000000e+00 1.98e+03 9.71e+15  -1.0 2.64e+05    -  1.00e+00 1.91e-03h  9
     257  0.0000000e+00 1.98e+03 9.69e+15  -1.0 1.89e+04    -  5.01e-01 1.91e-03h 10
     258  0.0000000e+00 1.98e+03 9.67e+15  -1.0 1.88e+04    -  1.00e+00 1.91e-03h 10
     259  0.0000000e+00 1.97e+03 9.65e+15  -1.0 1.87e+04    -  5.12e-01 1.91e-03h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     260  0.0000000e+00 1.97e+03 9.63e+15  -1.0 1.86e+04    -  1.00e+00 1.91e-03h 10
     261  0.0000000e+00 1.96e+03 9.61e+15  -1.0 1.84e+04    -  5.24e-01 1.91e-03h 10
     262  0.0000000e+00 3.56e+03 9.42e+15  -1.0 1.83e+04    -  3.05e-01 1.97e-02H  1
     263  0.0000000e+00 3.32e+03 8.79e+15  -1.0 1.44e+04    -  6.28e-03 6.14e-02f  5
     264  0.0000000e+00 3.21e+03 8.50e+15  -1.0 1.41e+04    -  1.00e+00 3.07e-02f  6
     265  0.0000000e+00 3.11e+03 8.23e+15  -1.0 1.39e+04    -  1.15e-01 3.07e-02h  6
     266  0.0000000e+00 3.01e+03 7.97e+15  -1.0 1.38e+04    -  1.00e+00 3.07e-02h  6
     267  0.0000000e+00 3.00e+03 7.94e+15  -1.0 1.36e+04    -  2.15e-01 3.84e-03h  9
     268  0.0000000e+00 2.99e+03 7.91e+15  -1.0 1.36e+04    -  1.00e+00 3.84e-03h  9
     269  0.0000000e+00 2.97e+03 7.88e+15  -1.0 1.35e+04    -  2.33e-01 3.84e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     270  0.0000000e+00 3.93e+03 7.79e+15  -1.0 1.35e+04    -  3.39e-01 1.12e-02H  1
     271  0.0000000e+00 3.67e+03 7.26e+15  -1.0 1.16e+04    -  2.22e-03 6.15e-02f  5
     272  0.0000000e+00 4.05e+03 7.23e+15  -1.0 1.15e+04    -  9.15e-01 4.30e-03H  1
     273  0.0000000e+00 3.78e+03 6.74e+15  -1.0 1.02e+04    -  7.00e-04 6.16e-02f  5
     274  0.0000000e+00 3.66e+03 6.53e+15  -1.0 1.02e+04    -  1.00e+00 3.08e-02f  6
     275  0.0000000e+00 3.60e+03 6.42e+15  -1.0 1.02e+04    -  1.15e-01 1.54e-02h  7
     276  0.0000000e+00 3.54e+03 6.32e+15  -1.0 1.02e+04    -  1.00e+00 1.54e-02h  7
     277  0.0000000e+00 3.49e+03 6.22e+15  -1.0 1.02e+04    -  1.62e-01 1.54e-02h  7
     278  0.0000000e+00 3.46e+03 6.18e+15  -1.0 1.02e+04    -  1.00e+00 7.70e-03h  8
     279  0.0000000e+00 4.29e+03 6.12e+15  -1.0 1.02e+04    -  3.70e-01 9.49e-03H  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     280  0.0000000e+00 4.00e+03 5.70e+15  -1.0 8.91e+03    -  1.45e-03 6.16e-02f  5
     281  0.0000000e+00 3.87e+03 5.51e+15  -1.0 8.98e+03    -  1.23e-01 3.08e-02f  6
     282  0.0000000e+00 4.42e+03 5.48e+15  -1.0 9.00e+03    -  1.00e+00 6.25e-03H  1
     283  0.0000000e+00 4.13e+03 5.11e+15  -1.0 7.76e+03    -  1.04e-03 6.16e-02f  5
     284  0.0000000e+00 4.00e+03 4.94e+15  -1.0 7.90e+03    -  1.00e+00 3.08e-02f  6
     285  0.0000000e+00 3.87e+03 4.78e+15  -1.0 7.96e+03    -  1.14e-01 3.08e-02h  6
     286  0.0000000e+00 3.81e+03 4.71e+15  -1.0 8.01e+03    -  1.00e+00 1.54e-02h  7
     287  0.0000000e+00 4.60e+03 4.66e+15  -1.0 8.04e+03    -  4.24e-01 9.37e-03H  1
     288  0.0000000e+00 4.30e+03 4.34e+15  -1.0 7.05e+03    -  1.39e-03 6.16e-02f  5
     289  0.0000000e+00 4.68e+03 4.32e+15  -1.0 7.19e+03    -  1.00e+00 4.46e-03H  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     290  0.0000000e+00 4.37e+03 3.96e+15  -1.0 6.34e+03    -  6.67e-04 6.16e-02f  5
     291  0.0000000e+00 4.74e+03 3.94e+15  -1.0 6.52e+03    -  5.82e-01 4.45e-03H  1
     292  0.0000000e+00 4.43e+03 3.59e+15  -1.0 5.93e+03    -  6.60e-04 6.16e-02f  5
     293  0.0000000e+00 4.80e+03 3.57e+15  -1.0 6.11e+03    -  1.00e+00 4.51e-03H  1
     294  0.0000000e+00 4.49e+03 3.14e+15  -1.0 5.37e+03    -  6.59e-04 6.16e-02f  5
     295  0.0000000e+00 4.86e+03 3.12e+15  -1.0 5.57e+03    -  6.18e-01 4.55e-03H  1
     296  0.0000000e+00 4.54e+03 2.65e+15  -1.0 5.05e+03    -  6.59e-04 6.16e-02f  5
     297  0.0000000e+00 4.91e+03 2.64e+15  -1.0 5.26e+03    -  1.00e+00 4.63e-03H  1
     298  0.0000000e+00 4.92e+03 5.89e+16  -1.0 4.85e+03    -  5.07e-01 4.52e-05H  1
     299  0.0000000e+00 4.86e+03 2.34e+16  -1.0 2.02e-05  16.0 6.03e-01 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     300  0.0000000e+00 1.36e+05 1.93e+16  -1.0 1.02e-04  15.5 1.74e-01 5.79e-01h  1
     301  0.0000000e+00 1.39e+05 3.86e+15  -1.0 1.24e-05  16.9 9.91e-01 4.77e-02h  1
     302  0.0000000e+00 1.39e+05 3.85e+15  -1.0 3.50e-06  17.3 1.69e-03 1.69e-03s 13
     303r 0.0000000e+00 1.39e+05 1.00e+03  -1.0 0.00e+00  16.8 0.00e+00 0.00e+00R  1
     304r 0.0000000e+00 1.31e+05 3.68e+03  -1.0 4.09e+02    -  1.99e-02 6.01e-03f  1
     305  0.0000000e+00 1.32e+05 5.01e+04  -1.0 3.97e+05    -  1.67e-01 6.31e-02f  1
     306  0.0000000e+00 5.61e+04 6.22e+05  -1.0 3.32e+04    -  6.78e-02 3.13e-01H  1
     307  0.0000000e+00 5.48e+04 3.65e+05  -1.0 5.60e+04    -  6.45e-01 2.42e-02h  2
     308  0.0000000e+00 4.52e+04 6.25e+05  -1.0 2.94e+04    -  9.90e-01 1.64e-01h  1
     309  0.0000000e+00 4.52e+04 1.20e+06  -1.0 4.90e+04    -  9.77e-01 4.22e-04h 12
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     310  0.0000000e+00 4.52e+04 1.87e+06  -1.0 4.89e+04    -  9.96e-01 2.10e-04h 13
     311  0.0000000e+00 4.52e+04 2.65e+06  -1.0 4.90e+04    -  1.00e+00 2.11e-04h 13
     312  0.0000000e+00 4.51e+04 3.54e+06  -1.0 4.90e+04    -  1.00e+00 2.11e-04h 13
     313  0.0000000e+00 4.51e+04 4.55e+06  -1.0 4.91e+04    -  1.00e+00 2.11e-04h 13
     314  0.0000000e+00 4.51e+04 5.70e+06  -1.0 4.92e+04    -  1.00e+00 2.11e-04h 13
     315  0.0000000e+00 4.51e+04 7.00e+06  -1.0 4.94e+04    -  1.00e+00 2.12e-04h 13
     316  0.0000000e+00 4.51e+04 8.43e+06  -1.0 4.97e+04    -  1.00e+00 2.13e-04h 13
     317  0.0000000e+00 4.48e+04 9.89e+06  -1.0 5.03e+04    -  1.00e+00 6.89e-03h  8
     318  0.0000000e+00 4.47e+04 1.12e+07  -1.0 5.23e+04    -  1.00e+00 1.78e-03h 10
     319  0.0000000e+00 1.41e+06 3.60e+06  -1.0 5.68e+04    -  1.00e+00 9.17e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     320  0.0000000e+00 9.89e+05 1.07e+07  -1.0 5.23e+04    -  8.52e-01 4.87e-01w  1
     321  0.0000000e+00 2.90e+06 6.67e+08  -1.0 2.15e+05    -  1.87e-01 3.91e-01w  1
     322  0.0000000e+00 4.47e+04 1.16e+07  -1.0 8.00e+04    -  1.00e+00 4.48e-04h 11
     323  0.0000000e+00 4.47e+04 1.27e+07  -1.0 5.40e+04    -  9.80e-01 4.59e-04h 12
     324  0.0000000e+00 1.60e+05 3.31e+06  -1.0 1.91e+05    -  3.28e-01 1.00e-01f  2
     325  0.0000000e+00 1.58e+05 3.75e+06  -1.0 7.11e+04    -  1.00e+00 3.70e-02h  5
     326  0.0000000e+00 1.57e+05 4.36e+06  -1.0 7.53e+04    -  1.00e+00 4.05e-03h  8
     327  0.0000000e+00 1.57e+05 5.02e+06  -1.0 7.73e+04    -  1.00e+00 6.07e-05h 14
     328  0.0000000e+00 1.57e+05 5.81e+06  -1.0 7.56e+04    -  1.00e+00 1.57e-05h 16
     329r 0.0000000e+00 1.57e+05 1.00e+03   1.6 0.00e+00    -  0.00e+00 2.53e-07R 22
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     330r 0.0000000e+00 3.94e+04 6.98e+03   1.6 5.27e+04    -  3.58e-02 9.98e-04f  1
     331  0.0000000e+00 2.65e+04 6.48e+03  -1.0 4.08e+04    -  7.02e-01 4.95e-01h  2
     332  0.0000000e+00 1.90e+04 6.56e+03  -1.0 1.14e+04    -  9.86e-01 1.99e-01h  2
     333  0.0000000e+00 1.79e+04 1.68e+04  -1.0 9.22e+03    -  9.90e-01 5.21e-02h  4
     334  0.0000000e+00 1.78e+04 3.16e+04  -1.0 3.08e+04    -  1.00e+00 4.65e-03h  8
     335  0.0000000e+00 1.78e+04 4.50e+04  -1.0 4.43e+04    -  1.00e+00 4.14e-03h  8
     336  0.0000000e+00 1.77e+04 5.94e+04  -1.0 4.74e+04    -  1.00e+00 3.93e-03h  8
     337  0.0000000e+00 1.76e+04 7.79e+04  -1.0 4.70e+04    -  1.00e+00 3.94e-03h  8
     338  0.0000000e+00 1.76e+04 1.28e+05  -1.0 4.53e+04    -  1.00e+00 4.05e-03h  8
     339  0.0000000e+00 1.75e+04 2.37e+05  -1.0 4.27e+04    -  1.00e+00 4.22e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     340  0.0000000e+00 1.74e+04 4.18e+05  -1.0 3.99e+04    -  1.00e+00 4.42e-03h  8
     341  0.0000000e+00 6.07e+05 1.20e+06  -1.0 3.77e+04    -  1.00e+00 5.88e-01w  1
     342  0.0000000e+00 5.63e+05 2.96e+05  -1.0 3.64e+04    -  1.00e+00 4.74e-01w  1
     343  0.0000000e+00 3.22e+05 1.51e+06  -1.0 1.05e+04    -  2.68e-01 4.97e-01w  1
     344  0.0000000e+00 1.73e+04 7.06e+05  -1.0 7.54e+03    -  1.00e+00 4.60e-03h  7
     345  0.0000000e+00 1.72e+04 1.16e+06  -1.0 3.62e+04    -  1.00e+00 4.72e-03h  8
     346  0.0000000e+00 1.72e+04 1.87e+06  -1.0 3.52e+04    -  1.00e+00 4.80e-03h  8
     347  0.0000000e+00 1.71e+04 2.98e+06  -1.0 3.46e+04    -  1.00e+00 4.85e-03h  8
     348  0.0000000e+00 1.70e+04 4.72e+06  -1.0 3.43e+04    -  1.00e+00 4.87e-03h  8
     349  0.0000000e+00 1.69e+04 7.48e+06  -1.0 3.40e+04    -  1.00e+00 4.86e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     350  0.0000000e+00 1.68e+04 1.18e+07  -1.0 3.39e+04    -  1.00e+00 4.84e-03h  8
     351  0.0000000e+00 1.67e+04 1.87e+07  -1.0 3.38e+04    -  1.00e+00 4.84e-03h  8
     352  0.0000000e+00 1.67e+04 2.96e+07  -1.0 3.37e+04    -  1.00e+00 4.83e-03h  8
     353  0.0000000e+00 1.66e+04 4.68e+07  -1.0 3.37e+04    -  1.00e+00 4.82e-03h  8
     354  0.0000000e+00 5.87e+05 1.41e+08  -1.0 3.36e+04    -  1.00e+00 6.17e-01w  1
     355  0.0000000e+00 5.29e+05 3.31e+07  -1.0 4.27e+04    -  1.00e+00 3.02e-01w  1
     356  0.0000000e+00 4.41e+05 8.78e+07  -1.0 2.52e+04    -  1.00e+00 2.11e-01w  1
     357  0.0000000e+00 1.65e+04 7.42e+07  -1.0 1.62e+04    -  1.00e+00 4.82e-03h  7
     358  0.0000000e+00 1.64e+04 1.18e+08  -1.0 3.36e+04    -  1.00e+00 4.81e-03h  8
     359  0.0000000e+00 1.63e+04 1.87e+08  -1.0 3.35e+04    -  1.00e+00 4.81e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     360  0.0000000e+00 1.63e+04 2.97e+08  -1.0 3.35e+04    -  1.00e+00 4.80e-03h  8
     361  0.0000000e+00 1.61e+04 4.71e+08  -1.0 3.35e+04    -  1.00e+00 9.60e-03h  7
     362  0.0000000e+00 1.60e+04 7.46e+08  -1.0 3.31e+04    -  1.00e+00 9.55e-03h  7
     363  0.0000000e+00 1.58e+04 1.19e+09  -1.0 3.31e+04    -  1.00e+00 9.54e-03h  7
     364  0.0000000e+00 1.56e+04 1.89e+09  -1.0 3.30e+04    -  1.00e+00 9.52e-03h  7
     365  0.0000000e+00 1.55e+04 3.02e+09  -1.0 3.29e+04    -  1.00e+00 9.51e-03h  7
     366  0.0000000e+00 1.54e+04 4.85e+09  -1.0 3.29e+04    -  1.00e+00 9.50e-03h  7
     367  0.0000000e+00 5.26e+05 1.35e+10  -1.0 3.28e+04    -  1.00e+00 6.07e-01w  1
     368  0.0000000e+00 4.79e+05 3.97e+09  -1.0 4.04e+04    -  1.00e+00 3.20e-01w  1
     369  0.0000000e+00 4.26e+05 1.17e+10  -1.0 2.35e+04    -  1.00e+00 1.26e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     370  0.0000000e+00 1.52e+04 7.79e+09  -1.0 1.43e+04    -  1.00e+00 9.48e-03h  6
     371  0.0000000e+00 1.51e+04 1.25e+10  -1.0 3.27e+04    -  1.00e+00 9.47e-03h  7
     372  0.0000000e+00 1.49e+04 2.02e+10  -1.0 3.26e+04    -  1.00e+00 9.45e-03h  7
     373  0.0000000e+00 1.48e+04 3.27e+10  -1.0 3.25e+04    -  1.00e+00 9.44e-03h  7
     374  0.0000000e+00 1.46e+04 5.30e+10  -1.0 3.24e+04    -  1.00e+00 9.42e-03h  7
     375  0.0000000e+00 1.45e+04 8.61e+10  -1.0 3.23e+04    -  1.00e+00 9.41e-03h  7
     376  0.0000000e+00 1.44e+04 1.40e+11  -1.0 3.22e+04    -  1.00e+00 9.39e-03h  7
     377  0.0000000e+00 1.42e+04 2.28e+11  -1.0 3.21e+04    -  1.00e+00 9.38e-03h  7
     378  0.0000000e+00 1.41e+04 3.72e+11  -1.0 3.20e+04    -  1.00e+00 9.36e-03h  7
     379  0.0000000e+00 1.40e+04 6.08e+11  -1.0 3.20e+04    -  1.00e+00 9.35e-03h  7
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     380  0.0000000e+00 4.63e+05 1.57e+12  -1.0 3.19e+04    -  1.00e+00 5.98e-01w  1
     381  0.0000000e+00 4.22e+05 5.65e+11  -1.0 3.68e+04    -  1.00e+00 3.36e-01w  1
     382  0.0000000e+00 4.14e+05 4.53e+12  -1.0 2.15e+04    -  6.11e-01 2.01e-02w  1
     383  0.0000000e+00 1.38e+04 9.96e+11  -1.0 1.15e+04    -  1.00e+00 9.34e-03h  6
     384  0.0000000e+00 1.37e+04 1.63e+12  -1.0 3.18e+04    -  1.00e+00 9.32e-03h  7
     385  0.0000000e+00 1.36e+04 2.68e+12  -1.0 3.17e+04    -  1.00e+00 9.31e-03h  7
     386  0.0000000e+00 1.34e+04 4.41e+12  -1.0 3.16e+04    -  1.00e+00 9.30e-03h  7
     387  0.0000000e+00 1.33e+04 7.24e+12  -1.0 3.17e+04    -  1.00e+00 9.32e-03h  7
     388  0.0000000e+00 1.32e+04 1.19e+13  -1.0 3.15e+04    -  1.00e+00 9.29e-03h  7
     389  0.0000000e+00 1.31e+04 1.96e+13  -1.0 3.17e+04    -  1.00e+00 9.31e-03h  7
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     390  0.0000000e+00 1.29e+04 3.09e+13  -1.0 3.47e+04    -  1.00e+00 9.68e-03h  7
     391  0.0000000e+00 1.28e+04 3.61e+13  -1.0 3.30e+04    -  1.00e+00 9.46e-03h  7
     392  0.0000000e+00 1.27e+04 4.63e+13  -1.0 4.45e+04    -  1.00e+00 8.11e-03h  7
     393  0.0000000e+00 4.27e+05 4.51e+13  -1.0 5.62e+04    -  1.00e+00 4.48e-01w  1
     394  0.0000000e+00 4.02e+05 4.18e+13  -1.0 2.38e+04    -  1.00e+00 4.86e-01w  1
     395  0.0000000e+00 3.59e+05 1.30e+14  -1.0 1.33e+04    -  1.00e+00 1.16e-01w  1
     396  0.0000000e+00 1.27e+04 5.49e+13  -1.0 4.61e+03    -  1.00e+00 3.50e-03h  7
     397  0.0000000e+00 1.26e+04 4.26e+13  -1.0 7.54e+04    -  1.00e+00 9.17e-03h  7
     398  0.0000000e+00 1.25e+04 5.58e+13  -1.0 4.78e+04    -  1.00e+00 3.39e-03h  8
     399  0.0000000e+00 1.24e+04 5.58e+13  -1.0 6.59e+04    -  1.00e+00 7.17e-03h  7
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     400  0.0000000e+00 1.20e+04 1.96e+12  -1.0 6.71e-06  16.3 9.65e-01 1.00e+00h  1
     401  0.0000000e+00 4.20e+04 3.57e+12  -1.0 2.23e-05  16.7 1.00e+00 1.00e+00h  1
     402  0.0000000e+00 4.29e+04 4.84e+12  -1.0 3.93e-05  17.2 1.00e+00 2.95e-01h  1
     403  0.0000000e+00 5.93e+04 7.51e+12  -1.0 1.88e-05  17.6 1.00e+00 5.75e-01h  1
     404r 0.0000000e+00 5.93e+04 1.00e+03   0.4 0.00e+00  18.0 0.00e+00 4.13e-07R 15
     405r 0.0000000e+00 5.15e+04 5.46e+04   0.4 3.26e+03    -  4.96e-02 9.91e-04f  1
     406  0.0000000e+00 2.60e+04 2.24e+04  -1.0 2.71e+04    -  6.61e-01 4.47e-01h  1
     407  0.0000000e+00 1.19e+04 9.55e+05  -1.0 3.63e+04    -  6.24e-01 4.97e-01h  1
     408  0.0000000e+00 1.03e+04 2.07e+08  -1.0 3.15e+04    -  9.85e-01 1.33e-01h  3
     409  0.0000000e+00 1.03e+04 7.38e+08  -1.0 3.79e+04    -  9.90e-01 4.20e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     410  0.0000000e+00 1.03e+04 1.73e+09  -1.0 3.92e+04    -  1.00e+00 2.10e-03h  9
     411  0.0000000e+00 1.03e+04 3.55e+09  -1.0 3.91e+04    -  1.00e+00 1.05e-03h 10
     412  0.0000000e+00 1.03e+04 6.91e+09  -1.0 3.89e+04    -  1.00e+00 1.05e-03h 10
     413  0.0000000e+00 1.02e+04 1.31e+10  -1.0 3.86e+04    -  1.00e+00 1.05e-03h 10
     414  0.0000000e+00 1.02e+04 2.45e+10  -1.0 3.79e+04    -  1.00e+00 1.05e-03h 10
     415  0.0000000e+00 1.02e+04 4.54e+10  -1.0 3.66e+04    -  1.00e+00 1.05e-03h 10
     416  0.0000000e+00 1.02e+04 8.39e+10  -1.0 3.45e+04    -  1.00e+00 2.10e-03h  9
     417  0.0000000e+00 1.02e+04 1.55e+11  -1.0 3.07e+04    -  1.00e+00 2.10e-03h  9
     418  0.0000000e+00 1.31e+05 4.39e+10  -1.0 2.41e+04    -  1.00e+00 5.38e-01w  1
     419  0.0000000e+00 1.42e+05 9.72e+11  -1.0 4.64e+04    -  1.64e-01 4.36e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     420  0.0000000e+00 1.59e+05 9.05e+11  -1.0 1.32e+05    -  5.57e-01 6.78e-02w  1
     421  0.0000000e+00 1.02e+04 2.85e+11  -1.0 6.07e+04    -  1.00e+00 2.10e-03h  8
     422  0.0000000e+00 1.01e+04 5.22e+11  -1.0 1.27e+04    -  1.00e+00 4.21e-03h  8
     423  0.0000000e+00 9.76e+03 9.29e+11  -1.0 5.70e+03    -  1.00e+00 3.37e-02h  5
     424  0.0000000e+00 9.41e+03 1.67e+12  -1.0 2.48e+04    -  1.00e+00 3.38e-02h  5
     425  0.0000000e+00 9.37e+03 3.11e+12  -1.0 4.18e+04    -  1.00e+00 4.24e-03h  8
     426  0.0000000e+00 9.36e+03 5.74e+12  -1.0 7.03e+04    -  1.00e+00 8.26e-04h 10
     427  0.0000000e+00 9.36e+03 1.05e+13  -1.0 8.23e+04    -  1.00e+00 7.36e-04h 10
     428  0.0000000e+00 9.35e+03 1.92e+13  -1.0 7.84e+04    -  1.00e+00 7.62e-04h 10
     429  0.0000000e+00 9.34e+03 3.51e+13  -1.0 6.98e+04    -  1.00e+00 8.29e-04h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     430  0.0000000e+00 9.33e+03 6.40e+13  -1.0 6.16e+04    -  1.00e+00 9.05e-04h 10
     431  0.0000000e+00 1.78e+05 5.28e+12  -1.0 5.62e+04    -  1.00e+00 4.93e-01w  1
     432  0.0000000e+00 1.35e+05 1.47e+14  -1.0 2.30e+04    -  7.85e-02 8.13e-01w  1
     433  0.0000000e+00 2.04e+05 3.84e+14  -1.0 4.97e+04    -  7.86e-02 3.19e-01w  1
     434  0.0000000e+00 9.32e+03 1.16e+14  -1.0 6.05e-06  17.6 1.00e+00 1.92e-03h  8
     435  0.0000000e+00 9.30e+03 2.12e+14  -1.0 5.31e+04    -  1.00e+00 2.00e-03h  9
     436  0.0000000e+00 9.28e+03 3.85e+14  -1.0 5.15e+04    -  1.00e+00 2.04e-03h  9
     437  0.0000000e+00 9.26e+03 7.00e+14  -1.0 5.06e+04    -  1.00e+00 2.06e-03h  9
     438  0.0000000e+00 9.24e+03 1.27e+15  -1.0 5.00e+04    -  1.00e+00 2.07e-03h  9
     439  0.0000000e+00 9.22e+03 2.32e+15  -1.0 4.96e+04    -  1.00e+00 2.08e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     440  0.0000000e+00 9.20e+03 4.21e+15  -1.0 4.93e+04    -  1.00e+00 2.09e-03h  9
     441  0.0000000e+00 9.18e+03 7.66e+15  -1.0 4.90e+04    -  1.00e+00 2.10e-03h  9
     442  0.0000000e+00 9.16e+03 1.39e+16  -1.0 4.87e+04    -  1.00e+00 2.10e-03h  9
     443  0.0000000e+00 9.14e+03 2.53e+16  -1.0 4.85e+04    -  1.00e+00 2.10e-03h  9
     444  0.0000000e+00 1.39e+05 6.88e+15  -1.0 4.83e+04    -  1.00e+00 5.40e-01w  1
     445  0.0000000e+00 1.41e+05 4.47e+16  -1.0 9.28e+04    -  3.76e-01 5.67e-02w  1
     446  0.0000000e+00 1.60e+05 2.96e+16  -1.0 5.03e+04    -  4.49e-01 2.43e-01w  1
     447  0.0000000e+00 9.12e+03 4.60e+16  -1.0 3.78e+04    -  1.00e+00 2.11e-03h  8
     448  0.0000000e+00 9.10e+03 8.35e+16  -1.0 4.82e+04    -  1.00e+00 2.11e-03h  9
     449  0.0000000e+00 9.08e+03 1.52e+17  -1.0 4.80e+04    -  1.00e+00 2.11e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     450  0.0000000e+00 9.06e+03 2.75e+17  -1.0 4.79e+04    -  1.00e+00 2.12e-03h  9
     451  0.0000000e+00 9.05e+03 5.00e+17  -1.0 4.77e+04    -  1.00e+00 2.12e-03h  9
     452  0.0000000e+00 9.03e+03 9.07e+17  -1.0 4.76e+04    -  1.00e+00 2.12e-03h  9
     453  0.0000000e+00 9.01e+03 1.65e+18  -1.0 4.75e+04    -  1.00e+00 2.12e-03h  9
     454  0.0000000e+00 8.99e+03 2.99e+18  -1.0 4.74e+04    -  1.00e+00 2.12e-03h  9
     455  0.0000000e+00 8.97e+03 4.24e+18  -1.0 4.73e+04    -  1.00e+00 2.13e-03h  9
     456  0.0000000e+00 8.89e+03 4.24e+18  -1.0 3.05e+04    -  1.00e+00 8.52e-03h  7
     457  0.0000000e+00 5.11e+03 1.39e+18  -1.0 1.22e+04    -  1.00e+00 5.46e-01w  1
     458  0.0000000e+00 8.90e+03 7.09e+18  -1.0 8.20e-06  18.9 9.90e-01 1.10e-01w  1
     459  0.0000000e+00 8.99e+03 1.16e+19  -1.0 1.94e-05  18.4 9.92e-01 4.55e-02w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     460  0.0000000e+00 5.34e+03 4.87e+18  -1.0 4.90e-05  17.9 1.00e+00 2.73e-01h  1
     461  0.0000000e+00 5.33e+03 5.60e+18  -1.0 1.76e+05    -  2.11e-01 1.94e-03h  9
     462  0.0000000e+00 9.33e+03 4.81e+18  -1.0 1.36e-05  17.4 1.39e-01 1.39e-01s 20
     463  0.0000000e+00 9.44e+03 4.44e+18  -1.0 7.19e-06  17.0 7.70e-02 7.70e-02s 20
     464  0.0000000e+00 9.44e+03 3.93e+18  -1.0 6.96e-06  17.4 1.13e-01 1.13e-01s 20
     465  0.0000000e+00 1.11e+04 3.18e+18  -1.0 8.25e-06  16.9 1.88e-01 1.88e-01s 20
     466  0.0000000e+00 1.55e+04 1.63e+18  -1.0 8.43e-06  17.3 4.54e-01 4.54e-01s 20
     467  0.0000000e+00 6.58e+04 3.60e+17  -1.0 1.36e-05  16.9 4.66e-01 4.66e-01s 20
     468  0.0000000e+00 6.67e+04 3.18e+17  -1.0 8.13e-06  18.2 7.88e-02 7.88e-02s 20
     469  0.0000000e+00 6.67e+04 3.16e+17  -1.0 1.21e-05  17.7 7.44e-03 7.44e-03s 20
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     470  0.0000000e+00 6.67e+04 3.14e+17  -1.0 1.44e-05  17.2 6.24e-03 6.24e-03s 20
     471r 0.0000000e+00 6.67e+04 1.00e+03  -0.6 0.00e+00  17.7 0.00e+00 0.00e+00R  1
     472r 0.0000000e+00 6.33e+04 1.28e+05  -0.6 7.20e+03    -  1.43e-02 9.96e-04f  1
     473  0.0000000e+00 3.72e+04 4.09e+05  -1.0 1.12e+04    -  9.89e-01 3.69e-01h  1
     474  0.0000000e+00 3.66e+04 2.27e+07  -1.0 2.43e+03    -  6.06e-01 1.68e-02h  6
     475  0.0000000e+00 3.65e+04 8.70e+07  -1.0 2.80e+04    -  9.90e-01 9.55e-04h 10
     476  0.0000000e+00 3.65e+04 1.96e+08  -1.0 2.03e+04    -  8.67e-01 9.81e-04h 10
     477  0.0000000e+00 3.65e+04 4.18e+08  -1.0 1.24e+04    -  1.00e+00 1.01e-03h 10
     478  0.0000000e+00 3.58e+04 7.71e+08  -1.0 8.26e+03    -  1.00e+00 1.74e-02h  6
     479  0.0000000e+00 3.56e+04 1.32e+09  -1.0 2.58e+04    -  1.00e+00 4.67e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     480  0.0000000e+00 3.56e+04 2.08e+09  -1.0 4.48e+04    -  1.00e+00 6.33e-04h 11
     481  0.0000000e+00 3.56e+04 3.02e+09  -1.0 5.92e+04    -  1.00e+00 3.38e-04h 12
     482  0.0000000e+00 3.56e+04 4.20e+09  -1.0 6.77e+04    -  1.00e+00 3.20e-04h 12
     483  0.0000000e+00 3.56e+04 5.69e+09  -1.0 7.26e+04    -  1.00e+00 2.90e-04h 12
     484  0.0000000e+00 4.44e+05 3.52e+09  -1.0 7.48e+04    -  1.00e+00 5.70e-01w  1
     485  0.0000000e+00 2.32e+05 6.20e+09  -1.0 1.51e+04    -  4.48e-02 5.39e-01w  1
     486  0.0000000e+00 3.41e+05 1.25e+10  -1.0 2.71e+04    -  5.97e-01 9.82e-01w  1
     487  0.0000000e+00 3.56e+04 7.60e+09  -1.0 2.45e-05  17.2 1.00e+00 1.39e-04h 12
     488  0.0000000e+00 3.56e+04 1.01e+10  -1.0 7.51e+04    -  1.00e+00 1.39e-04h 13
     489  0.0000000e+00 3.56e+04 1.35e+10  -1.0 7.40e+04    -  1.00e+00 1.41e-04h 13
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     490  0.0000000e+00 3.56e+04 1.83e+10  -1.0 7.21e+04    -  1.00e+00 2.92e-04h 12
     491  0.0000000e+00 3.56e+04 2.48e+10  -1.0 7.01e+04    -  1.00e+00 3.04e-04h 12
     492  0.0000000e+00 3.55e+04 3.41e+10  -1.0 6.84e+04    -  1.00e+00 3.14e-04h 12
     493  0.0000000e+00 3.55e+04 4.71e+10  -1.0 6.71e+04    -  1.00e+00 3.23e-04h 12
     494  0.0000000e+00 3.55e+04 6.53e+10  -1.0 6.61e+04    -  1.00e+00 3.29e-04h 12
     495  0.0000000e+00 3.55e+04 9.08e+10  -1.0 6.55e+04    -  1.00e+00 3.34e-04h 12
     496  0.0000000e+00 3.55e+04 1.27e+11  -1.0 6.50e+04    -  1.00e+00 3.37e-04h 12
     497  0.0000000e+00 4.87e+05 4.73e+10  -1.0 6.47e+04    -  1.00e+00 6.96e-01w  1
     498  0.0000000e+00 3.45e+05 1.64e+11  -1.0 2.98e+04    -  2.01e-01 4.86e-01w  1
     499  0.0000000e+00 3.44e+05 1.03e+11  -1.0 1.24e+05    -  1.00e+00 4.24e-02w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     500  0.0000000e+00 3.55e+04 1.77e+11  -1.0 1.57e+04    -  1.00e+00 3.40e-04h 11
     501  0.0000000e+00 3.55e+04 2.47e+11  -1.0 6.44e+04    -  1.00e+00 3.41e-04h 12
     502  0.0000000e+00 3.55e+04 3.45e+11  -1.0 6.43e+04    -  1.00e+00 3.42e-04h 12
     503  0.0000000e+00 3.54e+04 4.83e+11  -1.0 6.42e+04    -  1.00e+00 3.43e-04h 12
     504  0.0000000e+00 3.54e+04 6.77e+11  -1.0 6.41e+04    -  1.00e+00 3.44e-04h 12
     505  0.0000000e+00 3.54e+04 9.48e+11  -1.0 6.40e+04    -  1.00e+00 3.44e-04h 12
     506  0.0000000e+00 3.54e+04 1.33e+12  -1.0 6.40e+04    -  1.00e+00 3.44e-04h 12
     507  0.0000000e+00 3.54e+04 1.86e+12  -1.0 6.39e+04    -  1.00e+00 3.44e-04h 12
     508  0.0000000e+00 3.54e+04 2.61e+12  -1.0 6.39e+04    -  1.00e+00 3.44e-04h 12
     509  0.0000000e+00 3.54e+04 3.66e+12  -1.0 6.39e+04    -  1.00e+00 3.44e-04h 12
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     510  0.0000000e+00 4.83e+05 1.25e+12  -1.0 6.38e+04    -  1.00e+00 7.06e-01w  1
     511  0.0000000e+00 3.45e+05 5.05e+12  -1.0 3.02e+04    -  2.35e-01 5.02e-01w  1
     512  0.0000000e+00 3.45e+05 1.45e+12  -1.0 6.92e-05  16.7 2.63e-01 5.63e-03w  1
     513  0.0000000e+00 3.54e+04 5.13e+12  -1.0 2.09e-04  16.2 1.00e+00 3.44e-04h 11
     514  0.0000000e+00 3.53e+04 7.19e+12  -1.0 6.38e+04    -  1.00e+00 3.45e-04h 12
     515  0.0000000e+00 3.53e+04 1.01e+13  -1.0 6.38e+04    -  1.00e+00 3.45e-04h 12
     516  0.0000000e+00 3.53e+04 1.42e+13  -1.0 6.38e+04    -  1.00e+00 3.45e-04h 12
     517  0.0000000e+00 3.53e+04 1.99e+13  -1.0 6.37e+04    -  1.00e+00 3.44e-04h 12
     518  0.0000000e+00 3.53e+04 2.79e+13  -1.0 6.37e+04    -  1.00e+00 3.44e-04h 12
     519  0.0000000e+00 3.53e+04 3.91e+13  -1.0 6.37e+04    -  1.00e+00 3.44e-04h 12
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     520  0.0000000e+00 3.53e+04 5.49e+13  -1.0 6.37e+04    -  1.00e+00 3.44e-04h 12
     521  0.0000000e+00 3.53e+04 7.71e+13  -1.0 6.37e+04    -  1.00e+00 3.44e-04h 12
     522  0.0000000e+00 3.53e+04 1.08e+14  -1.0 6.37e+04    -  1.00e+00 3.44e-04h 12
     523  0.0000000e+00 4.80e+05 3.67e+13  -1.0 6.36e+04    -  1.00e+00 7.05e-01w  1
     524  0.0000000e+00 3.43e+05 1.50e+14  -1.0 3.01e+04    -  2.40e-01 5.05e-01w  1
     525  0.0000000e+00 3.43e+05 2.90e+13  -1.0 1.21e-06  18.5 3.87e-01 3.69e-03w  1
     526  0.0000000e+00 3.52e+04 1.52e+14  -1.0 3.64e-06  18.0 1.00e+00 3.44e-04h 11
     527  0.0000000e+00 3.52e+04 2.14e+14  -1.0 6.36e+04    -  1.00e+00 3.44e-04h 12
     528  0.0000000e+00 3.52e+04 3.00e+14  -1.0 6.36e+04    -  1.00e+00 3.44e-04h 12
     529  0.0000000e+00 3.52e+04 4.22e+14  -1.0 6.36e+04    -  1.00e+00 3.44e-04h 12
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     530  0.0000000e+00 3.52e+04 5.92e+14  -1.0 6.36e+04    -  1.00e+00 3.44e-04h 12
     531  0.0000000e+00 3.52e+04 8.32e+14  -1.0 6.36e+04    -  1.00e+00 3.44e-04h 12
     532  0.0000000e+00 3.52e+04 1.17e+15  -1.0 6.35e+04    -  1.00e+00 3.44e-04h 12
     533  0.0000000e+00 3.52e+04 1.64e+15  -1.0 6.35e+04    -  1.00e+00 3.44e-04h 12
     534  0.0000000e+00 3.51e+04 2.31e+15  -1.0 6.35e+04    -  1.00e+00 3.44e-04h 12
     535  0.0000000e+00 3.51e+04 3.25e+15  -1.0 6.35e+04    -  1.00e+00 3.44e-04h 12
     536  0.0000000e+00 4.77e+05 1.10e+15  -1.0 6.35e+04    -  1.00e+00 7.04e-01w  1
     537  0.0000000e+00 3.41e+05 4.51e+15  -1.0 3.01e+04    -  2.42e-01 5.05e-01w  1
    WARNING: Problem in step computation; switching to emergency mode.
     538  0.0000000e+00 3.51e+04 4.57e+15  -1.0 3.01e+04    -  1.00e+00 3.44e-04h 12
     539  0.0000000e+00 3.51e+04 6.43e+15  -1.0 6.35e+04    -  1.00e+00 3.44e-04h 12
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     540  0.0000000e+00 3.51e+04 9.04e+15  -1.0 6.35e+04    -  1.00e+00 3.44e-04h 12
     541  0.0000000e+00 3.51e+04 1.27e+16  -1.0 6.34e+04    -  1.00e+00 3.44e-04h 12
     542  0.0000000e+00 3.51e+04 1.79e+16  -1.0 6.34e+04    -  1.00e+00 3.44e-04h 12
     543  0.0000000e+00 3.51e+04 2.52e+16  -1.0 6.34e+04    -  1.00e+00 3.44e-04h 12
     544  0.0000000e+00 3.50e+04 3.54e+16  -1.0 6.34e+04    -  1.00e+00 3.43e-04h 12
     545  0.0000000e+00 3.50e+04 4.99e+16  -1.0 6.34e+04    -  1.00e+00 3.43e-04h 12
     546  0.0000000e+00 3.50e+04 7.02e+16  -1.0 6.34e+04    -  1.00e+00 3.43e-04h 12
     547  0.0000000e+00 3.50e+04 9.88e+16  -1.0 6.33e+04    -  1.00e+00 3.43e-04h 12
     548  0.0000000e+00 4.75e+05 3.32e+16  -1.0 6.34e+04    -  1.00e+00 7.03e-01w  1
     549  0.0000000e+00 3.39e+05 1.37e+17  -1.0 3.00e+04    -  2.42e-01 5.05e-01w  1
    WARNING: Problem in step computation; switching to emergency mode.
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     550  0.0000000e+00 3.50e+04 1.39e+17  -1.0 3.00e+04    -  1.00e+00 3.43e-04h 12
     551  0.0000000e+00 3.50e+04 1.96e+17  -1.0 6.34e+04    -  1.00e+00 3.43e-04h 12
     552  0.0000000e+00 3.50e+04 2.76e+17  -1.0 6.36e+04    -  1.00e+00 3.44e-04h 12
     553  0.0000000e+00 3.50e+04 3.87e+17  -1.0 6.38e+04    -  1.00e+00 3.44e-04h 12
     554  0.0000000e+00 3.49e+04 4.32e+17  -1.0 6.33e+04    -  1.00e+00 3.43e-04h 12
     555  0.0000000e+00 3.49e+04 4.32e+17  -1.0 5.57e+04    -  1.00e+00 3.31e-04h 12
     556  0.0000000e+00 3.49e+04 4.32e+17  -1.0 5.11e+04    -  1.00e+00 6.49e-04h 11
     557  0.0000000e+00 3.49e+04 4.32e+17  -1.0 5.17e+04    -  1.00e+00 6.51e-04h 11
     558  0.0000000e+00 3.49e+04 4.32e+17  -1.0 5.25e+04    -  1.00e+00 6.54e-04h 11
     559  0.0000000e+00 3.48e+04 4.32e+17  -1.0 5.06e+04    -  1.00e+00 6.47e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     560  0.0000000e+00 2.37e+05 1.27e+17  -1.0 5.53e+04    -  1.00e+00 6.79e-01w  1
     561  0.0000000e+00 8.95e+04 7.14e+16  -1.0 2.90e+03    -  1.00e+00 6.32e-01w  1
     562  0.0000000e+00 3.49e+05 1.56e+19  -1.0 2.60e+04    -  1.24e-01 9.90e-01w  1
    WARNING: Problem in step computation; switching to emergency mode.
     563  0.0000000e+00 3.48e+04 4.32e+17  -1.0 2.60e+04    -  1.00e+00 6.63e-04h 11
     564  0.0000000e+00 3.48e+04 4.32e+17  -1.0 5.18e+04    -  1.00e+00 6.51e-04h 11
     565  0.0000000e+00 3.48e+04 4.32e+17  -1.0 5.05e+04    -  1.00e+00 6.47e-04h 11
     566  0.0000000e+00 3.48e+04 4.32e+17  -1.0 5.02e+04    -  1.00e+00 6.46e-04h 11
     567  0.0000000e+00 3.47e+04 4.32e+17  -1.0 5.01e+04    -  1.00e+00 6.46e-04h 11
     568  0.0000000e+00 3.47e+04 4.32e+17  -1.0 5.21e+04    -  1.00e+00 6.52e-04h 11
     569  0.0000000e+00 3.47e+04 4.32e+17  -1.0 5.04e+04    -  1.00e+00 6.47e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     570  0.0000000e+00 3.47e+04 4.32e+17  -1.0 5.01e+04    -  1.00e+00 6.46e-04h 11
     571  0.0000000e+00 3.46e+04 4.32e+17  -1.0 5.09e+04    -  1.00e+00 6.48e-04h 11
     572  0.0000000e+00 3.46e+04 4.32e+17  -1.0 5.15e+04    -  1.00e+00 6.50e-04h 11
     573  0.0000000e+00 2.33e+05 1.19e+17  -1.0 5.19e+04    -  1.00e+00 6.67e-01w  1
     574  0.0000000e+00 9.64e+04 9.86e+16  -1.0 4.42e+03    -  1.00e+00 6.03e-01w  1
     575  0.0000000e+00 2.25e+05 2.81e+18  -1.0 5.20e+04    -  1.33e-01 3.79e-01w  1
     576  0.0000000e+00 3.46e+04 4.32e+17  -1.0 3.03e+04    -  1.00e+00 6.52e-04h 10
     577  0.0000000e+00 3.46e+04 4.32e+17  -1.0 5.02e+04    -  1.00e+00 6.46e-04h 11
     578  0.0000000e+00 3.46e+04 4.32e+17  -1.0 5.13e+04    -  1.00e+00 6.50e-04h 11
     579  0.0000000e+00 3.45e+04 4.32e+17  -1.0 5.00e+04    -  1.00e+00 6.45e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     580  0.0000000e+00 3.45e+04 4.32e+17  -1.0 4.98e+04    -  1.00e+00 6.45e-04h 11
     581  0.0000000e+00 3.45e+04 4.32e+17  -1.0 5.56e+04    -  1.00e+00 6.64e-04h 11
     582  0.0000000e+00 3.45e+04 4.32e+17  -1.0 5.23e+04    -  1.00e+00 6.52e-04h 11
     583  0.0000000e+00 3.44e+04 4.32e+17  -1.0 5.01e+04    -  1.00e+00 6.45e-04h 11
     584  0.0000000e+00 3.44e+04 4.32e+17  -1.0 4.96e+04    -  1.00e+00 6.44e-04h 11
     585  0.0000000e+00 3.44e+04 4.32e+17  -1.0 4.99e+04    -  1.00e+00 6.45e-04h 11
     586  0.0000000e+00 2.27e+05 1.14e+17  -1.0 4.99e+04    -  1.00e+00 6.60e-01w  1
     587  0.0000000e+00 9.83e+04 1.06e+17  -1.0 6.43e+03    -  1.00e+00 5.95e-01w  1
     588  0.0000000e+00 3.76e+05 3.37e+17  -1.0 1.53e+05    -  1.43e-01 1.69e-01w  1
     589  0.0000000e+00 3.44e+04 4.32e+17  -1.0 1.92e+04    -  1.00e+00 6.45e-04h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     590  0.0000000e+00 3.44e+04 4.32e+17  -1.0 5.38e+04    -  1.00e+00 6.58e-04h 11
     591  0.0000000e+00 3.43e+04 4.32e+17  -1.0 5.04e+04    -  1.00e+00 6.46e-04h 11
     592  0.0000000e+00 3.43e+04 4.32e+17  -1.0 5.23e+04    -  1.00e+00 6.53e-04h 11
     593  0.0000000e+00 3.43e+04 4.32e+17  -1.0 4.99e+04    -  1.00e+00 6.45e-04h 11
     594  0.0000000e+00 3.43e+04 4.32e+17  -1.0 5.03e+04    -  1.00e+00 6.46e-04h 11
     595  0.0000000e+00 3.42e+04 4.32e+17  -1.0 5.14e+04    -  1.00e+00 6.49e-04h 11
     596  0.0000000e+00 3.42e+04 4.32e+17  -1.0 4.96e+04    -  1.00e+00 6.44e-04h 11
     597  0.0000000e+00 3.42e+04 4.32e+17  -1.0 4.96e+04    -  1.00e+00 6.44e-04h 11
     598  0.0000000e+00 3.42e+04 4.32e+17  -1.0 5.33e+04    -  1.00e+00 6.56e-04h 11
     599  0.0000000e+00 2.31e+05 1.14e+17  -1.0 5.01e+04    -  1.00e+00 6.61e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     600  0.0000000e+00 1.01e+05 1.12e+17  -1.0 7.62e+03    -  1.00e+00 5.91e-01w  1
     601  0.0000000e+00 4.26e+05 3.38e+17  -1.0 2.59e+05    -  1.51e-01 1.07e-01w  1
     602  0.0000000e+00 3.42e+04 4.32e+17  -1.0 1.46e+04    -  1.00e+00 6.45e-04h 10
     603  0.0000000e+00 3.41e+04 4.32e+17  -1.0 4.94e+04    -  1.00e+00 6.43e-04h 11
     604  0.0000000e+00 3.41e+04 4.32e+17  -1.0 4.91e+04    -  1.00e+00 6.42e-04h 11
     605  0.0000000e+00 3.41e+04 4.32e+17  -1.0 4.96e+04    -  1.00e+00 6.43e-04h 11
     606  0.0000000e+00 3.41e+04 4.32e+17  -1.0 5.00e+04    -  1.00e+00 6.45e-04h 11
     607  0.0000000e+00 3.40e+04 4.32e+17  -1.0 4.92e+04    -  1.00e+00 6.42e-04h 11
     608  0.0000000e+00 3.40e+04 4.32e+17  -1.0 4.97e+04    -  1.00e+00 6.44e-04h 11
     609  0.0000000e+00 3.40e+04 4.32e+17  -1.0 5.10e+04    -  1.00e+00 6.48e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     610  0.0000000e+00 3.40e+04 4.32e+17  -1.0 5.36e+04    -  1.00e+00 6.57e-04h 11
     611  0.0000000e+00 3.40e+04 4.32e+17  -1.0 5.06e+04    -  1.00e+00 6.46e-04h 11
     612  0.0000000e+00 2.23e+05 1.14e+17  -1.0 4.95e+04    -  1.00e+00 6.58e-01w  1
     613  0.0000000e+00 9.66e+04 1.12e+17  -1.0 7.49e+03    -  1.00e+00 5.92e-01w  1
     614  0.0000000e+00 3.99e+05 2.36e+16  -1.0 1.96e+05    -  1.45e-01 1.36e-01w  1
     615  0.0000000e+00 3.39e+04 4.32e+17  -1.0 1.66e+04    -  1.00e+00 6.43e-04h 10
     616  0.0000000e+00 3.39e+04 4.32e+17  -1.0 4.98e+04    -  1.00e+00 6.44e-04h 11
     617  0.0000000e+00 3.39e+04 4.32e+17  -1.0 4.93e+04    -  1.00e+00 6.42e-04h 11
     618  0.0000000e+00 3.39e+04 4.32e+17  -1.0 5.36e+04    -  1.00e+00 6.56e-04h 11
     619  0.0000000e+00 3.38e+04 4.32e+17  -1.0 5.03e+04    -  1.00e+00 6.45e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     620  0.0000000e+00 3.38e+04 4.32e+17  -1.0 4.89e+04    -  1.00e+00 6.41e-04h 11
     621  0.0000000e+00 3.38e+04 4.32e+17  -1.0 4.95e+04    -  1.00e+00 6.43e-04h 11
     622  0.0000000e+00 3.38e+04 4.32e+17  -1.0 4.95e+04    -  1.00e+00 6.43e-04h 11
     623  0.0000000e+00 3.38e+04 4.32e+17  -1.0 4.86e+04    -  1.00e+00 6.40e-04h 11
     624  0.0000000e+00 3.37e+04 4.32e+17  -1.0 4.88e+04    -  1.00e+00 6.41e-04h 11
     625  0.0000000e+00 2.15e+05 1.15e+17  -1.0 4.85e+04    -  1.00e+00 6.55e-01w  1
     626  0.0000000e+00 9.33e+04 1.14e+17  -1.0 7.36e+03    -  1.00e+00 5.93e-01w  1
     627  0.0000000e+00 3.73e+05 2.99e+17  -1.0 1.62e+05    -  1.42e-01 1.64e-01w  1
     628  0.0000000e+00 3.37e+04 4.32e+17  -1.0 1.77e+04    -  1.00e+00 6.40e-04h 10
     629  0.0000000e+00 3.37e+04 4.32e+17  -1.0 4.92e+04    -  1.00e+00 6.42e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     630  0.0000000e+00 3.37e+04 4.32e+17  -1.0 5.19e+04    -  1.00e+00 6.51e-04h 11
     631  0.0000000e+00 3.36e+04 4.32e+17  -1.0 5.10e+04    -  1.00e+00 6.48e-04h 11
     632  0.0000000e+00 3.36e+04 4.32e+17  -1.0 5.07e+04    -  1.00e+00 6.47e-04h 11
     633  0.0000000e+00 3.36e+04 4.32e+17  -1.0 5.09e+04    -  1.00e+00 6.47e-04h 11
     634  0.0000000e+00 3.36e+04 4.32e+17  -1.0 4.92e+04    -  1.00e+00 6.42e-04h 11
     635  0.0000000e+00 3.36e+04 4.32e+17  -1.0 4.81e+04    -  1.00e+00 6.39e-04h 11
     636  0.0000000e+00 3.35e+04 4.32e+17  -1.0 5.37e+04    -  1.00e+00 6.56e-04h 11
     637  0.0000000e+00 3.35e+04 4.32e+17  -1.0 5.55e+04    -  1.00e+00 6.62e-04h 11
     638  0.0000000e+00 2.19e+05 1.30e+17  -1.0 4.93e+04    -  1.00e+00 6.57e-01w  1
     639  0.0000000e+00 9.14e+04 1.13e+17  -1.0 5.63e+03    -  1.00e+00 6.01e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     640  0.0000000e+00 2.40e+05 2.57e+18  -1.0 6.02e+04    -  1.32e-01 3.46e-01w  1
     641  0.0000000e+00 3.35e+04 4.33e+17  -1.0 3.02e+04    -  1.00e+00 6.42e-04h 10
     642  0.0000000e+00 3.35e+04 4.33e+17  -1.0 5.56e+04    -  1.00e+00 6.63e-04h 11
     643  0.0000000e+00 3.35e+04 4.33e+17  -1.0 5.65e+04    -  1.00e+00 6.65e-04h 11
     644  0.0000000e+00 3.34e+04 4.33e+17  -1.0 4.87e+04    -  1.00e+00 6.40e-04h 11
     645  0.0000000e+00 3.34e+04 4.33e+17  -1.0 4.85e+04    -  1.00e+00 6.40e-04h 11
     646  0.0000000e+00 3.34e+04 4.33e+17  -1.0 4.69e+04    -  1.00e+00 6.35e-04h 11
     647  0.0000000e+00 3.34e+04 4.33e+17  -1.0 4.63e+04    -  1.00e+00 6.33e-04h 11
     648  0.0000000e+00 3.33e+04 4.33e+17  -1.0 5.08e+04    -  1.00e+00 6.47e-04h 11
     649  0.0000000e+00 3.33e+04 4.33e+17  -1.0 4.97e+04    -  1.00e+00 6.43e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     650  0.0000000e+00 3.33e+04 4.33e+17  -1.0 4.82e+04    -  1.00e+00 6.38e-04h 11
     651  0.0000000e+00 2.10e+05 1.50e+17  -1.0 4.85e+04    -  1.00e+00 6.55e-01w  1
     652  0.0000000e+00 8.41e+04 1.12e+17  -1.0 3.92e+03    -  1.00e+00 6.10e-01w  1
     653  0.0000000e+00 1.44e+05 4.66e+18  -1.0 3.19e+04    -  1.28e-01 4.89e-01w  1
     654  0.0000000e+00 3.33e+04 4.33e+17  -1.0 4.79e+04    -  1.00e+00 6.40e-04h 10
     655  0.0000000e+00 3.33e+04 4.33e+17  -1.0 4.82e+04    -  1.00e+00 6.38e-04h 11
     656  0.0000000e+00 3.32e+04 4.33e+17  -1.0 5.07e+04    -  1.00e+00 6.46e-04h 11
     657  0.0000000e+00 3.32e+04 4.33e+17  -1.0 4.85e+04    -  1.00e+00 6.39e-04h 11
     658  0.0000000e+00 3.32e+04 4.33e+17  -1.0 5.22e+04    -  1.00e+00 6.51e-04h 11
     659  0.0000000e+00 3.32e+04 4.33e+17  -1.0 5.43e+04    -  1.00e+00 6.58e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     660  0.0000000e+00 3.32e+04 4.33e+17  -1.0 4.87e+04    -  1.00e+00 6.40e-04h 11
     661  0.0000000e+00 3.31e+04 4.33e+17  -1.0 5.46e+04    -  1.00e+00 6.59e-04h 11
     662  0.0000000e+00 3.31e+04 4.33e+17  -1.0 4.83e+04    -  1.00e+00 6.39e-04h 11
     663  0.0000000e+00 3.31e+04 4.33e+17  -1.0 5.48e+04    -  1.00e+00 6.59e-04h 11
     664  0.0000000e+00 2.17e+05 1.51e+17  -1.0 4.91e+04    -  1.00e+00 6.56e-01w  1
     665  0.0000000e+00 8.68e+04 1.15e+17  -1.0 4.17e+03    -  1.00e+00 6.09e-01w  1
     666  0.0000000e+00 1.54e+05 4.46e+18  -1.0 3.36e+04    -  1.27e-01 4.72e-01w  1
     667  0.0000000e+00 3.31e+04 4.33e+17  -1.0 4.78e+04    -  1.00e+00 6.41e-04h 10
     668  0.0000000e+00 3.30e+04 4.33e+17  -1.0 5.28e+04    -  1.00e+00 6.53e-04h 11
     669  0.0000000e+00 3.30e+04 4.33e+17  -1.0 5.20e+04    -  1.00e+00 6.50e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     670  0.0000000e+00 3.30e+04 4.33e+17  -1.0 4.84e+04    -  1.00e+00 6.39e-04h 11
     671  0.0000000e+00 3.30e+04 4.33e+17  -1.0 5.49e+04    -  1.00e+00 6.60e-04h 11
     672  0.0000000e+00 3.30e+04 4.33e+17  -1.0 4.88e+04    -  1.00e+00 6.40e-04h 11
     673  0.0000000e+00 3.29e+04 4.33e+17  -1.0 4.82e+04    -  1.00e+00 6.38e-04h 11
     674  0.0000000e+00 3.29e+04 4.33e+17  -1.0 4.91e+04    -  1.00e+00 6.41e-04h 11
     675  0.0000000e+00 3.29e+04 4.33e+17  -1.0 5.11e+04    -  1.00e+00 6.47e-04h 11
     676  0.0000000e+00 3.29e+04 4.33e+17  -1.0 5.12e+04    -  1.00e+00 6.47e-04h 11
     677  0.0000000e+00 2.09e+05 1.46e+17  -1.0 4.78e+04    -  1.00e+00 6.52e-01w  1
     678  0.0000000e+00 8.53e+04 1.20e+17  -1.0 5.30e+03    -  1.00e+00 6.05e-01w  1
     679  0.0000000e+00 1.86e+05 3.72e+18  -1.0 4.37e+04    -  1.28e-01 4.20e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     680  0.0000000e+00 3.29e+04 4.33e+17  -1.0 4.24e+04    -  1.00e+00 6.37e-04h 10
     681  0.0000000e+00 3.28e+04 4.33e+17  -1.0 6.49e+04    -  1.00e+00 6.94e-04h 11
     682  0.0000000e+00 3.28e+04 4.33e+17  -1.0 5.37e+04    -  1.00e+00 6.55e-04h 11
     683  0.0000000e+00 3.28e+04 4.33e+17  -1.0 4.80e+04    -  1.00e+00 6.37e-04h 11
     684  0.0000000e+00 3.28e+04 4.33e+17  -1.0 5.26e+04    -  1.00e+00 6.52e-04h 11
     685  0.0000000e+00 3.27e+04 4.33e+17  -1.0 4.88e+04    -  1.00e+00 6.40e-04h 11
     686  0.0000000e+00 3.27e+04 4.33e+17  -1.0 4.72e+04    -  1.00e+00 6.35e-04h 11
     687  0.0000000e+00 3.27e+04 4.33e+17  -1.0 4.83e+04    -  1.00e+00 6.38e-04h 11
     688  0.0000000e+00 3.27e+04 4.33e+17  -1.0 4.94e+04    -  1.00e+00 6.42e-04h 11
     689  0.0000000e+00 3.27e+04 4.33e+17  -1.0 5.95e+04    -  1.00e+00 6.75e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     690  0.0000000e+00 2.18e+05 1.51e+17  -1.0 4.95e+04    -  1.00e+00 6.57e-01w  1
     691  0.0000000e+00 8.70e+04 9.42e+16  -1.0 2.26e+03    -  1.00e+00 6.30e-01w  1
     692  0.0000000e+00 8.54e+04 7.61e+18  -1.0 2.24e+04    -  1.21e-01 6.16e-01w  1
     693  0.0000000e+00 3.26e+04 4.33e+17  -1.0 5.36e+04    -  1.00e+00 6.42e-04h 10
     694  0.0000000e+00 3.26e+04 4.33e+17  -1.0 4.92e+04    -  1.00e+00 6.41e-04h 11
     695  0.0000000e+00 3.26e+04 4.34e+17  -1.0 4.83e+04    -  1.00e+00 6.38e-04h 11
     696  0.0000000e+00 3.26e+04 4.34e+17  -1.0 4.73e+04    -  1.00e+00 6.35e-04h 11
     697  0.0000000e+00 3.26e+04 4.34e+17  -1.0 4.70e+04    -  1.00e+00 6.34e-04h 11
     698  0.0000000e+00 3.25e+04 4.34e+17  -1.0 4.97e+04    -  1.00e+00 6.43e-04h 11
     699  0.0000000e+00 3.25e+04 4.34e+17  -1.0 4.70e+04    -  1.00e+00 6.34e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     700  0.0000000e+00 3.25e+04 4.34e+17  -1.0 4.79e+04    -  1.00e+00 6.37e-04h 11
     701  0.0000000e+00 3.25e+04 4.34e+17  -1.0 4.99e+04    -  1.00e+00 6.43e-04h 11
     702  0.0000000e+00 3.25e+04 4.34e+17  -1.0 5.20e+04    -  1.00e+00 6.50e-04h 11
     703  0.0000000e+00 2.12e+05 1.74e+17  -1.0 5.93e+04    -  1.00e+00 6.89e-01w  1
     704  0.0000000e+00 6.24e+04 6.15e+15  -1.0 1.94e+04    -  1.00e+00 7.09e-01w  1
     705  0.0000000e+00 6.87e+06 6.38e+18  -1.0 2.64e+05    -  1.04e-01 4.35e-01w  1
     706  0.0000000e+00 3.24e+04 4.34e+17  -1.0 2.27e+04    -  1.00e+00 6.73e-04h 10
     707  0.0000000e+00 3.24e+04 4.34e+17  -1.0 5.43e+04    -  1.00e+00 6.57e-04h 11
     708  0.0000000e+00 3.24e+04 4.34e+17  -1.0 5.06e+04    -  1.00e+00 6.45e-04h 11
     709  0.0000000e+00 3.24e+04 4.34e+17  -1.0 4.70e+04    -  1.00e+00 6.34e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     710  0.0000000e+00 3.23e+04 4.34e+17  -1.0 4.72e+04    -  1.00e+00 6.35e-04h 11
     711  0.0000000e+00 3.23e+04 4.34e+17  -1.0 4.72e+04    -  1.00e+00 6.35e-04h 11
     712  0.0000000e+00 3.23e+04 4.34e+17  -1.0 5.41e+04    -  1.00e+00 6.56e-04h 11
     713  0.0000000e+00 3.23e+04 4.34e+17  -1.0 5.32e+04    -  1.00e+00 6.53e-04h 11
     714  0.0000000e+00 3.23e+04 4.34e+17  -1.0 4.72e+04    -  1.00e+00 6.34e-04h 11
     715  0.0000000e+00 3.22e+04 4.34e+17  -1.0 4.66e+04    -  1.00e+00 6.33e-04h 11
     716  0.0000000e+00 1.97e+05 1.45e+17  -1.0 4.83e+04    -  1.00e+00 6.53e-01w  1
     717  0.0000000e+00 7.77e+04 1.14e+17  -1.0 3.14e+03    -  1.00e+00 6.17e-01w  1
     718  0.0000000e+00 1.28e+05 5.48e+18  -1.0 2.83e+04    -  1.23e-01 5.12e-01w  1
     719  0.0000000e+00 3.22e+04 4.34e+17  -1.0 4.91e+04    -  1.00e+00 6.38e-04h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     720  0.0000000e+00 3.22e+04 4.34e+17  -1.0 4.84e+04    -  1.00e+00 6.38e-04h 11
     721  0.0000000e+00 3.22e+04 4.34e+17  -1.0 4.72e+04    -  1.00e+00 6.35e-04h 11
     722  0.0000000e+00 3.22e+04 4.34e+17  -1.0 4.89e+04    -  1.00e+00 6.40e-04h 11
     723  0.0000000e+00 3.21e+04 4.34e+17  -1.0 4.83e+04    -  1.00e+00 6.38e-04h 11
     724  0.0000000e+00 3.21e+04 4.34e+17  -1.0 4.72e+04    -  1.00e+00 6.35e-04h 11
     725  0.0000000e+00 3.21e+04 4.34e+17  -1.0 4.76e+04    -  1.00e+00 6.36e-04h 11
     726  0.0000000e+00 3.21e+04 4.34e+17  -1.0 4.75e+04    -  1.00e+00 6.35e-04h 11
     727  0.0000000e+00 3.21e+04 4.34e+17  -1.0 4.70e+04    -  1.00e+00 6.34e-04h 11
     728  0.0000000e+00 3.20e+04 4.34e+17  -1.0 4.69e+04    -  1.00e+00 6.34e-04h 11
     729  0.0000000e+00 1.98e+05 1.65e+17  -1.0 5.65e+04    -  1.00e+00 6.80e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     730  0.0000000e+00 6.64e+04 5.64e+16  -1.0 8.77e+03    -  1.00e+00 6.66e-01w  1
     731  0.0000000e+00 6.53e+06 3.77e+19  -1.0 1.30e+05    -  1.13e-01 9.05e-01w  1
     732  0.0000000e+00 3.20e+04 4.34e+17  -1.0 1.39e+05    -  1.00e+00 6.64e-04h 10
     733  0.0000000e+00 3.20e+04 4.35e+17  -1.0 4.78e+04    -  1.00e+00 6.36e-04h 11
     734  0.0000000e+00 3.20e+04 4.35e+17  -1.0 4.68e+04    -  1.00e+00 6.33e-04h 11
     735  0.0000000e+00 3.20e+04 4.35e+17  -1.0 4.74e+04    -  1.00e+00 6.35e-04h 11
     736  0.0000000e+00 3.19e+04 4.35e+17  -1.0 4.59e+04    -  1.00e+00 6.31e-04h 11
     737  0.0000000e+00 3.19e+04 4.35e+17  -1.0 4.60e+04    -  1.00e+00 6.31e-04h 11
     738  0.0000000e+00 3.19e+04 4.35e+17  -1.0 4.58e+04    -  1.00e+00 6.31e-04h 11
     739  0.0000000e+00 3.19e+04 4.35e+17  -1.0 4.64e+04    -  1.00e+00 6.32e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     740  0.0000000e+00 3.19e+04 4.35e+17  -1.0 4.88e+04    -  1.00e+00 6.39e-04h 11
     741  0.0000000e+00 3.18e+04 4.35e+17  -1.0 4.61e+04    -  1.00e+00 6.31e-04h 11
     742  0.0000000e+00 1.91e+05 1.37e+17  -1.0 4.56e+04    -  1.00e+00 6.45e-01w  1
     743  0.0000000e+00 7.98e+04 1.31e+17  -1.0 6.43e+03    -  1.00e+00 6.03e-01w  1
     744  0.0000000e+00 2.20e+05 3.14e+18  -1.0 5.47e+04    -  1.27e-01 3.60e-01w  1
     745  0.0000000e+00 3.18e+04 4.35e+17  -1.0 3.69e+04    -  1.00e+00 6.30e-04h 10
     746  0.0000000e+00 3.18e+04 4.35e+17  -1.0 4.78e+04    -  1.00e+00 6.36e-04h 11
     747  0.0000000e+00 3.18e+04 4.35e+17  -1.0 4.87e+04    -  1.00e+00 6.39e-04h 11
     748  0.0000000e+00 3.18e+04 4.35e+17  -1.0 5.11e+04    -  1.00e+00 6.46e-04h 11
     749  0.0000000e+00 3.17e+04 4.35e+17  -1.0 4.87e+04    -  1.00e+00 6.39e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     750  0.0000000e+00 3.17e+04 4.35e+17  -1.0 4.64e+04    -  1.00e+00 6.32e-04h 11
     751  0.0000000e+00 3.17e+04 4.35e+17  -1.0 4.58e+04    -  1.00e+00 6.30e-04h 11
     752  0.0000000e+00 3.17e+04 4.35e+17  -1.0 4.90e+04    -  1.00e+00 6.40e-04h 11
     753  0.0000000e+00 3.17e+04 4.35e+17  -1.0 4.91e+04    -  1.00e+00 6.40e-04h 11
     754  0.0000000e+00 3.16e+04 4.35e+17  -1.0 4.94e+04    -  1.00e+00 6.41e-04h 11
     755  0.0000000e+00 1.95e+05 1.51e+17  -1.0 5.16e+04    -  1.00e+00 6.63e-01w  1
     756  0.0000000e+00 7.11e+04 9.60e+16  -1.0 1.48e+03    -  1.00e+00 6.39e-01w  1
     757  0.0000000e+00 3.72e+04 1.50e+19  -1.0 1.21e+03    -  1.15e-01 8.77e-01w  1
     758  0.0000000e+00 3.16e+04 4.35e+17  -1.0 3.04e+04    -  1.00e+00 6.48e-04h 10
     759  0.0000000e+00 3.16e+04 4.35e+17  -1.0 4.96e+04    -  1.00e+00 6.41e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     760  0.0000000e+00 3.16e+04 4.35e+17  -1.0 4.61e+04    -  1.00e+00 6.31e-04h 11
     761  0.0000000e+00 3.16e+04 4.35e+17  -1.0 5.67e+04    -  1.00e+00 6.64e-04h 11
     762  0.0000000e+00 3.15e+04 4.35e+17  -1.0 5.45e+04    -  1.00e+00 6.56e-04h 11
     763  0.0000000e+00 3.15e+04 4.36e+17  -1.0 4.90e+04    -  1.00e+00 6.39e-04h 11
     764  0.0000000e+00 3.15e+04 4.36e+17  -1.0 4.61e+04    -  1.00e+00 6.31e-04h 11
     765  0.0000000e+00 3.15e+04 4.36e+17  -1.0 5.21e+04    -  1.00e+00 6.49e-04h 11
     766  0.0000000e+00 3.15e+04 4.36e+17  -1.0 4.95e+04    -  1.00e+00 6.41e-04h 11
     767  0.0000000e+00 3.14e+04 4.36e+17  -1.0 4.57e+04    -  1.00e+00 6.30e-04h 11
     768  0.0000000e+00 1.87e+05 1.35e+17  -1.0 4.56e+04    -  1.00e+00 6.45e-01w  1
     769  0.0000000e+00 7.79e+04 1.33e+17  -1.0 6.27e+03    -  1.00e+00 6.05e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     770  0.0000000e+00 2.18e+05 3.27e+18  -1.0 5.49e+04    -  1.25e-01 3.61e-01w  1
     771  0.0000000e+00 3.14e+04 4.36e+17  -1.0 3.68e+04    -  1.00e+00 6.30e-04h 10
     772  0.0000000e+00 3.14e+04 4.36e+17  -1.0 4.57e+04    -  1.00e+00 6.30e-04h 11
     773  0.0000000e+00 3.14e+04 4.36e+17  -1.0 5.00e+04    -  1.00e+00 6.43e-04h 11
     774  0.0000000e+00 3.14e+04 4.36e+17  -1.0 4.59e+04    -  1.00e+00 6.30e-04h 11
     775  0.0000000e+00 3.13e+04 4.36e+17  -1.0 5.12e+04    -  1.00e+00 6.46e-04h 11
     776  0.0000000e+00 3.13e+04 4.36e+17  -1.0 6.25e+04    -  1.00e+00 6.83e-04h 11
     777  0.0000000e+00 3.13e+04 4.36e+17  -1.0 4.77e+04    -  1.00e+00 6.35e-04h 11
     778  0.0000000e+00 3.13e+04 4.36e+17  -1.0 4.91e+04    -  1.00e+00 6.40e-04h 11
     779  0.0000000e+00 3.13e+04 4.36e+17  -1.0 4.61e+04    -  1.00e+00 6.31e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     780  0.0000000e+00 3.12e+04 4.36e+17  -1.0 4.69e+04    -  1.00e+00 6.33e-04h 11
     781  0.0000000e+00 1.87e+05 1.49e+17  -1.0 5.15e+04    -  1.00e+00 6.63e-01w  1
     782  0.0000000e+00 6.76e+04 9.63e+16  -1.0 1.83e+03    -  1.00e+00 6.42e-01w  1
     783  0.0000000e+00 4.50e+04 1.71e+19  -1.0 2.80e+03    -  1.13e-01 9.23e-01w  1
     784  0.0000000e+00 3.12e+04 4.36e+17  -1.0 6.39e+04    -  1.00e+00 6.47e-04h 10
     785  0.0000000e+00 3.12e+04 4.36e+17  -1.0 4.70e+04    -  1.00e+00 6.33e-04h 11
     786  0.0000000e+00 3.12e+04 4.36e+17  -1.0 4.85e+04    -  1.00e+00 6.38e-04h 11
     787  0.0000000e+00 3.12e+04 4.36e+17  -1.0 4.60e+04    -  1.00e+00 6.31e-04h 11
     788  0.0000000e+00 3.11e+04 4.36e+17  -1.0 6.07e+04    -  1.00e+00 6.76e-04h 11
     789  0.0000000e+00 3.11e+04 4.36e+17  -1.0 4.66e+04    -  1.00e+00 6.32e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     790  0.0000000e+00 3.11e+04 4.37e+17  -1.0 4.67e+04    -  1.00e+00 6.32e-04h 11
     791  0.0000000e+00 3.11e+04 4.37e+17  -1.0 4.62e+04    -  1.00e+00 6.31e-04h 11
     792  0.0000000e+00 3.11e+04 4.37e+17  -1.0 4.51e+04    -  1.00e+00 6.28e-04h 11
     793  0.0000000e+00 3.10e+04 4.37e+17  -1.0 4.60e+04    -  1.00e+00 6.31e-04h 11
     794  0.0000000e+00 1.82e+05 1.32e+17  -1.0 4.45e+04    -  1.00e+00 6.41e-01w  1
     795  0.0000000e+00 7.76e+04 1.35e+17  -1.0 6.27e+03    -  1.00e+00 6.07e-01w  1
     796  0.0000000e+00 2.38e+05 2.91e+18  -1.0 6.40e+04    -  1.25e-01 3.27e-01w  1
     797  0.0000000e+00 3.10e+04 4.37e+17  -1.0 3.35e+04    -  1.00e+00 6.26e-04h 10
     798  0.0000000e+00 3.10e+04 4.37e+17  -1.0 4.59e+04    -  1.00e+00 6.30e-04h 11
     799  0.0000000e+00 3.10e+04 4.37e+17  -1.0 4.57e+04    -  1.00e+00 6.30e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     800  0.0000000e+00 3.10e+04 4.37e+17  -1.0 4.50e+04    -  1.00e+00 6.28e-04h 11
     801  0.0000000e+00 3.09e+04 4.37e+17  -1.0 4.48e+04    -  1.00e+00 6.27e-04h 11
     802  0.0000000e+00 3.09e+04 4.37e+17  -1.0 4.43e+04    -  1.00e+00 6.26e-04h 11
     803  0.0000000e+00 3.09e+04 4.37e+17  -1.0 5.26e+04    -  1.00e+00 6.50e-04h 11
     804  0.0000000e+00 3.09e+04 4.37e+17  -1.0 5.26e+04    -  1.00e+00 6.50e-04h 11
     805  0.0000000e+00 3.09e+04 4.37e+17  -1.0 5.01e+04    -  1.00e+00 6.42e-04h 11
     806  0.0000000e+00 3.08e+04 4.37e+17  -1.0 4.52e+04    -  1.00e+00 6.28e-04h 11
     807  0.0000000e+00 1.81e+05 1.31e+17  -1.0 4.45e+04    -  1.00e+00 6.41e-01w  1
     808  0.0000000e+00 7.71e+04 1.38e+17  -1.0 6.65e+03    -  1.00e+00 6.06e-01w  1
     809  0.0000000e+00 2.42e+05 2.81e+18  -1.0 6.56e+04    -  1.25e-01 3.19e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     810  0.0000000e+00 3.08e+04 4.37e+17  -1.0 3.29e+04    -  1.00e+00 6.26e-04h 10
     811  0.0000000e+00 3.08e+04 4.37e+17  -1.0 4.44e+04    -  1.00e+00 6.26e-04h 11
     812  0.0000000e+00 3.08e+04 4.37e+17  -1.0 4.48e+04    -  1.00e+00 6.27e-04h 11
     813  0.0000000e+00 3.08e+04 4.37e+17  -1.0 4.51e+04    -  1.00e+00 6.28e-04h 11
     814  0.0000000e+00 3.07e+04 4.37e+17  -1.0 4.40e+04    -  1.00e+00 6.25e-04h 11
     815  0.0000000e+00 3.07e+04 4.38e+17  -1.0 4.48e+04    -  1.00e+00 6.27e-04h 11
     816  0.0000000e+00 3.07e+04 4.38e+17  -1.0 4.43e+04    -  1.00e+00 6.26e-04h 11
     817  0.0000000e+00 3.07e+04 4.38e+17  -1.0 4.66e+04    -  1.00e+00 6.32e-04h 11
     818  0.0000000e+00 3.07e+04 4.38e+17  -1.0 4.40e+04    -  1.00e+00 6.25e-04h 11
     819  0.0000000e+00 3.06e+04 4.38e+17  -1.0 5.59e+04    -  1.00e+00 6.60e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     820  0.0000000e+00 1.90e+05 1.33e+17  -1.0 4.54e+04    -  1.00e+00 6.43e-01w  1
     821  0.0000000e+00 8.13e+04 1.42e+17  -1.0 6.98e+03    -  1.00e+00 6.04e-01w  1
     822  0.0000000e+00 2.56e+05 2.56e+18  -1.0 7.02e+04    -  1.24e-01 3.00e-01w  1
     823  0.0000000e+00 3.06e+04 4.38e+17  -1.0 3.23e+04    -  1.00e+00 6.28e-04h 10
     824  0.0000000e+00 3.06e+04 4.38e+17  -1.0 4.56e+04    -  1.00e+00 6.29e-04h 11
     825  0.0000000e+00 3.06e+04 4.38e+17  -1.0 4.46e+04    -  1.00e+00 6.26e-04h 11
     826  0.0000000e+00 3.06e+04 4.38e+17  -1.0 4.38e+04    -  1.00e+00 6.24e-04h 11
     827  0.0000000e+00 3.05e+04 4.38e+17  -1.0 4.46e+04    -  1.00e+00 6.26e-04h 11
     828  0.0000000e+00 3.05e+04 4.38e+17  -1.0 4.45e+04    -  1.00e+00 6.26e-04h 11
     829  0.0000000e+00 3.05e+04 4.38e+17  -1.0 4.36e+04    -  1.00e+00 6.24e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     830  0.0000000e+00 3.05e+04 4.38e+17  -1.0 4.35e+04    -  1.00e+00 6.23e-04h 11
     831  0.0000000e+00 3.05e+04 4.38e+17  -1.0 4.69e+04    -  1.00e+00 6.33e-04h 11
     832  0.0000000e+00 3.05e+04 4.38e+17  -1.0 4.39e+04    -  1.00e+00 6.24e-04h 11
     833  0.0000000e+00 1.74e+05 1.29e+17  -1.0 4.41e+04    -  1.00e+00 6.40e-01w  1
     834  0.0000000e+00 7.40e+04 1.42e+17  -1.0 6.89e+03    -  1.00e+00 6.06e-01w  1
     835  0.0000000e+00 2.42e+05 2.87e+18  -1.0 6.55e+04    -  1.23e-01 3.17e-01w  1
     836  0.0000000e+00 3.04e+04 4.38e+17  -1.0 3.50e+04    -  1.00e+00 6.25e-04h 10
     837  0.0000000e+00 3.04e+04 4.38e+17  -1.0 4.49e+04    -  1.00e+00 6.27e-04h 11
     838  0.0000000e+00 3.04e+04 4.38e+17  -1.0 4.35e+04    -  1.00e+00 6.23e-04h 11
     839  0.0000000e+00 3.04e+04 4.39e+17  -1.0 4.37e+04    -  1.00e+00 6.24e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     840  0.0000000e+00 3.04e+04 4.39e+17  -1.0 4.34e+04    -  1.00e+00 6.23e-04h 11
     841  0.0000000e+00 3.03e+04 4.39e+17  -1.0 4.36e+04    -  1.00e+00 6.24e-04h 11
     842  0.0000000e+00 3.03e+04 4.39e+17  -1.0 4.37e+04    -  1.00e+00 6.24e-04h 11
     843  0.0000000e+00 3.03e+04 4.39e+17  -1.0 4.65e+04    -  1.00e+00 6.32e-04h 11
     844  0.0000000e+00 3.03e+04 4.39e+17  -1.0 7.39e+04    -  1.00e+00 7.22e-04h 11
     845  0.0000000e+00 3.03e+04 4.39e+17  -1.0 4.95e+04    -  1.00e+00 6.38e-04h 11
     846  0.0000000e+00 1.68e+05 1.25e+17  -1.0 4.30e+04    -  1.00e+00 6.37e-01w  1
     847  0.0000000e+00 7.33e+04 1.48e+17  -1.0 7.87e+03    -  1.00e+00 6.03e-01w  1
     848  0.0000000e+00 2.63e+05 2.32e+18  -1.0 7.83e+04    -  1.25e-01 2.80e-01w  1
     849  0.0000000e+00 3.02e+04 4.39e+17  -1.0 3.02e+04    -  1.00e+00 6.22e-04h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     850  0.0000000e+00 3.02e+04 4.39e+17  -1.0 4.36e+04    -  1.00e+00 6.24e-04h 11
     851  0.0000000e+00 3.02e+04 4.39e+17  -1.0 5.18e+04    -  1.00e+00 6.47e-04h 11
     852  0.0000000e+00 3.02e+04 4.39e+17  -1.0 4.79e+04    -  1.00e+00 6.36e-04h 11
     853  0.0000000e+00 3.02e+04 4.39e+17  -1.0 4.37e+04    -  1.00e+00 6.24e-04h 11
     854  0.0000000e+00 3.01e+04 4.39e+17  -1.0 4.33e+04    -  1.00e+00 6.23e-04h 11
     855  0.0000000e+00 3.01e+04 4.39e+17  -1.0 4.37e+04    -  1.00e+00 6.24e-04h 11
     856  0.0000000e+00 3.01e+04 4.39e+17  -1.0 4.73e+04    -  1.00e+00 6.34e-04h 11
     857  0.0000000e+00 3.01e+04 4.40e+17  -1.0 4.60e+04    -  1.00e+00 6.30e-04h 11
     858  0.0000000e+00 3.01e+04 4.40e+17  -1.0 4.61e+04    -  1.00e+00 6.30e-04h 11
     859  0.0000000e+00 1.73e+05 1.29e+17  -1.0 4.45e+04    -  1.00e+00 6.41e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     860  0.0000000e+00 7.34e+04 1.39e+17  -1.0 5.58e+03    -  1.00e+00 6.13e-01w  1
     861  0.0000000e+00 2.29e+05 3.35e+18  -1.0 6.05e+04    -  1.20e-01 3.35e-01w  1
     862  0.0000000e+00 3.01e+04 4.40e+17  -1.0 3.47e+04    -  1.00e+00 6.26e-04h 10
     863  0.0000000e+00 3.00e+04 4.40e+17  -1.0 4.45e+04    -  1.00e+00 6.26e-04h 11
     864  0.0000000e+00 3.00e+04 4.40e+17  -1.0 4.36e+04    -  1.00e+00 6.24e-04h 11
     865  0.0000000e+00 3.00e+04 4.40e+17  -1.0 4.28e+04    -  1.00e+00 6.21e-04h 11
     866  0.0000000e+00 3.00e+04 4.40e+17  -1.0 4.76e+04    -  1.00e+00 6.35e-04h 11
     867  0.0000000e+00 3.00e+04 4.40e+17  -1.0 4.46e+04    -  1.00e+00 6.26e-04h 11
     868  0.0000000e+00 2.99e+04 4.40e+17  -1.0 4.42e+04    -  1.00e+00 6.25e-04h 11
     869  0.0000000e+00 2.99e+04 4.40e+17  -1.0 4.56e+04    -  1.00e+00 6.29e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     870  0.0000000e+00 2.99e+04 4.40e+17  -1.0 4.33e+04    -  1.00e+00 6.23e-04h 11
     871  0.0000000e+00 2.99e+04 4.40e+17  -1.0 4.58e+04    -  1.00e+00 6.30e-04h 11
     872  0.0000000e+00 1.70e+05 1.25e+17  -1.0 4.34e+04    -  1.00e+00 6.38e-01w  1
     873  0.0000000e+00 7.39e+04 1.50e+17  -1.0 7.57e+03    -  1.00e+00 6.05e-01w  1
     874  0.0000000e+00 2.65e+05 2.37e+18  -1.0 7.88e+04    -  1.23e-01 2.76e-01w  1
     875  0.0000000e+00 2.99e+04 4.40e+17  -1.0 2.98e+04    -  1.00e+00 6.23e-04h 10
     876  0.0000000e+00 2.98e+04 4.40e+17  -1.0 4.41e+04    -  1.00e+00 6.25e-04h 11
     877  0.0000000e+00 2.98e+04 4.40e+17  -1.0 4.48e+04    -  1.00e+00 6.27e-04h 11
     878  0.0000000e+00 2.98e+04 4.41e+17  -1.0 4.27e+04    -  1.00e+00 6.21e-04h 11
     879  0.0000000e+00 2.98e+04 4.41e+17  -1.0 4.25e+04    -  1.00e+00 6.20e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     880  0.0000000e+00 2.98e+04 4.41e+17  -1.0 4.31e+04    -  1.00e+00 6.22e-04h 11
     881  0.0000000e+00 2.98e+04 4.41e+17  -1.0 5.32e+04    -  1.00e+00 6.51e-04h 11
     882  0.0000000e+00 2.97e+04 4.41e+17  -1.0 4.41e+04    -  1.00e+00 6.24e-04h 11
     883  0.0000000e+00 2.97e+04 4.41e+17  -1.0 4.36e+04    -  1.00e+00 6.23e-04h 11
     884  0.0000000e+00 2.97e+04 4.41e+17  -1.0 4.26e+04    -  1.00e+00 6.21e-04h 11
     885  0.0000000e+00 1.63e+05 1.49e+17  -1.0 5.46e+04    -  1.00e+00 6.71e-01w  1
     886  0.0000000e+00 4.97e+04 2.37e+16  -1.0 1.62e+04    -  1.00e+00 7.10e-01w  1
     887  0.0000000e+00 4.31e+06 2.31e+19  -1.0 8.88e+04    -  9.65e-02 9.90e-01w  1
     888  0.0000000e+00 2.97e+04 4.41e+17  -1.0 3.79e+04    -  1.00e+00 6.56e-04h 10
     889  0.0000000e+00 2.97e+04 4.41e+17  -1.0 4.35e+04    -  1.00e+00 6.23e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     890  0.0000000e+00 2.96e+04 4.41e+17  -1.0 4.74e+04    -  1.00e+00 6.34e-04h 11
     891  0.0000000e+00 2.96e+04 4.41e+17  -1.0 4.31e+04    -  1.00e+00 6.22e-04h 11
     892  0.0000000e+00 2.96e+04 4.41e+17  -1.0 4.28e+04    -  1.00e+00 6.21e-04h 11
     893  0.0000000e+00 2.96e+04 4.41e+17  -1.0 8.71e+04    -  1.00e+00 7.73e-04h 11
     894  0.0000000e+00 2.96e+04 4.41e+17  -1.0 6.02e+04    -  1.00e+00 3.34e-04h 12
     895  0.0000000e+00 2.95e+04 4.41e+17  -1.0 6.25e+04    -  1.00e+00 6.80e-04h 11
     896  0.0000000e+00 2.95e+04 4.42e+17  -1.0 4.64e+04    -  1.00e+00 6.30e-04h 11
     897  0.0000000e+00 2.95e+04 4.42e+17  -1.0 5.51e+04    -  1.00e+00 6.57e-04h 11
     898  0.0000000e+00 1.76e+05 1.26e+17  -1.0 4.38e+04    -  1.00e+00 6.38e-01w  1
     899  0.0000000e+00 7.73e+04 1.56e+17  -1.0 8.02e+03    -  1.00e+00 6.04e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     900  0.0000000e+00 2.81e+05 2.02e+18  -1.0 8.86e+04    -  1.23e-01 2.50e-01w  1
     901  0.0000000e+00 2.95e+04 4.42e+17  -1.0 2.69e+04    -  1.00e+00 6.24e-04h 10
     902  0.0000000e+00 2.95e+04 4.42e+17  -1.0 4.34e+04    -  1.00e+00 6.23e-04h 11
     903  0.0000000e+00 2.95e+04 4.42e+17  -1.0 4.21e+04    -  1.00e+00 6.19e-04h 11
     904  0.0000000e+00 2.94e+04 4.42e+17  -1.0 4.64e+04    -  1.00e+00 6.31e-04h 11
     905  0.0000000e+00 2.94e+04 4.42e+17  -1.0 4.53e+04    -  1.00e+00 6.28e-04h 11
     906  0.0000000e+00 2.94e+04 4.42e+17  -1.0 4.73e+04    -  1.00e+00 6.34e-04h 11
     907  0.0000000e+00 2.94e+04 4.42e+17  -1.0 4.28e+04    -  1.00e+00 6.21e-04h 11
     908  0.0000000e+00 2.94e+04 4.42e+17  -1.0 4.23e+04    -  1.00e+00 6.20e-04h 11
     909  0.0000000e+00 2.93e+04 4.42e+17  -1.0 4.61e+04    -  1.00e+00 6.30e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     910  0.0000000e+00 2.93e+04 4.42e+17  -1.0 4.22e+04    -  1.00e+00 6.19e-04h 11
     911  0.0000000e+00 1.60e+05 1.31e+17  -1.0 4.65e+04    -  1.00e+00 6.47e-01w  1
     912  0.0000000e+00 6.28e+04 1.30e+17  -1.0 2.49e+03    -  1.00e+00 6.31e-01w  1
     913  0.0000000e+00 1.40e+05 6.32e+18  -1.0 3.17e+04    -  1.12e-01 4.80e-01w  1
     914  0.0000000e+00 2.93e+04 4.42e+17  -1.0 4.90e+04    -  1.00e+00 6.31e-04h 10
     915  0.0000000e+00 2.93e+04 4.43e+17  -1.0 4.39e+04    -  1.00e+00 6.24e-04h 11
     916  0.0000000e+00 2.93e+04 4.43e+17  -1.0 4.55e+04    -  1.00e+00 6.29e-04h 11
     917  0.0000000e+00 2.93e+04 4.43e+17  -1.0 4.41e+04    -  1.00e+00 6.25e-04h 11
     918  0.0000000e+00 2.92e+04 4.43e+17  -1.0 4.20e+04    -  1.00e+00 6.19e-04h 11
     919  0.0000000e+00 2.92e+04 4.43e+17  -1.0 4.90e+04    -  1.00e+00 6.39e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     920  0.0000000e+00 2.92e+04 4.43e+17  -1.0 4.22e+04    -  1.00e+00 6.19e-04h 11
     921  0.0000000e+00 2.92e+04 4.43e+17  -1.0 4.17e+04    -  1.00e+00 6.18e-04h 11
     922  0.0000000e+00 2.92e+04 4.43e+17  -1.0 4.45e+04    -  1.00e+00 6.26e-04h 11
     923  0.0000000e+00 2.91e+04 4.43e+17  -1.0 4.18e+04    -  1.00e+00 6.18e-04h 11
     924  0.0000000e+00 1.58e+05 1.23e+17  -1.0 4.34e+04    -  1.00e+00 6.38e-01w  1
     925  0.0000000e+00 6.75e+04 1.50e+17  -1.0 6.24e+03    -  1.00e+00 6.15e-01w  1
     926  0.0000000e+00 2.40e+05 3.21e+18  -1.0 6.68e+04    -  1.18e-01 3.12e-01w  1
     927  0.0000000e+00 2.91e+04 4.43e+17  -1.0 3.39e+04    -  1.00e+00 6.23e-04h 10
     928  0.0000000e+00 2.91e+04 4.43e+17  -1.0 4.15e+04    -  1.00e+00 6.18e-04h 11
     929  0.0000000e+00 2.91e+04 4.43e+17  -1.0 4.18e+04    -  1.00e+00 6.18e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     930  0.0000000e+00 2.91e+04 4.43e+17  -1.0 4.45e+04    -  1.00e+00 6.26e-04h 11
     931  0.0000000e+00 2.91e+04 4.44e+17  -1.0 4.17e+04    -  1.00e+00 6.18e-04h 11
     932  0.0000000e+00 2.90e+04 4.44e+17  -1.0 4.34e+04    -  1.00e+00 6.23e-04h 11
     933  0.0000000e+00 2.90e+04 4.44e+17  -1.0 4.14e+04    -  1.00e+00 6.18e-04h 11
     934  0.0000000e+00 2.90e+04 4.44e+17  -1.0 4.12e+04    -  1.00e+00 6.17e-04h 11
     935  0.0000000e+00 2.90e+04 4.44e+17  -1.0 4.52e+04    -  1.00e+00 6.28e-04h 11
     936  0.0000000e+00 2.90e+04 4.44e+17  -1.0 4.20e+04    -  1.00e+00 6.19e-04h 11
     937  0.0000000e+00 1.54e+05 1.40e+17  -1.0 5.12e+04    -  1.00e+00 6.60e-01w  1
     938  0.0000000e+00 5.24e+04 9.35e+16  -1.0 4.76e+03    -  1.00e+00 6.64e-01w  1
     939  0.0000000e+00 1.43e+05 2.37e+19  -1.0 1.92e+04    -  1.03e-01 9.90e-01w  1
    WARNING: Problem in step computation; switching to emergency mode.
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     940  0.0000000e+00 2.89e+04 4.44e+17  -1.0 1.92e+04  20.0 1.00e+00 6.45e-04h 11
     941  0.0000000e+00 2.89e+04 4.44e+17  -1.0 4.42e+04    -  1.00e+00 6.25e-04h 11
     942  0.0000000e+00 2.89e+04 4.44e+17  -1.0 4.64e+04    -  1.00e+00 6.31e-04h 11
     943  0.0000000e+00 2.89e+04 4.44e+17  -1.0 4.17e+04    -  1.00e+00 6.18e-04h 11
     944  0.0000000e+00 2.89e+04 4.44e+17  -1.0 4.17e+04    -  1.00e+00 6.18e-04h 11
     945  0.0000000e+00 2.89e+04 4.44e+17  -1.0 4.11e+04    -  1.00e+00 6.17e-04h 11
     946  0.0000000e+00 2.88e+04 4.45e+17  -1.0 4.55e+04    -  1.00e+00 6.29e-04h 11
     947  0.0000000e+00 2.88e+04 4.45e+17  -1.0 4.22e+04    -  1.00e+00 6.19e-04h 11
     948  0.0000000e+00 2.88e+04 4.45e+17  -1.0 4.74e+04    -  1.00e+00 6.34e-04h 11
     949  0.0000000e+00 2.88e+04 4.45e+17  -1.0 4.15e+04    -  1.00e+00 6.18e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     950  0.0000000e+00 1.55e+05 1.26e+17  -1.0 4.50e+04    -  1.00e+00 6.42e-01w  1
     951  0.0000000e+00 6.32e+04 1.41e+17  -1.0 3.91e+03    -  1.00e+00 6.27e-01w  1
     952  0.0000000e+00 1.87e+05 4.95e+18  -1.0 4.52e+04    -  1.12e-01 3.97e-01w  1
     953  0.0000000e+00 2.88e+04 4.45e+17  -1.0 4.33e+04    -  1.00e+00 6.27e-04h 10
     954  0.0000000e+00 2.87e+04 4.45e+17  -1.0 4.27e+04    -  1.00e+00 6.21e-04h 11
     955  0.0000000e+00 2.87e+04 4.45e+17  -1.0 4.10e+04    -  1.00e+00 6.16e-04h 11
     956  0.0000000e+00 2.87e+04 4.45e+17  -1.0 4.11e+04    -  1.00e+00 6.17e-04h 11
     957  0.0000000e+00 2.87e+04 4.45e+17  -1.0 4.27e+04    -  1.00e+00 6.21e-04h 11
     958  0.0000000e+00 2.87e+04 4.45e+17  -1.0 4.08e+04    -  1.00e+00 6.16e-04h 11
     959  0.0000000e+00 2.87e+04 4.45e+17  -1.0 4.06e+04    -  1.00e+00 6.15e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     960  0.0000000e+00 2.86e+04 4.45e+17  -1.0 4.15e+04    -  1.00e+00 6.18e-04h 11
     961  0.0000000e+00 2.86e+04 4.46e+17  -1.0 4.06e+04    -  1.00e+00 6.15e-04h 11
     962  0.0000000e+00 2.86e+04 4.46e+17  -1.0 4.05e+04    -  1.00e+00 6.15e-04h 11
     963  0.0000000e+00 1.51e+05 1.16e+17  -1.0 4.09e+04    -  1.00e+00 6.31e-01w  1
     964  0.0000000e+00 6.94e+04 1.66e+17  -1.0 8.70e+03    -  1.00e+00 6.06e-01w  1
     965  0.0000000e+00 2.97e+05 1.48e+18  -1.0 1.12e+05    -  1.24e-01 2.07e-01w  1
     966  0.0000000e+00 2.86e+04 4.46e+17  -1.0 2.40e+04    -  1.00e+00 6.16e-04h 10
     967  0.0000000e+00 2.86e+04 4.46e+17  -1.0 4.19e+04    -  1.00e+00 6.19e-04h 11
     968  0.0000000e+00 2.86e+04 4.46e+17  -1.0 4.24e+04    -  1.00e+00 6.20e-04h 11
     969  0.0000000e+00 2.85e+04 4.46e+17  -1.0 4.62e+04    -  1.00e+00 6.31e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     970  0.0000000e+00 2.85e+04 4.46e+17  -1.0 4.10e+04    -  1.00e+00 6.16e-04h 11
     971  0.0000000e+00 2.85e+04 4.46e+17  -1.0 4.08e+04    -  1.00e+00 6.16e-04h 11
     972  0.0000000e+00 2.85e+04 4.46e+17  -1.0 4.04e+04    -  1.00e+00 6.15e-04h 11
     973  0.0000000e+00 2.85e+04 4.46e+17  -1.0 4.03e+04    -  1.00e+00 6.15e-04h 11
     974  0.0000000e+00 2.84e+04 4.46e+17  -1.0 4.06e+04    -  1.00e+00 6.15e-04h 11
     975  0.0000000e+00 2.84e+04 4.47e+17  -1.0 4.19e+04    -  1.00e+00 6.19e-04h 11
     976  0.0000000e+00 1.51e+05 1.15e+17  -1.0 4.04e+04    -  1.00e+00 6.30e-01w  1
     977  0.0000000e+00 7.05e+04 1.73e+17  -1.0 9.72e+03    -  1.00e+00 6.02e-01w  1
     978  0.0000000e+00 3.15e+05 9.10e+17  -1.0 1.36e+05    -  1.27e-01 1.74e-01w  1
     979  0.0000000e+00 2.84e+04 4.47e+17  -1.0 2.39e+04    -  1.00e+00 6.15e-04h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     980  0.0000000e+00 2.84e+04 4.47e+17  -1.0 4.06e+04    -  1.00e+00 6.15e-04h 11
     981  0.0000000e+00 2.84e+04 4.47e+17  -1.0 4.14e+04    -  1.00e+00 6.17e-04h 11
     982  0.0000000e+00 2.84e+04 4.47e+17  -1.0 4.29e+04    -  1.00e+00 6.22e-04h 11
     983  0.0000000e+00 2.83e+04 4.47e+17  -1.0 4.18e+04    -  1.00e+00 6.18e-04h 11
     984  0.0000000e+00 2.83e+04 4.47e+17  -1.0 4.22e+04    -  1.00e+00 6.20e-04h 11
     985  0.0000000e+00 2.83e+04 4.47e+17  -1.0 4.08e+04    -  1.00e+00 6.16e-04h 11
     986  0.0000000e+00 2.83e+04 4.47e+17  -1.0 4.07e+04    -  1.00e+00 6.16e-04h 11
     987  0.0000000e+00 2.83e+04 4.47e+17  -1.0 4.26e+04    -  1.00e+00 6.21e-04h 11
     988  0.0000000e+00 2.83e+04 4.47e+17  -1.0 4.22e+04    -  1.00e+00 6.20e-04h 11
     989  0.0000000e+00 1.48e+05 1.24e+17  -1.0 4.47e+04    -  1.00e+00 6.41e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     990  0.0000000e+00 6.03e+04 1.45e+17  -1.0 3.76e+03    -  1.00e+00 6.30e-01w  1
     991  0.0000000e+00 1.87e+05 5.17e+18  -1.0 4.56e+04    -  1.10e-01 3.95e-01w  1
     992  0.0000000e+00 2.82e+04 4.48e+17  -1.0 4.55e+04    -  1.00e+00 6.26e-04h 10
     993  0.0000000e+00 2.82e+04 4.48e+17  -1.0 4.46e+04    -  1.00e+00 6.26e-04h 11
     994  0.0000000e+00 2.82e+04 4.48e+17  -1.0 4.18e+04    -  1.00e+00 6.19e-04h 11
     995  0.0000000e+00 2.82e+04 4.48e+17  -1.0 4.17e+04    -  1.00e+00 6.18e-04h 11
     996  0.0000000e+00 2.82e+04 4.48e+17  -1.0 4.02e+04    -  1.00e+00 6.14e-04h 11
     997  0.0000000e+00 2.81e+04 4.48e+17  -1.0 4.26e+04    -  1.00e+00 6.21e-04h 11
     998  0.0000000e+00 2.81e+04 4.48e+17  -1.0 4.19e+04    -  1.00e+00 6.19e-04h 11
     999  0.0000000e+00 2.81e+04 4.48e+17  -1.0 5.04e+04    -  1.00e+00 6.42e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1000  0.0000000e+00 2.81e+04 4.48e+17  -1.0 4.22e+04    -  1.00e+00 6.19e-04h 11
    1001  0.0000000e+00 2.81e+04 4.48e+17  -1.0 4.48e+04    -  1.00e+00 6.27e-04h 11
    1002  0.0000000e+00 1.53e+05 1.15e+17  -1.0 4.09e+04    -  1.00e+00 6.31e-01w  1
    1003  0.0000000e+00 7.18e+04 1.77e+17  -1.0 9.51e+03    -  1.00e+00 6.04e-01w  1
    1004  0.0000000e+00 3.15e+05 9.08e+17  -1.0 1.39e+05    -  1.25e-01 1.71e-01w  1
    1005  0.0000000e+00 2.81e+04 4.49e+17  -1.0 2.07e+04    -  1.00e+00 6.16e-04h 10
    1006  0.0000000e+00 2.80e+04 4.49e+17  -1.0 4.18e+04    -  1.00e+00 6.19e-04h 11
    1007  0.0000000e+00 2.80e+04 4.49e+17  -1.0 4.00e+04    -  1.00e+00 6.14e-04h 11
    1008  0.0000000e+00 2.80e+04 4.49e+17  -1.0 4.09e+04    -  1.00e+00 6.16e-04h 11
    1009  0.0000000e+00 2.80e+04 4.49e+17  -1.0 4.00e+04    -  1.00e+00 6.14e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1010  0.0000000e+00 2.80e+04 4.49e+17  -1.0 4.10e+04    -  1.00e+00 6.17e-04h 11
    1011  0.0000000e+00 2.80e+04 4.49e+17  -1.0 3.97e+04    -  1.00e+00 6.13e-04h 11
    1012  0.0000000e+00 2.79e+04 4.49e+17  -1.0 4.04e+04    -  1.00e+00 6.15e-04h 11
    1013  0.0000000e+00 2.79e+04 4.49e+17  -1.0 3.99e+04    -  1.00e+00 6.14e-04h 11
    1014  0.0000000e+00 2.79e+04 4.49e+17  -1.0 3.99e+04    -  1.00e+00 6.14e-04h 11
    1015  0.0000000e+00 1.44e+05 1.12e+17  -1.0 3.97e+04    -  1.00e+00 6.28e-01w  1
    1016  0.0000000e+00 6.86e+04 1.79e+17  -1.0 9.73e+03    -  1.00e+00 6.04e-01w  1
    1017  0.0000000e+00 3.18e+05 7.28e+17  -1.0 1.49e+05    -  1.27e-01 1.60e-01w  1
    1018  0.0000000e+00 2.79e+04 4.49e+17  -1.0 1.93e+04    -  1.00e+00 6.13e-04h 10
    1019  0.0000000e+00 2.79e+04 4.50e+17  -1.0 3.99e+04    -  1.00e+00 6.14e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1020  0.0000000e+00 2.79e+04 4.50e+17  -1.0 3.98e+04    -  1.00e+00 6.14e-04h 11
    1021  0.0000000e+00 2.78e+04 4.50e+17  -1.0 3.96e+04    -  1.00e+00 6.13e-04h 11
    1022  0.0000000e+00 2.78e+04 4.50e+17  -1.0 3.96e+04    -  1.00e+00 6.13e-04h 11
    1023  0.0000000e+00 2.78e+04 4.50e+17  -1.0 3.97e+04    -  1.00e+00 6.13e-04h 11
    1024  0.0000000e+00 2.78e+04 4.50e+17  -1.0 3.93e+04    -  1.00e+00 6.12e-04h 11
    1025  0.0000000e+00 2.78e+04 4.50e+17  -1.0 4.15e+04    -  1.00e+00 6.18e-04h 11
    1026  0.0000000e+00 2.78e+04 4.50e+17  -1.0 3.94e+04    -  1.00e+00 6.13e-04h 11
    1027  0.0000000e+00 2.77e+04 4.50e+17  -1.0 4.00e+04    -  1.00e+00 6.14e-04h 11
    1028  0.0000000e+00 1.43e+05 1.11e+17  -1.0 3.92e+04    -  1.00e+00 6.27e-01w  1
    1029  0.0000000e+00 6.94e+04 1.83e+17  -1.0 1.01e+04    -  1.00e+00 6.03e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1030  0.0000000e+00 3.26e+05 3.77e+17  -1.0 1.71e+05    -  1.30e-01 1.41e-01w  1
    1031  0.0000000e+00 2.77e+04 4.50e+17  -1.0 1.76e+04    -  1.00e+00 6.12e-04h 10
    1032  0.0000000e+00 2.77e+04 4.51e+17  -1.0 3.94e+04    -  1.00e+00 6.13e-04h 11
    1033  0.0000000e+00 2.77e+04 4.51e+17  -1.0 3.91e+04    -  1.00e+00 6.12e-04h 11
    1034  0.0000000e+00 2.77e+04 4.51e+17  -1.0 4.09e+04    -  1.00e+00 6.17e-04h 11
    1035  0.0000000e+00 2.77e+04 4.51e+17  -1.0 3.92e+04    -  1.00e+00 6.12e-04h 11
    1036  0.0000000e+00 2.76e+04 4.51e+17  -1.0 3.91e+04    -  1.00e+00 6.12e-04h 11
    1037  0.0000000e+00 2.76e+04 4.51e+17  -1.0 3.93e+04    -  1.00e+00 6.12e-04h 11
    1038  0.0000000e+00 2.76e+04 4.51e+17  -1.0 3.93e+04    -  1.00e+00 6.12e-04h 11
    1039  0.0000000e+00 2.76e+04 4.51e+17  -1.0 3.94e+04    -  1.00e+00 6.12e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1040  0.0000000e+00 2.76e+04 4.51e+17  -1.0 3.89e+04    -  1.00e+00 6.11e-04h 11
    1041  0.0000000e+00 1.38e+05 1.17e+17  -1.0 4.25e+04    -  1.00e+00 6.36e-01w  1
    1042  0.0000000e+00 5.90e+04 1.60e+17  -1.0 5.46e+03    -  1.00e+00 6.26e-01w  1
    1043  0.0000000e+00 2.33e+05 3.81e+18  -1.0 6.51e+04    -  1.11e-01 3.12e-01w  1
    1044  0.0000000e+00 2.75e+04 4.52e+17  -1.0 3.92e+04    -  1.00e+00 6.21e-04h 10
    1045  0.0000000e+00 2.75e+04 4.52e+17  -1.0 3.92e+04    -  1.00e+00 6.12e-04h 11
    1046  0.0000000e+00 2.75e+04 4.52e+17  -1.0 3.90e+04    -  1.00e+00 6.11e-04h 11
    1047  0.0000000e+00 2.75e+04 4.52e+17  -1.0 3.89e+04    -  1.00e+00 6.11e-04h 11
    1048  0.0000000e+00 2.75e+04 4.52e+17  -1.0 3.88e+04    -  1.00e+00 6.11e-04h 11
    1049  0.0000000e+00 2.75e+04 4.52e+17  -1.0 3.91e+04    -  1.00e+00 6.12e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1050  0.0000000e+00 2.74e+04 4.52e+17  -1.0 3.91e+04    -  1.00e+00 6.12e-04h 11
    1051  0.0000000e+00 2.74e+04 4.52e+17  -1.0 4.23e+04    -  1.00e+00 6.20e-04h 11
    1052  0.0000000e+00 2.74e+04 4.52e+17  -1.0 4.05e+04    -  1.00e+00 6.15e-04h 11
    1053  0.0000000e+00 2.74e+04 4.52e+17  -1.0 3.93e+04    -  1.00e+00 6.12e-04h 11
    1054  0.0000000e+00 1.40e+05 1.10e+17  -1.0 3.91e+04    -  1.00e+00 6.27e-01w  1
    1055  0.0000000e+00 6.80e+04 1.88e+17  -1.0 1.03e+04    -  1.00e+00 6.04e-01w  1
    1056  0.0000000e+00 3.28e+05 2.75e+17  -1.0 1.81e+05    -  1.30e-01 1.34e-01w  1
    1057  0.0000000e+00 2.74e+04 4.53e+17  -1.0 1.65e+04    -  1.00e+00 6.12e-04h 10
    1058  0.0000000e+00 2.74e+04 4.53e+17  -1.0 3.90e+04    -  1.00e+00 6.12e-04h 11
    1059  0.0000000e+00 2.73e+04 4.53e+17  -1.0 3.97e+04    -  1.00e+00 6.13e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1060  0.0000000e+00 2.73e+04 4.53e+17  -1.0 3.90e+04    -  1.00e+00 6.12e-04h 11
    1061  0.0000000e+00 2.73e+04 4.53e+17  -1.0 3.90e+04    -  1.00e+00 6.12e-04h 11
    1062  0.0000000e+00 2.73e+04 4.53e+17  -1.0 3.89e+04    -  1.00e+00 6.11e-04h 11
    1063  0.0000000e+00 2.73e+04 4.53e+17  -1.0 3.89e+04    -  1.00e+00 6.11e-04h 11
    1064  0.0000000e+00 2.73e+04 4.53e+17  -1.0 4.19e+04    -  1.00e+00 6.19e-04h 11
    1065  0.0000000e+00 2.72e+04 4.53e+17  -1.0 4.02e+04    -  1.00e+00 6.15e-04h 11
    1066  0.0000000e+00 2.72e+04 4.54e+17  -1.0 3.91e+04    -  1.00e+00 6.12e-04h 11
    1067  0.0000000e+00 1.37e+05 1.13e+17  -1.0 4.06e+04    -  1.00e+00 6.31e-01w  1
    1068  0.0000000e+00 6.30e+04 1.77e+17  -1.0 7.94e+03    -  1.00e+00 6.16e-01w  1
    1069  0.0000000e+00 2.88e+05 1.87e+18  -1.0 1.08e+05    -  1.16e-01 2.11e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1070  0.0000000e+00 2.72e+04 4.54e+17  -1.0 2.58e+04    -  1.00e+00 6.16e-04h 10
    1071  0.0000000e+00 2.72e+04 4.54e+17  -1.0 3.89e+04    -  1.00e+00 6.11e-04h 11
    1072  0.0000000e+00 2.72e+04 4.54e+17  -1.0 3.87e+04    -  1.00e+00 6.11e-04h 11
    1073  0.0000000e+00 2.72e+04 4.54e+17  -1.0 3.84e+04    -  1.00e+00 6.10e-04h 11
    1074  0.0000000e+00 2.71e+04 4.54e+17  -1.0 3.87e+04    -  1.00e+00 6.11e-04h 11
    1075  0.0000000e+00 2.71e+04 4.54e+17  -1.0 3.86e+04    -  1.00e+00 6.11e-04h 11
    1076  0.0000000e+00 2.71e+04 4.54e+17  -1.0 3.97e+04    -  1.00e+00 6.14e-04h 11
    1077  0.0000000e+00 2.71e+04 4.54e+17  -1.0 3.85e+04    -  1.00e+00 6.11e-04h 11
    1078  0.0000000e+00 2.71e+04 4.55e+17  -1.0 3.85e+04    -  1.00e+00 6.10e-04h 11
    1079  0.0000000e+00 2.71e+04 4.55e+17  -1.0 3.85e+04    -  1.00e+00 6.10e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1080  0.0000000e+00 1.35e+05 1.08e+17  -1.0 3.85e+04    -  1.00e+00 6.25e-01w  1
    1081  0.0000000e+00 6.69e+04 1.93e+17  -1.0 1.05e+04    -  1.00e+00 6.05e-01w  1
    1082  0.0000000e+00 3.33e+05 3.89e+16  -1.0 2.09e+05    -  1.33e-01 1.16e-01w  1
    1083  0.0000000e+00 2.70e+04 4.55e+17  -1.0 1.51e+04    -  1.00e+00 6.10e-04h 10
    1084  0.0000000e+00 2.70e+04 4.55e+17  -1.0 3.81e+04    -  1.00e+00 6.10e-04h 11
    1085  0.0000000e+00 2.70e+04 4.55e+17  -1.0 3.94e+04    -  1.00e+00 6.13e-04h 11
    1086  0.0000000e+00 2.70e+04 4.55e+17  -1.0 3.85e+04    -  1.00e+00 6.10e-04h 11
    1087  0.0000000e+00 2.70e+04 4.55e+17  -1.0 3.84e+04    -  1.00e+00 6.10e-04h 11
    1088  0.0000000e+00 2.70e+04 4.55e+17  -1.0 3.95e+04    -  1.00e+00 6.13e-04h 11
    1089  0.0000000e+00 2.69e+04 4.55e+17  -1.0 3.85e+04    -  1.00e+00 6.10e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1090  0.0000000e+00 2.69e+04 4.56e+17  -1.0 3.90e+04    -  1.00e+00 6.12e-04h 11
    1091  0.0000000e+00 2.69e+04 4.56e+17  -1.0 3.90e+04    -  1.00e+00 6.12e-04h 11
    1092  0.0000000e+00 2.69e+04 4.56e+17  -1.0 3.80e+04    -  1.00e+00 6.09e-04h 11
    1093  0.0000000e+00 1.34e+05 1.06e+17  -1.0 3.79e+04    -  1.00e+00 6.24e-01w  1
    1094  0.0000000e+00 6.82e+04 1.95e+17  -1.0 1.07e+04    -  1.00e+00 6.05e-01w  1
    1095  0.0000000e+00 3.35e+05 3.38e+17  -1.0 2.40e+05    -  1.40e-01 1.02e-01w  1
    1096  0.0000000e+00 2.69e+04 4.56e+17  -1.0 1.54e+04    -  1.00e+00 6.09e-04h 10
    1097  0.0000000e+00 2.69e+04 4.56e+17  -1.0 3.94e+04    -  1.00e+00 6.13e-04h 11
    1098  0.0000000e+00 2.68e+04 4.56e+17  -1.0 3.83e+04    -  1.00e+00 6.10e-04h 11
    1099  0.0000000e+00 2.68e+04 4.56e+17  -1.0 3.79e+04    -  1.00e+00 6.09e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1100  0.0000000e+00 2.68e+04 4.56e+17  -1.0 3.93e+04    -  1.00e+00 6.13e-04h 11
    1101  0.0000000e+00 2.68e+04 4.57e+17  -1.0 3.79e+04    -  1.00e+00 6.09e-04h 11
    1102  0.0000000e+00 2.68e+04 4.57e+17  -1.0 3.92e+04    -  1.00e+00 6.12e-04h 11
    1103  0.0000000e+00 2.68e+04 4.57e+17  -1.0 3.87e+04    -  1.00e+00 6.11e-04h 11
    1104  0.0000000e+00 2.68e+04 4.57e+17  -1.0 3.80e+04    -  1.00e+00 6.09e-04h 11
    1105  0.0000000e+00 2.67e+04 4.57e+17  -1.0 3.80e+04    -  1.00e+00 6.09e-04h 11
    1106  0.0000000e+00 1.32e+05 1.06e+17  -1.0 3.79e+04    -  1.00e+00 6.24e-01w  1
    1107  0.0000000e+00 6.66e+04 1.99e+17  -1.0 1.10e+04    -  1.00e+00 6.04e-01w  1
    1108  0.0000000e+00 3.39e+05 4.44e+17  -1.0 2.55e+05    -  1.40e-01 9.62e-02w  1
    1109  0.0000000e+00 2.67e+04 4.57e+17  -1.0 1.34e+04    -  1.00e+00 6.09e-04h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1110  0.0000000e+00 2.67e+04 4.57e+17  -1.0 4.50e+04    -  1.00e+00 6.28e-04h 11
    1111  0.0000000e+00 2.67e+04 4.57e+17  -1.0 3.85e+04    -  1.00e+00 6.11e-04h 11
    1112  0.0000000e+00 2.67e+04 4.57e+17  -1.0 3.80e+04    -  1.00e+00 6.09e-04h 11
    1113  0.0000000e+00 2.67e+04 4.58e+17  -1.0 3.85e+04    -  1.00e+00 6.11e-04h 11
    1114  0.0000000e+00 2.66e+04 4.58e+17  -1.0 3.77e+04    -  1.00e+00 6.09e-04h 11
    1115  0.0000000e+00 2.66e+04 4.58e+17  -1.0 3.81e+04    -  1.00e+00 6.10e-04h 11
    1116  0.0000000e+00 2.66e+04 4.58e+17  -1.0 3.79e+04    -  1.00e+00 6.09e-04h 11
    1117  0.0000000e+00 2.66e+04 4.58e+17  -1.0 3.80e+04    -  1.00e+00 6.10e-04h 11
    1118  0.0000000e+00 2.66e+04 4.58e+17  -1.0 3.85e+04    -  1.00e+00 6.11e-04h 11
    1119  0.0000000e+00 1.29e+05 1.08e+17  -1.0 3.92e+04    -  1.00e+00 6.27e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1120  0.0000000e+00 6.25e+04 1.89e+17  -1.0 8.70e+03    -  1.00e+00 6.16e-01w  1
    1121  0.0000000e+00 3.06e+05 1.06e+18  -1.0 1.42e+05    -  1.19e-01 1.64e-01w  1
    1122  0.0000000e+00 2.66e+04 4.58e+17  -1.0 2.02e+04    -  1.00e+00 6.13e-04h 10
    1123  0.0000000e+00 2.65e+04 4.58e+17  -1.0 3.79e+04    -  1.00e+00 6.09e-04h 11
    1124  0.0000000e+00 2.65e+04 4.59e+17  -1.0 3.76e+04    -  1.00e+00 6.09e-04h 11
    1125  0.0000000e+00 2.65e+04 4.59e+17  -1.0 3.80e+04    -  1.00e+00 6.09e-04h 11
    1126  0.0000000e+00 2.65e+04 4.59e+17  -1.0 3.74e+04    -  1.00e+00 6.08e-04h 11
    1127  0.0000000e+00 2.65e+04 4.59e+17  -1.0 3.96e+04    -  1.00e+00 6.14e-04h 11
    1128  0.0000000e+00 2.65e+04 4.59e+17  -1.0 4.02e+04    -  1.00e+00 6.15e-04h 11
    1129  0.0000000e+00 2.64e+04 4.59e+17  -1.0 3.80e+04    -  1.00e+00 6.09e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1130  0.0000000e+00 2.64e+04 4.59e+17  -1.0 3.87e+04    -  1.00e+00 6.11e-04h 11
    1131  0.0000000e+00 2.64e+04 4.59e+17  -1.0 3.74e+04    -  1.00e+00 6.08e-04h 11
    1132  0.0000000e+00 1.29e+05 1.04e+17  -1.0 3.73e+04    -  1.00e+00 6.22e-01w  1
    1133  0.0000000e+00 6.72e+04 2.07e+17  -1.0 1.15e+04    -  1.00e+00 6.03e-01w  1
    1134  0.0000000e+00 3.46e+05 1.00e+18  -1.0 3.48e+05    -  1.59e-01 7.06e-02w  1
    1135  0.0000000e+00 2.64e+04 4.60e+17  -1.0 1.21e+04    -  1.00e+00 6.08e-04h 10
    1136  0.0000000e+00 2.64e+04 4.60e+17  -1.0 3.72e+04    -  1.00e+00 6.08e-04h 11
    1137  0.0000000e+00 2.64e+04 4.60e+17  -1.0 3.80e+04    -  1.00e+00 6.10e-04h 11
    1138  0.0000000e+00 2.63e+04 4.60e+17  -1.0 3.74e+04    -  1.00e+00 6.08e-04h 11
    1139  0.0000000e+00 2.63e+04 4.60e+17  -1.0 3.72e+04    -  1.00e+00 6.08e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1140  0.0000000e+00 2.63e+04 4.60e+17  -1.0 3.70e+04    -  1.00e+00 6.07e-04h 11
    1141  0.0000000e+00 2.63e+04 4.60e+17  -1.0 3.74e+04    -  1.00e+00 6.08e-04h 11
    1142  0.0000000e+00 2.63e+04 4.60e+17  -1.0 3.87e+04    -  1.00e+00 6.12e-04h 11
    1143  0.0000000e+00 2.63e+04 4.61e+17  -1.0 3.71e+04    -  1.00e+00 6.07e-04h 11
    1144  0.0000000e+00 2.62e+04 4.61e+17  -1.0 3.72e+04    -  1.00e+00 6.08e-04h 11
    1145  0.0000000e+00 1.27e+05 1.03e+17  -1.0 3.69e+04    -  1.00e+00 6.22e-01w  1
    1146  0.0000000e+00 6.67e+04 2.11e+17  -1.0 1.18e+04    -  1.00e+00 6.03e-01w  1
    1147  0.0000000e+00 3.47e+05 1.27e+18  -1.0 4.16e+05    -  1.78e-01 5.89e-02w  1
    1148  0.0000000e+00 2.62e+04 4.61e+17  -1.0 1.19e+04    -  1.00e+00 6.07e-04h 10
    1149  0.0000000e+00 2.62e+04 4.61e+17  -1.0 3.70e+04    -  1.00e+00 6.07e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1150  0.0000000e+00 2.62e+04 4.61e+17  -1.0 3.77e+04    -  1.00e+00 6.09e-04h 11
    1151  0.0000000e+00 2.62e+04 4.61e+17  -1.0 3.77e+04    -  1.00e+00 6.09e-04h 11
    1152  0.0000000e+00 2.62e+04 4.61e+17  -1.0 4.03e+04    -  1.00e+00 6.16e-04h 11
    1153  0.0000000e+00 2.62e+04 4.61e+17  -1.0 3.75e+04    -  1.00e+00 6.08e-04h 11
    1154  0.0000000e+00 2.61e+04 4.62e+17  -1.0 3.77e+04    -  1.00e+00 6.09e-04h 11
    1155  0.0000000e+00 2.61e+04 4.62e+17  -1.0 3.77e+04    -  1.00e+00 6.09e-04h 11
    1156  0.0000000e+00 2.61e+04 4.62e+17  -1.0 3.71e+04    -  1.00e+00 6.07e-04h 11
    1157  0.0000000e+00 2.61e+04 4.62e+17  -1.0 3.70e+04    -  1.00e+00 6.07e-04h 11
    1158  0.0000000e+00 1.25e+05 1.03e+17  -1.0 3.67e+04    -  1.00e+00 6.21e-01w  1
    1159  0.0000000e+00 6.66e+04 2.14e+17  -1.0 1.20e+04    -  1.00e+00 6.03e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1160  0.0000000e+00 3.44e+05 1.49e+18  -1.0 4.96e+05    -  2.05e-01 4.91e-02w  1
    1161  0.0000000e+00 2.61e+04 4.62e+17  -1.0 1.17e+04    -  1.00e+00 6.07e-04h 10
    1162  0.0000000e+00 2.61e+04 4.62e+17  -1.0 3.84e+04    -  1.00e+00 6.11e-04h 11
    1163  0.0000000e+00 2.60e+04 4.62e+17  -1.0 3.83e+04    -  1.00e+00 6.11e-04h 11
    1164  0.0000000e+00 2.60e+04 4.63e+17  -1.0 3.80e+04    -  1.00e+00 6.10e-04h 11
    1165  0.0000000e+00 2.60e+04 4.63e+17  -1.0 3.67e+04    -  1.00e+00 6.07e-04h 11
    1166  0.0000000e+00 2.60e+04 4.63e+17  -1.0 3.68e+04    -  1.00e+00 6.07e-04h 11
    1167  0.0000000e+00 2.60e+04 4.63e+17  -1.0 3.68e+04    -  1.00e+00 6.07e-04h 11
    1168  0.0000000e+00 2.60e+04 4.63e+17  -1.0 3.83e+04    -  1.00e+00 6.11e-04h 11
    1169  0.0000000e+00 2.59e+04 4.63e+17  -1.0 3.67e+04    -  1.00e+00 6.07e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1170  0.0000000e+00 2.59e+04 4.63e+17  -1.0 3.77e+04    -  1.00e+00 6.09e-04h 11
    1171  0.0000000e+00 1.23e+05 1.02e+17  -1.0 3.67e+04    -  1.00e+00 6.21e-01w  1
    1172  0.0000000e+00 6.59e+04 2.16e+17  -1.0 1.18e+04    -  1.00e+00 6.05e-01w  1
    1173  0.0000000e+00 3.39e+05 1.05e+18  -1.0 3.27e+05    -  1.56e-01 6.82e-02w  1
    1174  0.0000000e+00 2.59e+04 4.63e+17  -1.0 1.38e+04    -  1.00e+00 6.07e-04h 10
    1175  0.0000000e+00 2.59e+04 4.64e+17  -1.0 3.64e+04    -  1.00e+00 6.06e-04h 11
    1176  0.0000000e+00 2.59e+04 4.64e+17  -1.0 3.72e+04    -  1.00e+00 6.08e-04h 11
    1177  0.0000000e+00 2.59e+04 4.64e+17  -1.0 3.65e+04    -  1.00e+00 6.06e-04h 11
    1178  0.0000000e+00 2.59e+04 4.64e+17  -1.0 3.75e+04    -  1.00e+00 6.09e-04h 11
    1179  0.0000000e+00 2.58e+04 4.64e+17  -1.0 3.63e+04    -  1.00e+00 6.06e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1180  0.0000000e+00 2.58e+04 4.64e+17  -1.0 3.72e+04    -  1.00e+00 6.08e-04h 11
    1181  0.0000000e+00 2.58e+04 4.64e+17  -1.0 3.64e+04    -  1.00e+00 6.06e-04h 11
    1182  0.0000000e+00 2.58e+04 4.65e+17  -1.0 3.64e+04    -  1.00e+00 6.06e-04h 11
    1183  0.0000000e+00 2.58e+04 4.65e+17  -1.0 3.73e+04    -  1.00e+00 6.08e-04h 11
    1184  0.0000000e+00 1.21e+05 1.01e+17  -1.0 3.63e+04    -  1.00e+00 6.20e-01w  1
    1185  0.0000000e+00 6.58e+04 2.17e+17  -1.0 1.17e+04    -  1.00e+00 6.06e-01w  1
    1186  0.0000000e+00 3.33e+05 1.57e+18  -1.0 5.19e+05    -  2.17e-01 4.61e-02w  1
    1187  0.0000000e+00 2.58e+04 4.65e+17  -1.0 1.21e+04    -  1.00e+00 6.06e-04h 10
    1188  0.0000000e+00 2.57e+04 4.65e+17  -1.0 3.63e+04    -  1.00e+00 6.06e-04h 11
    1189  0.0000000e+00 2.57e+04 4.65e+17  -1.0 3.77e+04    -  1.00e+00 6.09e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1190  0.0000000e+00 2.57e+04 4.65e+17  -1.0 3.63e+04    -  1.00e+00 6.06e-04h 11
    1191  0.0000000e+00 2.57e+04 4.65e+17  -1.0 3.73e+04    -  1.00e+00 6.08e-04h 11
    1192  0.0000000e+00 2.57e+04 4.66e+17  -1.0 3.60e+04    -  1.00e+00 6.05e-04h 11
    1193  0.0000000e+00 2.57e+04 4.66e+17  -1.0 3.62e+04    -  1.00e+00 6.06e-04h 11
    1194  0.0000000e+00 2.57e+04 4.66e+17  -1.0 3.59e+04    -  1.00e+00 6.05e-04h 11
    1195  0.0000000e+00 2.56e+04 4.66e+17  -1.0 3.61e+04    -  1.00e+00 6.06e-04h 11
    1196  0.0000000e+00 2.56e+04 4.66e+17  -1.0 3.61e+04    -  1.00e+00 6.05e-04h 11
    1197  0.0000000e+00 1.19e+05 1.00e+17  -1.0 3.59e+04    -  1.00e+00 6.19e-01w  1
    1198  0.0000000e+00 6.53e+04 2.23e+17  -1.0 1.23e+04    -  1.00e+00 6.04e-01w  1
    1199  0.0000000e+00 3.32e+05 3.93e+17  -1.0 9.33e+05    -  1.33e-01 2.55e-02w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1200  0.0000000e+00 2.56e+04 4.66e+17  -1.0 1.19e+04    -  1.00e+00 6.05e-04h 10
    1201  0.0000000e+00 2.56e+04 4.66e+17  -1.0 3.69e+04    -  1.00e+00 6.08e-04h 11
    1202  0.0000000e+00 2.56e+04 4.67e+17  -1.0 3.58e+04    -  1.00e+00 6.05e-04h 11
    1203  0.0000000e+00 2.56e+04 4.67e+17  -1.0 3.63e+04    -  1.00e+00 6.06e-04h 11
    1204  0.0000000e+00 2.55e+04 4.67e+17  -1.0 3.76e+04    -  1.00e+00 6.09e-04h 11
    1205  0.0000000e+00 2.55e+04 4.67e+17  -1.0 3.57e+04    -  1.00e+00 6.05e-04h 11
    1206  0.0000000e+00 2.55e+04 4.67e+17  -1.0 3.86e+04    -  1.00e+00 6.12e-04h 11
    1207  0.0000000e+00 2.55e+04 4.67e+17  -1.0 3.63e+04    -  1.00e+00 6.06e-04h 11
    1208  0.0000000e+00 2.55e+04 4.67e+17  -1.0 3.63e+04    -  1.00e+00 6.06e-04h 11
    1209  0.0000000e+00 2.55e+04 4.68e+17  -1.0 3.59e+04    -  1.00e+00 6.05e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1210  0.0000000e+00 1.17e+05 1.00e+17  -1.0 3.59e+04    -  1.00e+00 6.20e-01w  1
    1211  0.0000000e+00 6.46e+04 2.24e+17  -1.0 1.20e+04    -  1.00e+00 6.07e-01w  1
    1212  0.0000000e+00 3.37e+05 1.36e+18  -1.0 3.95e+05    -  1.75e-01 5.49e-02w  1
    1213  0.0000000e+00 2.54e+04 4.68e+17  -1.0 1.41e+04    -  1.00e+00 6.05e-04h 10
    1214  0.0000000e+00 2.54e+04 4.68e+17  -1.0 3.59e+04    -  1.00e+00 6.05e-04h 11
    1215  0.0000000e+00 2.54e+04 4.68e+17  -1.0 3.74e+04    -  1.00e+00 6.09e-04h 11
    1216  0.0000000e+00 2.54e+04 4.68e+17  -1.0 3.55e+04    -  1.00e+00 6.04e-04h 11
    1217  0.0000000e+00 2.54e+04 4.68e+17  -1.0 3.89e+04    -  1.00e+00 1.23e-03h 10
    1218  0.0000000e+00 2.54e+04 4.69e+17  -1.0 3.75e+04    -  1.00e+00 6.09e-04h 11
    1219  0.0000000e+00 2.53e+04 4.69e+17  -1.0 3.54e+04    -  1.00e+00 6.04e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1220  0.0000000e+00 2.53e+04 4.69e+17  -1.0 3.61e+04    -  1.00e+00 6.06e-04h 11
    1221  0.0000000e+00 2.53e+04 4.69e+17  -1.0 3.54e+04    -  1.00e+00 6.04e-04h 11
    1222  0.0000000e+00 2.53e+04 4.69e+17  -1.0 3.53e+04    -  1.00e+00 6.04e-04h 11
    1223  0.0000000e+00 1.14e+05 9.87e+16  -1.0 3.54e+04    -  1.00e+00 6.19e-01w  1
    1224  0.0000000e+00 6.42e+04 2.29e+17  -1.0 1.24e+04    -  1.00e+00 6.06e-01w  1
    1225  0.0000000e+00 3.27e+05 2.01e+17  -1.0 1.39e+06    -  8.78e-02 1.67e-02w  1
    1226  0.0000000e+00 2.53e+04 4.69e+17  -1.0 1.23e+04    -  1.00e+00 6.04e-04h 10
    1227  0.0000000e+00 2.53e+04 4.69e+17  -1.0 3.51e+04    -  1.00e+00 6.03e-04h 11
    1228  0.0000000e+00 2.52e+04 4.70e+17  -1.0 3.59e+04    -  1.00e+00 6.05e-04h 11
    1229  0.0000000e+00 2.52e+04 4.70e+17  -1.0 3.50e+04    -  1.00e+00 6.03e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1230  0.0000000e+00 2.52e+04 4.70e+17  -1.0 3.51e+04    -  1.00e+00 6.04e-04h 11
    1231  0.0000000e+00 2.52e+04 4.70e+17  -1.0 3.51e+04    -  1.00e+00 6.04e-04h 11
    1232  0.0000000e+00 2.52e+04 4.70e+17  -1.0 3.49e+04    -  1.00e+00 6.03e-04h 11
    1233  0.0000000e+00 2.52e+04 4.70e+17  -1.0 3.56e+04    -  1.00e+00 6.05e-04h 11
    1234  0.0000000e+00 2.52e+04 4.70e+17  -1.0 3.51e+04    -  1.00e+00 6.04e-04h 11
    1235  0.0000000e+00 2.51e+04 4.71e+17  -1.0 3.50e+04    -  1.00e+00 6.03e-04h 11
    1236  0.0000000e+00 1.10e+05 9.87e+16  -1.0 3.55e+04    -  1.00e+00 6.19e-01w  1
    1237  0.0000000e+00 6.16e+04 2.25e+17  -1.0 1.15e+04    -  1.00e+00 6.12e-01w  1
    1238  0.0000000e+00 2.69e+05 1.66e+18  -1.0 4.85e+04    -  1.07e-01 1.87e-01w  1
    1239  0.0000000e+00 2.51e+04 4.71e+17  -1.0 3.79e+04    -  1.00e+00 1.21e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1240  0.0000000e+00 2.51e+04 4.71e+17  -1.0 3.46e+04    -  1.00e+00 1.21e-03h 10
    1241  0.0000000e+00 2.51e+04 4.72e+17  -1.0 3.44e+04    -  1.00e+00 1.20e-03h 10
    1242  0.0000000e+00 2.50e+04 4.72e+17  -1.0 3.44e+04    -  1.00e+00 1.20e-03h 10
    1243  0.0000000e+00 2.50e+04 4.72e+17  -1.0 3.49e+04    -  1.00e+00 1.21e-03h 10
    1244  0.0000000e+00 2.50e+04 4.72e+17  -1.0 3.47e+04    -  1.00e+00 1.21e-03h 10
    1245  0.0000000e+00 2.49e+04 4.73e+17  -1.0 3.41e+04    -  1.00e+00 1.20e-03h 10
    1246  0.0000000e+00 2.49e+04 4.73e+17  -1.0 3.40e+04    -  1.00e+00 1.20e-03h 10
    1247  0.0000000e+00 2.49e+04 4.73e+17  -1.0 3.36e+04    -  1.00e+00 1.20e-03h 10
    1248  0.0000000e+00 2.48e+04 4.74e+17  -1.0 3.42e+04    -  1.00e+00 1.20e-03h 10
    1249  0.0000000e+00 1.01e+05 9.38e+16  -1.0 3.34e+04    -  1.00e+00 6.15e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1250  0.0000000e+00 6.05e+04 2.39e+17  -1.0 1.32e+04    -  1.00e+00 6.07e-01w  1
    WARNING: Problem in step computation; switching to emergency mode.
    1251  0.0000000e+00 2.48e+04 4.74e+17  -1.0 1.32e+04  20.0 1.00e+00 1.20e-03h 10
    1252  0.0000000e+00 2.48e+04 4.74e+17  -1.0 3.33e+04    -  1.00e+00 1.20e-03h 10
    1253  0.0000000e+00 2.48e+04 4.75e+17  -1.0 3.35e+04    -  1.00e+00 1.20e-03h 10
    1254  0.0000000e+00 2.47e+04 4.75e+17  -1.0 3.50e+04    -  1.00e+00 1.21e-03h 10
    1255  0.0000000e+00 2.47e+04 4.75e+17  -1.0 3.32e+04    -  1.00e+00 1.20e-03h 10
    1256  0.0000000e+00 2.47e+04 4.76e+17  -1.0 3.33e+04    -  1.00e+00 1.20e-03h 10
    1257  0.0000000e+00 2.46e+04 4.76e+17  -1.0 3.29e+04    -  1.00e+00 1.20e-03h 10
    1258  0.0000000e+00 2.46e+04 4.76e+17  -1.0 3.25e+04    -  1.00e+00 1.20e-03h 10
    1259  0.0000000e+00 2.46e+04 4.77e+17  -1.0 3.51e+04    -  1.00e+00 1.21e-03h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1260  0.0000000e+00 2.45e+04 4.77e+17  -1.0 3.23e+04    -  1.00e+00 1.20e-03h 10
    1261  0.0000000e+00 9.07e+04 9.04e+16  -1.0 3.21e+04    -  1.00e+00 6.12e-01w  1
    1262  0.0000000e+00 5.78e+04 2.45e+17  -1.0 1.35e+04    -  1.00e+00 6.08e-01w  1
    WARNING: Problem in step computation; switching to emergency mode.
    1263  0.0000000e+00 2.45e+04 4.77e+17  -1.0 1.35e+04    -  1.00e+00 1.20e-03h 10
    1264  0.0000000e+00 2.45e+04 4.78e+17  -1.0 3.22e+04    -  1.00e+00 1.20e-03h 10
    1265  0.0000000e+00 2.45e+04 4.78e+17  -1.0 3.23e+04    -  1.00e+00 1.20e-03h 10
    1266  0.0000000e+00 2.44e+04 4.78e+17  -1.0 3.28e+04    -  1.00e+00 1.20e-03h 10
    1267  0.0000000e+00 2.44e+04 4.79e+17  -1.0 3.14e+04    -  1.00e+00 1.19e-03h 10
    1268  0.0000000e+00 2.44e+04 4.79e+17  -1.0 3.14e+04    -  1.00e+00 1.19e-03h 10
    1269  0.0000000e+00 2.43e+04 4.79e+17  -1.0 3.12e+04    -  1.00e+00 1.19e-03h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1270  0.0000000e+00 2.43e+04 4.80e+17  -1.0 3.28e+04    -  1.00e+00 1.20e-03h 10
    1271  0.0000000e+00 2.43e+04 4.80e+17  -1.0 3.13e+04    -  1.00e+00 1.19e-03h 10
    1272  0.0000000e+00 2.43e+04 4.80e+17  -1.0 3.06e+04    -  1.00e+00 1.19e-03h 10
    1273  0.0000000e+00 7.77e+04 8.78e+16  -1.0 3.11e+04    -  1.00e+00 6.11e-01w  1
    1274  0.0000000e+00 5.14e+04 2.42e+17  -1.0 1.26e+04    -  1.00e+00 6.16e-01w  1
    WARNING: Problem in step computation; switching to emergency mode.
    1275  0.0000000e+00 2.42e+04 4.81e+17  -1.0 1.26e+04  20.0 1.00e+00 1.19e-03h 10
    1276  0.0000000e+00 2.42e+04 4.81e+17  -1.0 3.14e+04    -  1.00e+00 1.19e-03h 10
    1277  0.0000000e+00 2.42e+04 4.81e+17  -1.0 3.13e+04    -  1.00e+00 1.19e-03h 10
    1278  0.0000000e+00 2.41e+04 4.82e+17  -1.0 3.01e+04    -  1.00e+00 1.19e-03h 10
    1279  0.0000000e+00 2.41e+04 4.82e+17  -1.0 2.97e+04    -  1.00e+00 1.19e-03h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1280  0.0000000e+00 2.41e+04 4.82e+17  -1.0 2.97e+04    -  1.00e+00 1.19e-03h 10
    1281  0.0000000e+00 2.41e+04 4.83e+17  -1.0 2.93e+04    -  1.00e+00 1.18e-03h 10
    1282  0.0000000e+00 2.40e+04 4.83e+17  -1.0 2.97e+04    -  1.00e+00 1.19e-03h 10
    1283  0.0000000e+00 2.40e+04 4.83e+17  -1.0 2.91e+04    -  1.00e+00 1.18e-03h 10
    1284  0.0000000e+00 2.40e+04 4.84e+17  -1.0 2.88e+04    -  1.00e+00 1.18e-03h 10
    1285  0.0000000e+00 6.56e+04 8.22e+16  -1.0 2.87e+04    -  1.00e+00 6.06e-01w  1
    1286  0.0000000e+00 4.94e+04 2.53e+17  -1.0 1.40e+04    -  1.00e+00 6.14e-01w  1
    WARNING: Problem in step computation; switching to emergency mode.
    1287  0.0000000e+00 2.39e+04 4.84e+17  -1.0 1.40e+04  20.0 1.00e+00 1.18e-03h 10
    1288  0.0000000e+00 2.39e+04 4.84e+17  -1.0 2.84e+04    -  1.00e+00 1.18e-03h 10
    1289  0.0000000e+00 2.39e+04 4.85e+17  -1.0 2.95e+04    -  1.00e+00 1.19e-03h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1290  0.0000000e+00 2.39e+04 4.85e+17  -1.0 2.95e+04    -  1.00e+00 1.19e-03h 10
    1291  0.0000000e+00 2.38e+04 4.85e+17  -1.0 2.79e+04    -  1.00e+00 1.18e-03h 10
    1292  0.0000000e+00 2.38e+04 4.86e+17  -1.0 2.77e+04    -  1.00e+00 1.18e-03h 10
    1293  0.0000000e+00 2.38e+04 4.86e+17  -1.0 2.74e+04    -  1.00e+00 1.18e-03h 10
    1294  0.0000000e+00 2.37e+04 4.86e+17  -1.0 2.77e+04    -  1.00e+00 1.18e-03h 10
    1295  0.0000000e+00 2.37e+04 4.87e+17  -1.0 2.69e+04    -  1.00e+00 1.18e-03h 10
    1296  0.0000000e+00 2.37e+04 4.87e+17  -1.0 2.70e+04    -  1.00e+00 1.18e-03h 10
    1297  0.0000000e+00 5.38e+04 7.80e+16  -1.0 2.69e+04    -  1.00e+00 6.03e-01w  1
    1298  0.0000000e+00 4.57e+04 2.57e+17  -1.0 1.42e+04    -  1.00e+00 6.19e-01w  1
    WARNING: Problem in step computation; switching to emergency mode.
    1299  0.0000000e+00 2.37e+04 4.88e+17  -1.0 1.42e+04    -  1.00e+00 1.18e-03h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1300  0.0000000e+00 2.36e+04 4.88e+17  -1.0 2.77e+04    -  1.00e+00 2.36e-03h  9
    1301  0.0000000e+00 2.35e+04 4.89e+17  -1.0 2.61e+04    -  1.00e+00 2.35e-03h  9
    1302  0.0000000e+00 2.35e+04 4.90e+17  -1.0 2.62e+04    -  1.00e+00 2.35e-03h  9
    1303  0.0000000e+00 2.34e+04 4.90e+17  -1.0 2.64e+04    -  1.00e+00 2.35e-03h  9
    1304  0.0000000e+00 2.34e+04 4.91e+17  -1.0 2.50e+04    -  1.00e+00 2.34e-03h  9
    1305  0.0000000e+00 2.33e+04 4.92e+17  -1.0 2.53e+04    -  1.00e+00 2.34e-03h  9
    1306  0.0000000e+00 2.33e+04 4.93e+17  -1.0 2.43e+04    -  1.00e+00 2.34e-03h  9
    1307  0.0000000e+00 2.32e+04 4.93e+17  -1.0 2.43e+04    -  1.00e+00 2.34e-03h  9
    1308  0.0000000e+00 2.32e+04 4.94e+17  -1.0 2.35e+04    -  1.00e+00 2.33e-03h  9
    1309  0.0000000e+00 3.65e+04 7.19e+16  -1.0 2.44e+04    -  1.00e+00 5.98e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1310  0.0000000e+00 3.91e+04 2.60e+17  -1.0 1.37e+04    -  1.00e+00 6.30e-01w  1
    WARNING: Problem in step computation; switching to emergency mode.
    1311  0.0000000e+00 2.31e+04 4.95e+17  -1.0 1.37e+04  20.0 1.00e+00 2.34e-03h  9
    1312  0.0000000e+00 2.31e+04 4.96e+17  -1.0 2.35e+04    -  1.00e+00 2.33e-03h  9
    1313  0.0000000e+00 2.30e+04 4.97e+17  -1.0 2.33e+04    -  1.00e+00 2.33e-03h  9
    1314  0.0000000e+00 2.29e+04 4.97e+17  -1.0 2.22e+04    -  1.00e+00 2.32e-03h  9
    1315  0.0000000e+00 2.29e+04 4.98e+17  -1.0 2.29e+04    -  1.00e+00 2.33e-03h  9
    1316  0.0000000e+00 2.28e+04 4.99e+17  -1.0 2.16e+04    -  1.00e+00 2.32e-03h  9
    1317  0.0000000e+00 2.28e+04 5.00e+17  -1.0 2.19e+04    -  1.00e+00 2.32e-03h  9
    1318  0.0000000e+00 2.27e+04 5.01e+17  -1.0 2.13e+04    -  1.00e+00 2.32e-03h  9
    1319  0.0000000e+00 2.27e+04 5.01e+17  -1.0 2.13e+04    -  1.00e+00 2.32e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1320  0.0000000e+00 2.26e+04 5.02e+17  -1.0 2.05e+04    -  1.00e+00 2.31e-03h  9
    1321  0.0000000e+00 2.44e+04 6.46e+16  -1.0 2.11e+04    -  1.00e+00 5.93e-01w  1
    1322  0.0000000e+00 3.94e+04 2.75e+17  -1.0 1.51e+04    -  1.00e+00 6.34e-01w  1
    WARNING: Problem in step computation; switching to emergency mode.
    1323  0.0000000e+00 2.25e+04 5.04e+17  -1.0 1.51e+04  20.0 1.00e+00 4.63e-03h  8
    1324  0.0000000e+00 2.24e+04 5.06e+17  -1.0 1.96e+04    -  1.00e+00 4.61e-03h  8
    1325  0.0000000e+00 2.23e+04 5.07e+17  -1.0 1.91e+04    -  1.00e+00 4.60e-03h  8
    1326  0.0000000e+00 2.22e+04 5.09e+17  -1.0 1.88e+04    -  1.00e+00 4.60e-03h  8
    1327  0.0000000e+00 2.21e+04 5.11e+17  -1.0 1.90e+04    -  1.00e+00 4.61e-03h  8
    1328  0.0000000e+00 2.20e+04 5.12e+17  -1.0 1.91e+04    -  1.00e+00 4.61e-03h  8
    1329  0.0000000e+00 2.19e+04 5.14e+17  -1.0 1.76e+04    -  1.00e+00 4.59e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1330  0.0000000e+00 2.18e+04 5.16e+17  -1.0 1.69e+04    -  1.00e+00 4.58e-03h  8
    1331  0.0000000e+00 2.17e+04 5.18e+17  -1.0 1.68e+04    -  1.00e+00 4.58e-03h  8
    1332  0.0000000e+00 2.16e+04 5.20e+17  -1.0 1.64e+04    -  1.00e+00 4.58e-03h  8
    1333  0.0000000e+00 1.25e+04 5.42e+16  -1.0 1.64e+04    -  1.00e+00 5.86e-01w  1
    1334  0.0000000e+00 4.94e+04 3.16e+17  -1.0 1.84e+04    -  1.00e+00 6.40e-01w  1
    WARNING: Problem in step computation; switching to emergency mode.
    1335  0.0000000e+00 2.14e+04 5.24e+17  -1.0 1.84e+04    -  1.00e+00 9.16e-03h  7
    1336  0.0000000e+00 2.12e+04 5.28e+17  -1.0 1.59e+04    -  1.00e+00 9.15e-03h  7
    1337  0.0000000e+00 2.10e+04 5.32e+17  -1.0 1.64e+04    -  1.00e+00 9.17e-03h  7
    1338  0.0000000e+00 2.09e+04 5.34e+17  -1.0 1.63e+04    -  1.00e+00 4.59e-03h  8
    1339  0.0000000e+00 2.08e+04 5.36e+17  -1.0 1.67e+04    -  1.00e+00 4.59e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1340  0.0000000e+00 2.07e+04 5.38e+17  -1.0 1.64e+04    -  1.00e+00 4.59e-03h  8
    1341  0.0000000e+00 2.06e+04 5.40e+17  -1.0 1.66e+04    -  1.00e+00 4.60e-03h  8
    1342  0.0000000e+00 2.05e+04 5.42e+17  -1.0 1.64e+04    -  1.00e+00 4.59e-03h  8
    1343  0.0000000e+00 2.04e+04 5.44e+17  -1.0 1.71e+04    -  1.00e+00 4.61e-03h  8
    1344  0.0000000e+00 2.04e+04 5.47e+17  -1.0 1.70e+04    -  1.00e+00 4.61e-03h  8
    1345  0.0000000e+00 1.53e+04 5.46e+16  -1.0 1.71e+04    -  1.00e+00 5.90e-01w  1
    1346  0.0000000e+00 6.99e+04 3.82e+17  -1.0 2.11e+04    -  1.00e+00 6.39e-01w  1
    WARNING: Problem in step computation; switching to emergency mode.
    1347  0.0000000e+00 2.03e+04 5.49e+17  -1.0 2.11e+04  20.0 1.00e+00 4.61e-03h  8
    1348  0.0000000e+00 2.02e+04 5.51e+17  -1.0 1.69e+04    -  1.00e+00 4.61e-03h  8
    1349  0.0000000e+00 2.01e+04 5.53e+17  -1.0 1.69e+04    -  1.00e+00 4.61e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1350  0.0000000e+00 2.00e+04 5.56e+17  -1.0 1.66e+04    -  1.00e+00 4.61e-03h  8
    1351  0.0000000e+00 1.99e+04 5.58e+17  -1.0 1.77e+04    -  1.00e+00 4.63e-03h  8
    1352  0.0000000e+00 1.98e+04 5.60e+17  -1.0 1.71e+04    -  1.00e+00 4.62e-03h  8
    1353  0.0000000e+00 1.97e+04 5.63e+17  -1.0 1.73e+04    -  1.00e+00 4.62e-03h  8
    1354  0.0000000e+00 1.96e+04 5.65e+17  -1.0 1.72e+04    -  1.00e+00 4.62e-03h  8
    1355  0.0000000e+00 1.95e+04 5.68e+17  -1.0 1.69e+04    -  1.00e+00 4.62e-03h  8
    1356  0.0000000e+00 1.94e+04 5.70e+17  -1.0 1.69e+04    -  1.00e+00 4.62e-03h  8
    1357  0.0000000e+00 1.72e+04 5.37e+16  -1.0 1.68e+04    -  1.00e+00 5.92e-01w  1
    1358  0.0000000e+00 8.86e+04 3.63e+17  -1.0 2.44e+04    -  1.00e+00 6.17e-01w  1
    WARNING: Problem in step computation; switching to emergency mode.
    1359  0.0000000e+00 1.93e+04 5.73e+17  -1.0 2.44e+04    -  1.00e+00 4.62e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1360  0.0000000e+00 1.93e+04 5.75e+17  -1.0 1.67e+04    -  1.00e+00 4.62e-03h  8
    1361  0.0000000e+00 1.92e+04 5.78e+17  -1.0 1.68e+04    -  1.00e+00 4.62e-03h  8
    1362  0.0000000e+00 1.91e+04 5.80e+17  -1.0 1.70e+04    -  1.00e+00 4.63e-03h  8
    1363  0.0000000e+00 1.90e+04 5.83e+17  -1.0 1.69e+04    -  1.00e+00 4.63e-03h  8
    1364  0.0000000e+00 1.89e+04 5.86e+17  -1.0 1.69e+04    -  1.00e+00 4.63e-03h  8
    1365  0.0000000e+00 1.88e+04 5.88e+17  -1.0 1.68e+04    -  1.00e+00 4.63e-03h  8
    1366  0.0000000e+00 1.87e+04 5.91e+17  -1.0 1.68e+04    -  1.00e+00 4.63e-03h  8
    1367  0.0000000e+00 1.86e+04 5.94e+17  -1.0 1.69e+04    -  1.00e+00 4.64e-03h  8
    1368  0.0000000e+00 1.85e+04 5.97e+17  -1.0 1.69e+04    -  1.00e+00 4.64e-03h  8
    1369  0.0000000e+00 1.90e+04 5.37e+16  -1.0 1.70e+04    -  1.00e+00 5.94e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1370  0.0000000e+00 9.95e+04 2.55e+17  -1.0 2.68e+04    -  1.00e+00 5.84e-01w  1
    WARNING: Problem in step computation; switching to emergency mode.
    1371  0.0000000e+00 1.85e+04 5.99e+17  -1.0 2.68e+04  20.0 1.00e+00 4.64e-03h  8
    1372  0.0000000e+00 1.84e+04 6.02e+17  -1.0 1.70e+04    -  1.00e+00 4.64e-03h  8
    1373  0.0000000e+00 1.83e+04 6.05e+17  -1.0 1.70e+04    -  1.00e+00 4.64e-03h  8
    1374  0.0000000e+00 1.82e+04 6.08e+17  -1.0 1.71e+04    -  1.00e+00 4.65e-03h  8
    1375  0.0000000e+00 1.81e+04 6.11e+17  -1.0 1.73e+04    -  1.00e+00 4.65e-03h  8
    1376  0.0000000e+00 1.80e+04 6.14e+17  -1.0 1.72e+04    -  1.00e+00 4.65e-03h  8
    1377  0.0000000e+00 1.80e+04 6.17e+17  -1.0 1.72e+04    -  1.00e+00 4.65e-03h  8
    1378  0.0000000e+00 1.79e+04 6.20e+17  -1.0 1.72e+04    -  1.00e+00 4.66e-03h  8
    1379  0.0000000e+00 1.78e+04 6.23e+17  -1.0 1.72e+04    -  1.00e+00 4.66e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1380  0.0000000e+00 1.77e+04 6.26e+17  -1.0 1.72e+04    -  1.00e+00 4.66e-03h  8
    1381  0.0000000e+00 2.14e+04 5.39e+16  -1.0 1.72e+04    -  1.00e+00 5.97e-01w  1
    1382  0.0000000e+00 1.11e+05 1.51e+17  -1.0 2.94e+04    -  1.00e+00 5.52e-01w  1
    WARNING: Problem in step computation; switching to emergency mode.
    1383  0.0000000e+00 1.76e+04 6.29e+17  -1.0 2.94e+04    -  1.00e+00 4.66e-03h  8
    1384  0.0000000e+00 1.75e+04 6.33e+17  -1.0 1.73e+04    -  1.00e+00 4.66e-03h  8
    1385  0.0000000e+00 1.75e+04 6.36e+17  -1.0 1.73e+04    -  1.00e+00 4.67e-03h  8
    1386  0.0000000e+00 1.74e+04 6.39e+17  -1.0 1.73e+04    -  1.00e+00 4.67e-03h  8
    1387  0.0000000e+00 1.73e+04 6.42e+17  -1.0 1.73e+04    -  1.00e+00 4.67e-03h  8
    1388  0.0000000e+00 1.72e+04 6.46e+17  -1.0 1.74e+04    -  1.00e+00 4.67e-03h  8
    1389  0.0000000e+00 1.71e+04 6.49e+17  -1.0 1.74e+04    -  1.00e+00 4.67e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1390  0.0000000e+00 1.71e+04 6.52e+17  -1.0 1.74e+04    -  1.00e+00 4.68e-03h  8
    1391  0.0000000e+00 1.70e+04 6.56e+17  -1.0 1.75e+04    -  1.00e+00 4.68e-03h  8
    1392  0.0000000e+00 1.69e+04 6.59e+17  -1.0 1.75e+04    -  1.00e+00 4.68e-03h  8
    1393  0.0000000e+00 2.38e+04 5.43e+16  -1.0 1.75e+04    -  1.00e+00 6.00e-01w  1
    1394  0.0000000e+00 1.21e+05 4.70e+16  -1.0 3.19e+04    -  1.00e+00 5.23e-01w  1
    WARNING: Problem in step computation; switching to emergency mode.
    1395  0.0000000e+00 1.68e+04 6.63e+17  -1.0 3.19e+04    -  1.00e+00 4.68e-03h  8
    1396  0.0000000e+00 1.67e+04 6.66e+17  -1.0 1.75e+04    -  1.00e+00 4.69e-03h  8
    1397  0.0000000e+00 1.67e+04 6.70e+17  -1.0 1.76e+04    -  1.00e+00 4.69e-03h  8
    1398  0.0000000e+00 1.66e+04 6.72e+17  -1.0 1.76e+04    -  1.00e+00 2.35e-03h  9
    1399  0.0000000e+00 1.66e+04 6.73e+17  -1.0 1.77e+04    -  1.00e+00 2.35e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1400  0.0000000e+00 1.65e+04 6.75e+17  -1.0 1.77e+04    -  1.00e+00 2.35e-03h  9
    1401  0.0000000e+00 1.65e+04 6.77e+17  -1.0 1.77e+04    -  1.00e+00 2.35e-03h  9
    1402  0.0000000e+00 1.65e+04 6.79e+17  -1.0 1.75e+04    -  1.00e+00 2.35e-03h  9
    1403  0.0000000e+00 1.64e+04 6.81e+17  -1.0 1.76e+04    -  1.00e+00 2.35e-03h  9
    1404  0.0000000e+00 1.64e+04 6.83e+17  -1.0 1.77e+04    -  1.00e+00 2.35e-03h  9
    1405  0.0000000e+00 2.59e+04 5.48e+16  -1.0 1.78e+04    -  1.00e+00 6.02e-01w  1
    1406  0.0000000e+00 1.28e+05 2.07e+16  -1.0 3.35e+04    -  1.00e+00 5.06e-01w  1
    WARNING: Problem in step computation; switching to emergency mode.
    1407  0.0000000e+00 1.63e+04 6.84e+17  -1.0 3.35e+04    -  1.00e+00 2.35e-03h  9
    1408  0.0000000e+00 1.63e+04 6.86e+17  -1.0 1.77e+04    -  1.00e+00 2.35e-03h  9
    1409  0.0000000e+00 1.63e+04 6.88e+17  -1.0 1.78e+04    -  1.00e+00 2.35e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1410  0.0000000e+00 1.62e+04 6.90e+17  -1.0 1.78e+04    -  1.00e+00 2.35e-03h  9
    1411  0.0000000e+00 1.62e+04 6.92e+17  -1.0 1.79e+04    -  1.00e+00 2.35e-03h  9
    1412  0.0000000e+00 1.62e+04 6.94e+17  -1.0 1.79e+04    -  1.00e+00 2.35e-03h  9
    1413  0.0000000e+00 1.61e+04 6.96e+17  -1.0 1.78e+04    -  1.00e+00 2.35e-03h  9
    1414  0.0000000e+00 1.61e+04 6.98e+17  -1.0 1.79e+04    -  1.00e+00 2.35e-03h  9
    1415  0.0000000e+00 1.60e+04 7.00e+17  -1.0 1.79e+04    -  1.00e+00 2.35e-03h  9
    1416  0.0000000e+00 1.60e+04 7.02e+17  -1.0 1.79e+04    -  1.00e+00 2.36e-03h  9
    1417  0.0000000e+00 2.72e+04 5.51e+16  -1.0 1.79e+04    -  1.00e+00 6.03e-01w  1
    1418  0.0000000e+00 1.32e+05 7.44e+16  -1.0 3.46e+04    -  1.00e+00 4.94e-01w  1
    WARNING: Problem in step computation; switching to emergency mode.
    1419  0.0000000e+00 1.60e+04 7.04e+17  -1.0 3.46e+04  20.0 1.00e+00 2.36e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1420  0.0000000e+00 1.59e+04 7.05e+17  -1.0 1.79e+04    -  1.00e+00 2.36e-03h  9
    1421  0.0000000e+00 1.59e+04 7.07e+17  -1.0 1.79e+04    -  1.00e+00 2.36e-03h  9
    1422  0.0000000e+00 1.59e+04 7.09e+17  -1.0 1.80e+04    -  1.00e+00 2.36e-03h  9
    1423  0.0000000e+00 1.58e+04 7.11e+17  -1.0 1.80e+04    -  1.00e+00 2.36e-03h  9
    1424  0.0000000e+00 1.58e+04 7.13e+17  -1.0 1.80e+04    -  1.00e+00 2.36e-03h  9
    1425  0.0000000e+00 1.57e+04 7.15e+17  -1.0 1.80e+04    -  1.00e+00 2.36e-03h  9
    1426  0.0000000e+00 1.57e+04 7.17e+17  -1.0 1.81e+04    -  1.00e+00 2.36e-03h  9
    1427  0.0000000e+00 1.57e+04 7.19e+17  -1.0 1.80e+04    -  1.00e+00 2.36e-03h  9
    1428  0.0000000e+00 1.56e+04 7.22e+17  -1.0 1.80e+04    -  1.00e+00 2.36e-03h  9
    1429  0.0000000e+00 2.84e+04 5.53e+16  -1.0 1.81e+04    -  1.00e+00 6.05e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1430  0.0000000e+00 1.36e+05 1.29e+17  -1.0 3.58e+04    -  1.00e+00 4.82e-01w  1
    WARNING: Problem in step computation; switching to emergency mode.
    1431  0.0000000e+00 1.56e+04 7.24e+17  -1.0 3.58e+04    -  1.00e+00 2.36e-03h  9
    1432  0.0000000e+00 1.56e+04 7.26e+17  -1.0 1.80e+04    -  1.00e+00 2.36e-03h  9
    1433  0.0000000e+00 1.55e+04 7.28e+17  -1.0 1.80e+04    -  1.00e+00 2.36e-03h  9
    1434  0.0000000e+00 1.55e+04 7.30e+17  -1.0 1.81e+04    -  1.00e+00 2.36e-03h  9
    1435  0.0000000e+00 1.54e+04 7.32e+17  -1.0 1.79e+04    -  1.00e+00 2.36e-03h  9
    1436  0.0000000e+00 1.54e+04 7.34e+17  -1.0 1.80e+04    -  1.00e+00 2.36e-03h  9
    1437  0.0000000e+00 1.54e+04 7.36e+17  -1.0 1.81e+04    -  1.00e+00 2.36e-03h  9
    1438  0.0000000e+00 1.53e+04 7.38e+17  -1.0 1.81e+04    -  1.00e+00 2.37e-03h  9
    1439  0.0000000e+00 1.53e+04 7.40e+17  -1.0 1.81e+04    -  1.00e+00 2.37e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1440  0.0000000e+00 1.53e+04 7.43e+17  -1.0 1.80e+04    -  1.00e+00 2.37e-03h  9
    1441  0.0000000e+00 2.93e+04 5.55e+16  -1.0 1.81e+04    -  1.00e+00 6.06e-01w  1
    1442  0.0000000e+00 1.40e+05 1.85e+17  -1.0 3.69e+04    -  1.00e+00 4.70e-01w  1
    WARNING: Problem in step computation; switching to emergency mode.
    1443  0.0000000e+00 1.52e+04 7.45e+17  -1.0 3.69e+04  20.0 1.00e+00 2.37e-03h  9
    1444  0.0000000e+00 1.52e+04 7.47e+17  -1.0 1.82e+04    -  1.00e+00 2.37e-03h  9
    1445  0.0000000e+00 1.52e+04 7.49e+17  -1.0 1.81e+04    -  1.00e+00 2.37e-03h  9
    1446  0.0000000e+00 1.51e+04 7.51e+17  -1.0 1.78e+04    -  1.00e+00 2.37e-03h  9
    1447  0.0000000e+00 1.51e+04 7.53e+17  -1.0 1.79e+04    -  1.00e+00 2.37e-03h  9
    1448  0.0000000e+00 1.51e+04 7.56e+17  -1.0 1.81e+04    -  1.00e+00 2.37e-03h  9
    1449  0.0000000e+00 1.50e+04 7.58e+17  -1.0 1.81e+04    -  1.00e+00 2.37e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1450  0.0000000e+00 1.50e+04 7.60e+17  -1.0 1.82e+04    -  1.00e+00 2.37e-03h  9
    1451  0.0000000e+00 1.49e+04 7.62e+17  -1.0 1.82e+04    -  1.00e+00 2.37e-03h  9
    1452  0.0000000e+00 1.49e+04 7.65e+17  -1.0 1.83e+04    -  1.00e+00 2.37e-03h  9
    1453  0.0000000e+00 3.09e+04 5.60e+16  -1.0 1.83e+04    -  1.00e+00 6.08e-01w  1
    1454  0.0000000e+00 1.44e+05 2.43e+17  -1.0 3.80e+04    -  1.00e+00 4.59e-01w  1
    1455  0.0000000e+00 1.44e+05 9.54e+20  -1.0 4.85e-04  20.0 9.93e-01 7.97e-04w  1
    1456  0.0000000e+00 1.49e+04 7.67e+17  -1.0 4.50e-06  19.5 1.00e+00 2.37e-03h  8
    1457  0.0000000e+00 1.48e+04 7.69e+17  -1.0 1.83e+04    -  1.00e+00 2.37e-03h  9
    1458  0.0000000e+00 1.48e+04 7.71e+17  -1.0 1.83e+04    -  1.00e+00 2.37e-03h  9
    1459  0.0000000e+00 1.48e+04 7.74e+17  -1.0 1.83e+04    -  1.00e+00 2.37e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1460  0.0000000e+00 1.47e+04 7.76e+17  -1.0 1.83e+04    -  1.00e+00 2.38e-03h  9
    1461  0.0000000e+00 1.47e+04 7.78e+17  -1.0 1.83e+04    -  1.00e+00 2.38e-03h  9
    1462  0.0000000e+00 1.47e+04 7.81e+17  -1.0 1.83e+04    -  1.00e+00 2.38e-03h  9
    1463  0.0000000e+00 1.46e+04 7.83e+17  -1.0 1.82e+04    -  1.00e+00 2.38e-03h  9
    1464  0.0000000e+00 1.46e+04 7.85e+17  -1.0 1.83e+04    -  1.00e+00 2.38e-03h  9
    1465  0.0000000e+00 1.46e+04 7.88e+17  -1.0 1.83e+04    -  1.00e+00 2.38e-03h  9
    1466  0.0000000e+00 3.21e+04 5.63e+16  -1.0 1.84e+04    -  1.00e+00 6.09e-01w  1
    1467  0.0000000e+00 1.48e+05 3.03e+17  -1.0 3.92e+04    -  1.00e+00 4.48e-01w  1
    1468  0.0000000e+00 1.48e+05 1.48e+19  -1.0 7.08e+05    -  1.41e-01 9.13e-04w  1
    1469  0.0000000e+00 1.45e+04 7.90e+17  -1.0 1.06e+04    -  1.00e+00 2.38e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1470  0.0000000e+00 1.45e+04 7.92e+17  -1.0 1.84e+04    -  1.00e+00 2.38e-03h  9
    1471  0.0000000e+00 1.45e+04 7.95e+17  -1.0 1.84e+04    -  1.00e+00 2.38e-03h  9
    1472  0.0000000e+00 1.44e+04 7.97e+17  -1.0 1.84e+04    -  1.00e+00 2.38e-03h  9
    1473  0.0000000e+00 1.44e+04 8.00e+17  -1.0 1.84e+04    -  1.00e+00 2.38e-03h  9
    1474  0.0000000e+00 1.44e+04 8.02e+17  -1.0 1.84e+04    -  1.00e+00 2.38e-03h  9
    1475  0.0000000e+00 1.43e+04 8.05e+17  -1.0 1.84e+04    -  1.00e+00 2.38e-03h  9
    1476  0.0000000e+00 1.43e+04 8.07e+17  -1.0 1.84e+04    -  1.00e+00 2.38e-03h  9
    1477  0.0000000e+00 1.42e+04 8.09e+17  -1.0 1.85e+04    -  1.00e+00 2.38e-03h  9
    1478  0.0000000e+00 1.42e+04 8.12e+17  -1.0 1.85e+04    -  1.00e+00 2.38e-03h  9
    1479  0.0000000e+00 3.35e+04 5.65e+16  -1.0 1.84e+04    -  1.00e+00 6.10e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1480  0.0000000e+00 1.52e+05 3.69e+17  -1.0 4.05e+04    -  1.00e+00 4.37e-01w  1
    1481  0.0000000e+00 1.52e+05 1.20e+20  -1.0 2.99e-05  19.9 9.93e-01 1.31e-02w  1
    1482  0.0000000e+00 1.42e+04 8.14e+17  -1.0 4.62e-06  19.5 1.00e+00 2.38e-03h  8
    1483  0.0000000e+00 1.41e+04 8.17e+17  -1.0 1.84e+04    -  1.00e+00 2.38e-03h  9
    1484  0.0000000e+00 1.41e+04 8.19e+17  -1.0 1.85e+04    -  1.00e+00 2.39e-03h  9
    1485  0.0000000e+00 1.41e+04 8.22e+17  -1.0 1.85e+04    -  1.00e+00 2.39e-03h  9
    1486  0.0000000e+00 1.40e+04 8.25e+17  -1.0 1.84e+04    -  1.00e+00 2.39e-03h  9
    1487  0.0000000e+00 1.40e+04 8.27e+17  -1.0 1.84e+04    -  1.00e+00 2.39e-03h  9
    1488  0.0000000e+00 1.40e+04 8.30e+17  -1.0 1.84e+04    -  1.00e+00 2.39e-03h  9
    1489  0.0000000e+00 1.39e+04 8.32e+17  -1.0 1.84e+04    -  1.00e+00 2.39e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1490  0.0000000e+00 1.39e+04 8.35e+17  -1.0 1.85e+04    -  1.00e+00 2.39e-03h  9
    1491  0.0000000e+00 1.39e+04 8.37e+17  -1.0 1.85e+04    -  1.00e+00 2.39e-03h  9
    1492  0.0000000e+00 3.47e+04 5.68e+16  -1.0 1.85e+04    -  1.00e+00 6.12e-01w  1
    1493  0.0000000e+00 1.55e+05 4.34e+17  -1.0 4.16e+04    -  1.00e+00 4.27e-01w  1
    1494  0.0000000e+00 1.55e+05 1.05e+20  -1.0 2.55e-05  19.9 9.94e-01 1.56e-02w  1
    1495  0.0000000e+00 1.38e+04 8.40e+17  -1.0 4.67e-06  19.4 1.00e+00 2.39e-03h  8
    1496  0.0000000e+00 1.38e+04 8.43e+17  -1.0 1.85e+04    -  1.00e+00 2.39e-03h  9
    1497  0.0000000e+00 1.38e+04 8.45e+17  -1.0 1.85e+04    -  1.00e+00 2.39e-03h  9
    1498  0.0000000e+00 1.37e+04 8.48e+17  -1.0 1.86e+04    -  1.00e+00 2.39e-03h  9
    1499  0.0000000e+00 1.37e+04 8.51e+17  -1.0 1.86e+04    -  1.00e+00 2.39e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1500  0.0000000e+00 1.37e+04 8.53e+17  -1.0 1.85e+04    -  1.00e+00 2.39e-03h  9
    1501  0.0000000e+00 1.36e+04 8.56e+17  -1.0 1.84e+04    -  1.00e+00 2.39e-03h  9
    1502  0.0000000e+00 1.36e+04 8.59e+17  -1.0 1.85e+04    -  1.00e+00 2.39e-03h  9
    1503  0.0000000e+00 1.36e+04 8.61e+17  -1.0 1.86e+04    -  1.00e+00 2.39e-03h  9
    1504  0.0000000e+00 1.36e+04 8.64e+17  -1.0 1.86e+04    -  1.00e+00 2.40e-03h  9
    1505  0.0000000e+00 3.56e+04 5.74e+16  -1.0 1.87e+04    -  1.00e+00 6.14e-01w  1
    1506  0.0000000e+00 1.57e+05 4.93e+17  -1.0 4.24e+04    -  1.00e+00 4.18e-01w  1
    1507  0.0000000e+00 1.57e+05 2.47e+19  -1.0 1.94e+05    -  8.17e-01 4.13e-03w  1
    1508  0.0000000e+00 1.35e+04 8.67e+17  -1.0 1.05e+04    -  1.00e+00 2.40e-03h  8
    1509  0.0000000e+00 1.35e+04 8.70e+17  -1.0 1.84e+04    -  1.00e+00 2.39e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1510  0.0000000e+00 1.35e+04 8.72e+17  -1.0 1.85e+04    -  1.00e+00 2.40e-03h  9
    1511  0.0000000e+00 1.34e+04 8.75e+17  -1.0 1.86e+04    -  1.00e+00 2.40e-03h  9
    1512  0.0000000e+00 1.34e+04 8.78e+17  -1.0 1.87e+04    -  1.00e+00 2.40e-03h  9
    1513  0.0000000e+00 1.34e+04 8.81e+17  -1.0 1.87e+04    -  1.00e+00 2.40e-03h  9
    1514  0.0000000e+00 1.33e+04 8.84e+17  -1.0 1.87e+04    -  1.00e+00 2.40e-03h  9
    1515  0.0000000e+00 1.33e+04 8.86e+17  -1.0 1.86e+04    -  1.00e+00 2.40e-03h  9
    1516  0.0000000e+00 1.33e+04 8.89e+17  -1.0 1.85e+04    -  1.00e+00 2.40e-03h  9
    1517  0.0000000e+00 1.32e+04 8.92e+17  -1.0 1.86e+04    -  1.00e+00 2.40e-03h  9
    1518  0.0000000e+00 3.66e+04 5.77e+16  -1.0 1.87e+04    -  1.00e+00 6.15e-01w  1
    1519  0.0000000e+00 1.59e+05 5.59e+17  -1.0 4.35e+04    -  1.00e+00 4.09e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1520  0.0000000e+00 1.59e+05 1.82e+19  -1.0 1.18e+05    -  1.00e+00 6.37e-03w  1
    1521  0.0000000e+00 1.32e+04 8.95e+17  -1.0 1.24e+04    -  1.00e+00 2.40e-03h  8
    1522  0.0000000e+00 1.32e+04 8.98e+17  -1.0 1.87e+04    -  1.00e+00 2.40e-03h  9
    1523  0.0000000e+00 1.31e+04 9.01e+17  -1.0 1.88e+04    -  1.00e+00 2.40e-03h  9
    1524  0.0000000e+00 1.31e+04 9.04e+17  -1.0 1.88e+04    -  1.00e+00 2.40e-03h  9
    1525  0.0000000e+00 1.31e+04 9.07e+17  -1.0 1.87e+04    -  1.00e+00 2.40e-03h  9
    1526  0.0000000e+00 1.30e+04 9.09e+17  -1.0 1.87e+04    -  1.00e+00 2.41e-03h  9
    1527  0.0000000e+00 1.30e+04 9.12e+17  -1.0 1.86e+04    -  1.00e+00 2.40e-03h  9
    1528  0.0000000e+00 1.30e+04 9.15e+17  -1.0 1.87e+04    -  1.00e+00 2.41e-03h  9
    1529  0.0000000e+00 1.29e+04 9.18e+17  -1.0 1.88e+04    -  1.00e+00 2.41e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1530  0.0000000e+00 1.29e+04 9.21e+17  -1.0 1.88e+04    -  1.00e+00 2.41e-03h  9
    1531  0.0000000e+00 3.81e+04 5.82e+16  -1.0 1.88e+04    -  1.00e+00 6.17e-01w  1
    1532  0.0000000e+00 1.62e+05 6.31e+17  -1.0 4.45e+04    -  1.00e+00 3.99e-01w  1
    1533  0.0000000e+00 1.61e+05 1.90e+19  -1.0 1.19e+05    -  1.00e+00 7.10e-03w  1
    1534  0.0000000e+00 1.29e+04 9.24e+17  -1.0 1.40e+04    -  1.00e+00 2.41e-03h  8
    1535  0.0000000e+00 1.29e+04 9.27e+17  -1.0 1.89e+04    -  1.00e+00 2.41e-03h  9
    1536  0.0000000e+00 1.28e+04 9.30e+17  -1.0 1.88e+04    -  1.00e+00 2.41e-03h  9
    1537  0.0000000e+00 1.28e+04 9.33e+17  -1.0 1.89e+04    -  1.00e+00 2.41e-03h  9
    1538  0.0000000e+00 1.28e+04 9.37e+17  -1.0 1.89e+04    -  1.00e+00 2.41e-03h  9
    1539  0.0000000e+00 1.27e+04 9.40e+17  -1.0 1.88e+04    -  1.00e+00 2.41e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1540  0.0000000e+00 1.27e+04 9.43e+17  -1.0 1.86e+04    -  1.00e+00 2.41e-03h  9
    1541  0.0000000e+00 1.27e+04 9.46e+17  -1.0 1.88e+04    -  1.00e+00 2.41e-03h  9
    1542  0.0000000e+00 1.26e+04 9.49e+17  -1.0 1.87e+04    -  1.00e+00 2.41e-03h  9
    1543  0.0000000e+00 1.26e+04 9.52e+17  -1.0 1.88e+04    -  1.00e+00 2.41e-03h  9
    1544  0.0000000e+00 3.99e+04 5.80e+16  -1.0 1.85e+04    -  1.00e+00 6.17e-01w  1
    1545  0.0000000e+00 1.65e+05 7.41e+17  -1.0 4.66e+04    -  1.00e+00 3.85e-01w  1
    1546  0.0000000e+00 1.64e+05 1.44e+19  -1.0 8.70e+04    -  1.00e+00 9.79e-03w  1
    1547  0.0000000e+00 1.26e+04 9.55e+17  -1.0 1.24e+04    -  1.00e+00 2.41e-03h  8
    1548  0.0000000e+00 1.25e+04 9.58e+17  -1.0 1.87e+04    -  1.00e+00 2.41e-03h  9
    1549  0.0000000e+00 1.25e+04 9.61e+17  -1.0 1.87e+04    -  1.00e+00 2.41e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1550  0.0000000e+00 1.25e+04 9.65e+17  -1.0 1.88e+04    -  1.00e+00 2.42e-03h  9
    1551  0.0000000e+00 1.25e+04 9.68e+17  -1.0 1.88e+04    -  1.00e+00 2.42e-03h  9
    1552  0.0000000e+00 1.24e+04 9.71e+17  -1.0 1.88e+04    -  1.00e+00 2.42e-03h  9
    1553  0.0000000e+00 1.24e+04 9.74e+17  -1.0 1.89e+04    -  1.00e+00 2.42e-03h  9
    1554  0.0000000e+00 1.24e+04 9.78e+17  -1.0 1.88e+04    -  1.00e+00 2.42e-03h  9
    1555  0.0000000e+00 1.23e+04 9.81e+17  -1.0 1.86e+04    -  1.00e+00 2.42e-03h  9
    1556  0.0000000e+00 1.23e+04 9.84e+17  -1.0 1.88e+04    -  1.00e+00 2.42e-03h  9
    1557  0.0000000e+00 4.00e+04 5.89e+16  -1.0 1.89e+04    -  1.00e+00 6.20e-01w  1
    1558  0.0000000e+00 1.65e+05 7.76e+17  -1.0 4.65e+04    -  1.00e+00 3.81e-01w  1
    1559  0.0000000e+00 1.64e+05 1.35e+19  -1.0 8.24e+04    -  1.00e+00 1.02e-02w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1560  0.0000000e+00 1.23e+04 9.87e+17  -1.0 1.16e+04    -  1.00e+00 2.42e-03h  8
    1561  0.0000000e+00 1.22e+04 9.91e+17  -1.0 1.89e+04    -  1.00e+00 2.42e-03h  9
    1562  0.0000000e+00 1.22e+04 9.94e+17  -1.0 1.89e+04    -  1.00e+00 2.42e-03h  9
    1563  0.0000000e+00 1.22e+04 9.97e+17  -1.0 1.90e+04    -  1.00e+00 2.42e-03h  9
    1564  0.0000000e+00 1.22e+04 1.00e+18  -1.0 1.90e+04    -  1.00e+00 2.42e-03h  9
    1565  0.0000000e+00 1.21e+04 1.00e+18  -1.0 1.90e+04    -  1.00e+00 2.42e-03h  9
    1566  0.0000000e+00 1.21e+04 1.01e+18  -1.0 1.89e+04    -  1.00e+00 2.42e-03h  9
    1567  0.0000000e+00 1.21e+04 1.01e+18  -1.0 1.89e+04    -  1.00e+00 2.42e-03h  9
    1568  0.0000000e+00 1.20e+04 1.01e+18  -1.0 1.85e+04    -  1.00e+00 2.42e-03h  9
    1569  0.0000000e+00 1.20e+04 1.02e+18  -1.0 1.87e+04    -  1.00e+00 2.42e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1570  0.0000000e+00 4.10e+04 5.91e+16  -1.0 1.87e+04    -  1.00e+00 6.21e-01w  1
    1571  0.0000000e+00 1.67e+05 8.63e+17  -1.0 4.79e+04    -  1.00e+00 3.71e-01w  1
    1572  0.0000000e+00 1.65e+05 1.17e+19  -1.0 7.06e+04    -  1.00e+00 1.19e-02w  1
    1573  0.0000000e+00 1.20e+04 1.02e+18  -1.0 1.21e+04    -  1.00e+00 2.43e-03h  8
    1574  0.0000000e+00 1.20e+04 1.02e+18  -1.0 1.89e+04    -  1.00e+00 2.43e-03h  9
    1575  0.0000000e+00 1.19e+04 1.03e+18  -1.0 1.89e+04    -  1.00e+00 2.43e-03h  9
    1576  0.0000000e+00 1.19e+04 1.03e+18  -1.0 1.90e+04    -  1.00e+00 2.43e-03h  9
    1577  0.0000000e+00 1.19e+04 1.04e+18  -1.0 1.90e+04    -  1.00e+00 2.43e-03h  9
    1578  0.0000000e+00 1.18e+04 1.04e+18  -1.0 1.90e+04    -  1.00e+00 2.43e-03h  9
    1579  0.0000000e+00 1.18e+04 1.04e+18  -1.0 1.90e+04    -  1.00e+00 2.43e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1580  0.0000000e+00 1.18e+04 1.05e+18  -1.0 1.90e+04    -  1.00e+00 2.43e-03h  9
    1581  0.0000000e+00 1.18e+04 1.05e+18  -1.0 1.87e+04    -  1.00e+00 2.43e-03h  9
    1582  0.0000000e+00 1.17e+04 1.05e+18  -1.0 1.89e+04    -  1.00e+00 2.43e-03h  9
    1583  0.0000000e+00 4.24e+04 5.98e+16  -1.0 1.90e+04    -  1.00e+00 6.23e-01w  1
    1584  0.0000000e+00 1.67e+05 9.31e+17  -1.0 4.85e+04    -  1.00e+00 3.64e-01w  1
    1585  0.0000000e+00 1.65e+05 1.06e+19  -1.0 6.37e+04    -  1.00e+00 1.31e-02w  1
    1586  0.0000000e+00 1.17e+04 1.06e+18  -1.0 3.01e+04    -  1.00e+00 2.43e-03h  8
    1587  0.0000000e+00 1.17e+04 1.06e+18  -1.0 1.91e+04    -  1.00e+00 2.43e-03h  9
    1588  0.0000000e+00 1.16e+04 1.06e+18  -1.0 1.91e+04    -  1.00e+00 2.43e-03h  9
    1589  0.0000000e+00 1.16e+04 1.07e+18  -1.0 1.91e+04    -  1.00e+00 2.44e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1590  0.0000000e+00 1.16e+04 1.07e+18  -1.0 1.86e+04    -  1.00e+00 2.43e-03h  9
    1591  0.0000000e+00 1.16e+04 1.07e+18  -1.0 1.89e+04    -  1.00e+00 2.44e-03h  9
    1592  0.0000000e+00 1.15e+04 1.08e+18  -1.0 1.87e+04    -  1.00e+00 2.43e-03h  9
    1593  0.0000000e+00 1.15e+04 1.08e+18  -1.0 1.89e+04    -  1.00e+00 2.44e-03h  9
    1594  0.0000000e+00 1.15e+04 1.09e+18  -1.0 1.90e+04    -  1.00e+00 2.44e-03h  9
    1595  0.0000000e+00 1.14e+04 1.09e+18  -1.0 1.90e+04    -  1.00e+00 2.44e-03h  9
    1596  0.0000000e+00 4.35e+04 6.02e+16  -1.0 1.91e+04    -  1.00e+00 6.25e-01w  1
    1597  0.0000000e+00 1.68e+05 1.01e+18  -1.0 4.94e+04    -  1.00e+00 3.55e-01w  1
    1598  0.0000000e+00 1.66e+05 9.47e+18  -1.0 5.68e+04    -  1.00e+00 1.45e-02w  1
    1599  0.0000000e+00 1.14e+04 1.09e+18  -1.0 1.69e+04    -  1.00e+00 2.44e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1600  0.0000000e+00 1.14e+04 1.10e+18  -1.0 1.91e+04    -  1.00e+00 2.44e-03h  9
    1601  0.0000000e+00 1.14e+04 1.10e+18  -1.0 1.91e+04    -  1.00e+00 2.44e-03h  9
    1602  0.0000000e+00 1.13e+04 1.10e+18  -1.0 1.92e+04    -  1.00e+00 2.44e-03h  9
    1603  0.0000000e+00 1.13e+04 1.11e+18  -1.0 1.90e+04    -  1.00e+00 2.44e-03h  9
    1604  0.0000000e+00 1.13e+04 1.11e+18  -1.0 1.91e+04    -  1.00e+00 2.44e-03h  9
    1605  0.0000000e+00 1.12e+04 1.12e+18  -1.0 1.90e+04    -  1.00e+00 2.44e-03h  9
    1606  0.0000000e+00 1.12e+04 1.12e+18  -1.0 1.91e+04    -  1.00e+00 2.44e-03h  9
    1607  0.0000000e+00 1.12e+04 1.12e+18  -1.0 1.91e+04    -  1.00e+00 2.45e-03h  9
    1608  0.0000000e+00 1.12e+04 1.13e+18  -1.0 1.92e+04    -  1.00e+00 2.45e-03h  9
    1609  0.0000000e+00 4.53e+04 6.05e+16  -1.0 1.91e+04    -  1.00e+00 6.26e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1610  0.0000000e+00 1.69e+05 1.12e+18  -1.0 5.06e+04    -  1.00e+00 3.45e-01w  1
    1611  0.0000000e+00 1.66e+05 8.50e+18  -1.0 5.06e+04    -  1.00e+00 1.62e-02w  1
    1612  0.0000000e+00 1.11e+04 1.13e+18  -1.0 1.48e+04    -  1.00e+00 2.45e-03h  8
    1613  0.0000000e+00 1.11e+04 1.14e+18  -1.0 1.91e+04    -  1.00e+00 2.45e-03h  9
    1614  0.0000000e+00 1.11e+04 1.14e+18  -1.0 1.91e+04    -  1.00e+00 2.45e-03h  9
    1615  0.0000000e+00 1.11e+04 1.14e+18  -1.0 1.91e+04    -  1.00e+00 2.45e-03h  9
    1616  0.0000000e+00 1.10e+04 1.15e+18  -1.0 1.90e+04    -  1.00e+00 2.45e-03h  9
    1617  0.0000000e+00 1.10e+04 1.15e+18  -1.0 1.91e+04    -  1.00e+00 2.45e-03h  9
    1618  0.0000000e+00 1.10e+04 1.16e+18  -1.0 1.92e+04    -  1.00e+00 2.45e-03h  9
    1619  0.0000000e+00 1.09e+04 1.16e+18  -1.0 1.92e+04    -  1.00e+00 2.45e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1620  0.0000000e+00 1.09e+04 1.16e+18  -1.0 1.92e+04    -  1.00e+00 2.45e-03h  9
    1621  0.0000000e+00 1.09e+04 1.17e+18  -1.0 1.91e+04    -  1.00e+00 1.23e-03h 10
    1622  0.0000000e+00 4.60e+04 6.10e+16  -1.0 1.92e+04    -  1.00e+00 6.28e-01w  1
    1623  0.0000000e+00 1.68e+05 1.19e+18  -1.0 5.12e+04    -  1.00e+00 3.38e-01w  1
    1624  0.0000000e+00 1.66e+05 8.46e+18  -1.0 4.99e+04    -  1.00e+00 1.69e-02w  1
    1625  0.0000000e+00 1.09e+04 1.17e+18  -1.0 2.06e+04    -  1.00e+00 1.23e-03h  9
    1626  0.0000000e+00 1.09e+04 1.17e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    1627  0.0000000e+00 1.09e+04 1.17e+18  -1.0 1.91e+04    -  1.00e+00 1.23e-03h 10
    1628  0.0000000e+00 1.09e+04 1.17e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    1629  0.0000000e+00 1.08e+04 1.18e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1630  0.0000000e+00 1.08e+04 1.18e+18  -1.0 1.93e+04    -  1.00e+00 1.23e-03h 10
    1631  0.0000000e+00 1.08e+04 1.18e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    1632  0.0000000e+00 1.08e+04 1.18e+18  -1.0 1.93e+04    -  1.00e+00 1.23e-03h 10
    1633  0.0000000e+00 1.08e+04 1.18e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    1634  0.0000000e+00 1.08e+04 1.19e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    1635  0.0000000e+00 4.72e+04 6.11e+16  -1.0 1.91e+04    -  1.00e+00 6.28e-01w  1
    1636  0.0000000e+00 1.68e+05 1.26e+18  -1.0 5.20e+04    -  1.00e+00 3.32e-01w  1
    1637  0.0000000e+00 1.66e+05 8.02e+18  -1.0 4.70e+04    -  1.00e+00 1.79e-02w  1
    1638  0.0000000e+00 1.08e+04 1.19e+18  -1.0 1.39e+04    -  1.00e+00 1.23e-03h  9
    1639  0.0000000e+00 1.07e+04 1.19e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1640  0.0000000e+00 1.07e+04 1.19e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    1641  0.0000000e+00 1.07e+04 1.19e+18  -1.0 1.90e+04    -  1.00e+00 1.23e-03h 10
    1642  0.0000000e+00 1.07e+04 1.20e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    1643  0.0000000e+00 1.07e+04 1.20e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    1644  0.0000000e+00 1.07e+04 1.20e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    1645  0.0000000e+00 1.07e+04 1.20e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    1646  0.0000000e+00 1.07e+04 1.21e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    1647  0.0000000e+00 1.06e+04 1.21e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    1648  0.0000000e+00 4.73e+04 6.14e+16  -1.0 1.93e+04    -  1.00e+00 6.30e-01w  1
    1649  0.0000000e+00 1.68e+05 1.28e+18  -1.0 5.21e+04    -  1.00e+00 3.30e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1650  0.0000000e+00 1.65e+05 7.76e+18  -1.0 4.60e+04    -  1.00e+00 1.80e-02w  1
    1651  0.0000000e+00 1.06e+04 1.21e+18  -1.0 1.41e+04    -  1.00e+00 1.23e-03h  9
    1652  0.0000000e+00 1.06e+04 1.21e+18  -1.0 1.93e+04    -  1.00e+00 1.23e-03h 10
    1653  0.0000000e+00 1.06e+04 1.21e+18  -1.0 1.93e+04    -  1.00e+00 1.23e-03h 10
    1654  0.0000000e+00 1.06e+04 1.22e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    1655  0.0000000e+00 1.06e+04 1.22e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    1656  0.0000000e+00 1.06e+04 1.22e+18  -1.0 1.91e+04    -  1.00e+00 1.23e-03h 10
    1657  0.0000000e+00 1.06e+04 1.22e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    1658  0.0000000e+00 1.05e+04 1.22e+18  -1.0 1.93e+04    -  1.00e+00 1.23e-03h 10
    1659  0.0000000e+00 1.05e+04 1.23e+18  -1.0 1.93e+04    -  1.00e+00 1.23e-03h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1660  0.0000000e+00 1.05e+04 1.23e+18  -1.0 1.93e+04    -  1.00e+00 1.23e-03h 10
    1661  0.0000000e+00 4.80e+04 6.17e+16  -1.0 1.93e+04    -  1.00e+00 6.31e-01w  1
    1662  0.0000000e+00 1.67e+05 1.32e+18  -1.0 5.25e+04    -  1.00e+00 3.25e-01w  1
    1663  0.0000000e+00 1.65e+05 7.01e+18  -1.0 4.22e+04    -  1.00e+00 1.87e-02w  1
    1664  0.0000000e+00 1.05e+04 1.23e+18  -1.0 1.43e+04    -  1.00e+00 1.23e-03h  9
    1665  0.0000000e+00 1.05e+04 1.23e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    1666  0.0000000e+00 1.05e+04 1.24e+18  -1.0 1.93e+04    -  1.00e+00 1.23e-03h 10
    1667  0.0000000e+00 1.05e+04 1.24e+18  -1.0 1.91e+04    -  1.00e+00 1.23e-03h 10
    1668  0.0000000e+00 1.04e+04 1.24e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    1669  0.0000000e+00 1.04e+04 1.24e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1670  0.0000000e+00 1.04e+04 1.24e+18  -1.0 1.88e+04    -  1.00e+00 1.23e-03h 10
    1671  0.0000000e+00 1.04e+04 1.25e+18  -1.0 1.90e+04    -  1.00e+00 1.23e-03h 10
    1672  0.0000000e+00 1.04e+04 1.25e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    1673  0.0000000e+00 1.04e+04 1.25e+18  -1.0 1.91e+04    -  1.00e+00 1.23e-03h 10
    1674  0.0000000e+00 4.80e+04 6.18e+16  -1.0 1.92e+04    -  1.00e+00 6.31e-01w  1
    1675  0.0000000e+00 1.67e+05 1.36e+18  -1.0 5.29e+04    -  1.00e+00 3.22e-01w  1
    1676  0.0000000e+00 1.64e+05 7.61e+18  -1.0 4.49e+04    -  1.00e+00 1.87e-02w  1
    1677  0.0000000e+00 1.04e+04 1.25e+18  -1.0 1.61e+04    -  1.00e+00 1.23e-03h  9
    1678  0.0000000e+00 1.04e+04 1.26e+18  -1.0 1.89e+04    -  1.00e+00 1.23e-03h 10
    1679  0.0000000e+00 1.03e+04 1.26e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1680  0.0000000e+00 1.03e+04 1.26e+18  -1.0 1.89e+04    -  1.00e+00 1.23e-03h 10
    1681  0.0000000e+00 1.03e+04 1.26e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    1682  0.0000000e+00 1.03e+04 1.26e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    1683  0.0000000e+00 1.03e+04 1.27e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    1684  0.0000000e+00 1.03e+04 1.27e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    1685  0.0000000e+00 1.03e+04 1.27e+18  -1.0 1.88e+04    -  1.00e+00 1.23e-03h 10
    1686  0.0000000e+00 1.03e+04 1.27e+18  -1.0 1.91e+04    -  1.00e+00 1.23e-03h 10
    1687  0.0000000e+00 4.84e+04 6.19e+16  -1.0 1.92e+04    -  1.00e+00 6.32e-01w  1
    1688  0.0000000e+00 1.67e+05 1.41e+18  -1.0 5.35e+04    -  1.00e+00 3.17e-01w  1
    1689  0.0000000e+00 1.64e+05 6.87e+18  -1.0 4.13e+04    -  1.00e+00 1.97e-02w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1690  0.0000000e+00 1.02e+04 1.28e+18  -1.0 1.82e+04    -  1.00e+00 1.23e-03h  9
    1691  0.0000000e+00 1.02e+04 1.28e+18  -1.0 1.93e+04    -  1.00e+00 1.23e-03h 10
    1692  0.0000000e+00 1.02e+04 1.28e+18  -1.0 1.92e+04    -  1.00e+00 1.23e-03h 10
    1693  0.0000000e+00 1.02e+04 1.28e+18  -1.0 1.93e+04    -  1.00e+00 1.24e-03h 10
    1694  0.0000000e+00 1.02e+04 1.29e+18  -1.0 1.93e+04    -  1.00e+00 1.24e-03h 10
    1695  0.0000000e+00 1.02e+04 1.29e+18  -1.0 1.93e+04    -  1.00e+00 1.24e-03h 10
    1696  0.0000000e+00 1.02e+04 1.29e+18  -1.0 1.93e+04    -  1.00e+00 1.24e-03h 10
    1697  0.0000000e+00 1.02e+04 1.29e+18  -1.0 1.92e+04    -  1.00e+00 1.24e-03h 10
    1698  0.0000000e+00 1.01e+04 1.29e+18  -1.0 1.93e+04    -  1.00e+00 1.24e-03h 10
    1699  0.0000000e+00 1.01e+04 1.30e+18  -1.0 1.92e+04    -  1.00e+00 1.24e-03h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1700  0.0000000e+00 4.90e+04 6.20e+16  -1.0 1.92e+04    -  1.00e+00 6.33e-01w  1
    1701  0.0000000e+00 1.66e+05 1.47e+18  -1.0 5.42e+04    -  1.00e+00 3.11e-01w  1
    1702  0.0000000e+00 1.63e+05 6.69e+18  -1.0 4.03e+04    -  1.00e+00 2.10e-02w  1
    1703  0.0000000e+00 1.01e+04 1.30e+18  -1.0 1.89e+04    -  1.00e+00 1.24e-03h  9
    1704  0.0000000e+00 1.01e+04 1.30e+18  -1.0 1.92e+04    -  1.00e+00 1.24e-03h 10
    1705  0.0000000e+00 1.01e+04 1.30e+18  -1.0 1.93e+04    -  1.00e+00 1.24e-03h 10
    1706  0.0000000e+00 1.01e+04 1.31e+18  -1.0 1.92e+04    -  1.00e+00 1.24e-03h 10
    1707  0.0000000e+00 1.01e+04 1.31e+18  -1.0 1.92e+04    -  1.00e+00 1.24e-03h 10
    1708  0.0000000e+00 1.01e+04 1.31e+18  -1.0 1.92e+04    -  1.00e+00 1.24e-03h 10
    1709  0.0000000e+00 1.00e+04 1.31e+18  -1.0 1.89e+04    -  1.00e+00 1.24e-03h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1710  0.0000000e+00 1.00e+04 1.32e+18  -1.0 1.89e+04    -  1.00e+00 1.24e-03h 10
    1711  0.0000000e+00 1.00e+04 1.32e+18  -1.0 1.90e+04    -  1.00e+00 1.24e-03h 10
    1712  0.0000000e+00 1.00e+04 1.32e+18  -1.0 1.88e+04    -  1.00e+00 1.24e-03h 10
    1713  0.0000000e+00 4.75e+04 6.14e+16  -1.0 1.89e+04    -  1.00e+00 6.33e-01w  1
    1714  0.0000000e+00 1.67e+05 1.56e+18  -1.0 5.58e+04    -  1.00e+00 3.04e-01w  1
    1715  0.0000000e+00 1.63e+05 5.24e+18  -1.0 3.34e+04    -  1.00e+00 2.67e-02w  1
    1716  0.0000000e+00 9.99e+03 1.32e+18  -1.0 1.80e+04    -  1.00e+00 1.24e-03h  9
    1717  0.0000000e+00 9.98e+03 1.33e+18  -1.0 1.89e+04    -  1.00e+00 1.24e-03h 10
    1718  0.0000000e+00 9.97e+03 1.33e+18  -1.0 1.89e+04    -  1.00e+00 1.24e-03h 10
    1719  0.0000000e+00 9.96e+03 1.33e+18  -1.0 1.89e+04    -  1.00e+00 1.24e-03h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1720  0.0000000e+00 9.95e+03 1.33e+18  -1.0 1.88e+04    -  1.00e+00 1.24e-03h 10
    1721  0.0000000e+00 9.93e+03 1.33e+18  -1.0 1.87e+04    -  1.00e+00 1.24e-03h 10
    1722  0.0000000e+00 9.92e+03 1.34e+18  -1.0 1.85e+04    -  1.00e+00 1.24e-03h 10
    1723  0.0000000e+00 9.91e+03 1.34e+18  -1.0 1.85e+04    -  1.00e+00 1.24e-03h 10
    1724  0.0000000e+00 9.90e+03 1.34e+18  -1.0 1.82e+04    -  1.00e+00 1.24e-03h 10
    1725  0.0000000e+00 9.87e+03 1.35e+18  -1.0 1.82e+04    -  1.00e+00 2.47e-03h  9
    1726  0.0000000e+00 4.29e+04 5.90e+16  -1.0 1.80e+04    -  1.00e+00 6.33e-01w  1
    1727  0.0000000e+00 1.73e+05 1.83e+18  -1.0 6.19e+04    -  1.00e+00 2.82e-01w  1
    1728  0.0000000e+00 1.65e+05 4.39e+18  -1.0 2.83e+04    -  1.00e+00 4.86e-02w  1
    1729  0.0000000e+00 9.85e+03 1.35e+18  -1.0 1.72e+04    -  1.00e+00 2.47e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1730  0.0000000e+00 9.82e+03 1.36e+18  -1.0 1.79e+04    -  1.00e+00 2.47e-03h  9
    1731  0.0000000e+00 9.80e+03 1.36e+18  -1.0 1.78e+04    -  1.00e+00 2.47e-03h  9
    1732  0.0000000e+00 9.77e+03 1.37e+18  -1.0 1.75e+04    -  1.00e+00 2.47e-03h  9
    1733  0.0000000e+00 9.75e+03 1.37e+18  -1.0 1.73e+04    -  1.00e+00 2.47e-03h  9
    1734  0.0000000e+00 9.73e+03 1.38e+18  -1.0 1.70e+04    -  1.00e+00 2.47e-03h  9
    1735  0.0000000e+00 9.70e+03 1.38e+18  -1.0 1.68e+04    -  1.00e+00 2.47e-03h  9
    1736  0.0000000e+00 9.68e+03 1.39e+18  -1.0 1.68e+04    -  1.00e+00 2.47e-03h  9
    1737  0.0000000e+00 9.65e+03 1.39e+18  -1.0 1.66e+04    -  1.00e+00 2.47e-03h  9
    1738  0.0000000e+00 9.63e+03 1.40e+18  -1.0 1.64e+04    -  1.00e+00 2.47e-03h  9
    1739  0.0000000e+00 3.34e+04 5.42e+16  -1.0 1.62e+04    -  1.00e+00 6.32e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1740  0.0000000e+00 1.89e+05 2.58e+18  -1.0 7.97e+04    -  1.00e+00 2.35e-01w  1
    1741  0.0000000e+00 1.71e+05 3.39e+18  -1.0 2.08e+04    -  1.00e+00 1.10e-01w  1
    1742  0.0000000e+00 9.61e+03 1.40e+18  -1.0 1.56e+04    -  1.00e+00 2.47e-03h  8
    1743  0.0000000e+00 9.58e+03 1.41e+18  -1.0 1.60e+04    -  1.00e+00 2.47e-03h  9
    1744  0.0000000e+00 9.56e+03 1.41e+18  -1.0 1.58e+04    -  1.00e+00 2.47e-03h  9
    1745  0.0000000e+00 9.54e+03 1.42e+18  -1.0 1.57e+04    -  1.00e+00 2.47e-03h  9
    1746  0.0000000e+00 9.51e+03 1.42e+18  -1.0 1.55e+04    -  1.00e+00 2.47e-03h  9
    1747  0.0000000e+00 9.49e+03 1.43e+18  -1.0 1.54e+04    -  1.00e+00 2.47e-03h  9
    1748  0.0000000e+00 9.47e+03 1.43e+18  -1.0 1.51e+04    -  1.00e+00 2.47e-03h  9
    1749  0.0000000e+00 9.44e+03 1.44e+18  -1.0 1.50e+04    -  1.00e+00 2.47e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1750  0.0000000e+00 9.42e+03 1.44e+18  -1.0 1.49e+04    -  1.00e+00 2.47e-03h  9
    1751  0.0000000e+00 9.40e+03 1.45e+18  -1.0 1.48e+04    -  1.00e+00 2.47e-03h  9
    1752  0.0000000e+00 2.63e+04 4.98e+16  -1.0 1.47e+04    -  1.00e+00 6.32e-01w  1
    1753  0.0000000e+00 2.12e+05 3.98e+18  -1.0 1.11e+05    -  1.00e+00 1.81e-01w  1
    1754  0.0000000e+00 1.77e+05 3.78e+18  -1.0 1.86e+04    -  1.00e+00 1.99e-01w  1
    1755  0.0000000e+00 9.37e+03 1.45e+18  -1.0 1.31e+04    -  1.00e+00 2.47e-03h  8
    1756  0.0000000e+00 9.35e+03 1.46e+18  -1.0 1.45e+04    -  1.00e+00 2.47e-03h  9
    1757  0.0000000e+00 9.33e+03 1.46e+18  -1.0 1.43e+04    -  1.00e+00 2.47e-03h  9
    1758  0.0000000e+00 9.30e+03 1.47e+18  -1.0 1.42e+04    -  1.00e+00 2.47e-03h  9
    1759  0.0000000e+00 9.28e+03 1.47e+18  -1.0 1.41e+04    -  1.00e+00 2.47e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1760  0.0000000e+00 9.26e+03 1.48e+18  -1.0 1.37e+04    -  1.00e+00 2.47e-03h  9
    1761  0.0000000e+00 9.24e+03 1.49e+18  -1.0 1.38e+04    -  1.00e+00 2.47e-03h  9
    1762  0.0000000e+00 9.21e+03 1.49e+18  -1.0 1.39e+04    -  1.00e+00 2.47e-03h  9
    1763  0.0000000e+00 9.19e+03 1.50e+18  -1.0 1.36e+04    -  1.00e+00 2.47e-03h  9
    1764  0.0000000e+00 9.17e+03 1.50e+18  -1.0 1.37e+04    -  1.00e+00 2.47e-03h  9
    1765  0.0000000e+00 2.27e+04 4.78e+16  -1.0 1.37e+04    -  1.00e+00 6.32e-01w  1
    1766  0.0000000e+00 2.28e+05 4.66e+18  -1.0 1.42e+05    -  8.95e-01 1.49e-01w  1
    1767  0.0000000e+00 1.81e+05 3.69e+18  -1.0 1.71e+04    -  1.00e+00 2.38e-01w  1
    1768  0.0000000e+00 9.14e+03 1.51e+18  -1.0 1.22e+04    -  1.00e+00 2.47e-03h  8
    1769  0.0000000e+00 9.12e+03 1.51e+18  -1.0 1.37e+04    -  1.00e+00 2.47e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1770  0.0000000e+00 9.10e+03 1.52e+18  -1.0 1.35e+04    -  1.00e+00 2.47e-03h  9
    1771  0.0000000e+00 9.08e+03 1.52e+18  -1.0 1.35e+04    -  1.00e+00 2.47e-03h  9
    1772  0.0000000e+00 9.05e+03 1.53e+18  -1.0 1.34e+04    -  1.00e+00 2.47e-03h  9
    1773  0.0000000e+00 9.03e+03 1.54e+18  -1.0 1.34e+04    -  1.00e+00 2.47e-03h  9
    1774  0.0000000e+00 9.01e+03 1.54e+18  -1.0 1.32e+04    -  1.00e+00 2.47e-03h  9
    1775  0.0000000e+00 8.97e+03 1.55e+18  -1.0 1.33e+04    -  1.00e+00 4.95e-03h  8
    1776  0.0000000e+00 8.92e+03 1.56e+18  -1.0 1.32e+04    -  1.00e+00 4.95e-03h  8
    1777  0.0000000e+00 8.88e+03 1.58e+18  -1.0 1.32e+04    -  1.00e+00 4.95e-03h  8
    1778  0.0000000e+00 2.08e+04 4.72e+16  -1.0 1.31e+04    -  1.00e+00 6.34e-01w  1
    1779  0.0000000e+00 2.33e+05 4.66e+18  -1.0 1.58e+05    -  8.01e-01 1.36e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1780  0.0000000e+00 1.85e+05 3.36e+18  -1.0 1.67e+04    -  1.00e+00 2.40e-01w  1
    1781  0.0000000e+00 8.83e+03 1.59e+18  -1.0 1.17e+04    -  1.00e+00 4.95e-03h  7
    1782  0.0000000e+00 8.79e+03 1.60e+18  -1.0 1.30e+04    -  1.00e+00 4.96e-03h  8
    1783  0.0000000e+00 8.75e+03 1.61e+18  -1.0 1.30e+04    -  1.00e+00 4.96e-03h  8
    1784  0.0000000e+00 8.70e+03 1.62e+18  -1.0 1.29e+04    -  1.00e+00 4.96e-03h  8
    1785  0.0000000e+00 8.66e+03 1.63e+18  -1.0 1.28e+04    -  1.00e+00 4.96e-03h  8
    1786  0.0000000e+00 8.62e+03 1.65e+18  -1.0 1.27e+04    -  1.00e+00 4.97e-03h  8
    1787  0.0000000e+00 8.57e+03 1.66e+18  -1.0 1.26e+04    -  1.00e+00 4.97e-03h  8
    1788  0.0000000e+00 8.53e+03 1.67e+18  -1.0 1.26e+04    -  1.00e+00 4.97e-03h  8
    1789  0.0000000e+00 8.49e+03 1.68e+18  -1.0 1.26e+04    -  1.00e+00 4.97e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1790  0.0000000e+00 8.45e+03 1.70e+18  -1.0 1.25e+04    -  1.00e+00 4.98e-03h  8
    1791  0.0000000e+00 1.90e+04 4.71e+16  -1.0 1.24e+04    -  1.00e+00 6.37e-01w  1
    1792  0.0000000e+00 2.34e+05 4.68e+18  -1.0 1.64e+05    -  7.71e-01 1.31e-01w  1
    1793  0.0000000e+00 1.86e+05 2.87e+18  -1.0 1.64e+04    -  1.00e+00 2.38e-01w  1
    1794  0.0000000e+00 8.40e+03 1.71e+18  -1.0 1.15e+04    -  1.00e+00 4.98e-03h  7
    1795  0.0000000e+00 8.36e+03 1.72e+18  -1.0 1.24e+04    -  1.00e+00 4.98e-03h  8
    1796  0.0000000e+00 8.32e+03 1.74e+18  -1.0 1.21e+04    -  1.00e+00 4.98e-03h  8
    1797  0.0000000e+00 8.28e+03 1.75e+18  -1.0 1.22e+04    -  1.00e+00 4.99e-03h  8
    1798  0.0000000e+00 8.24e+03 1.76e+18  -1.0 1.22e+04    -  1.00e+00 4.99e-03h  8
    1799  0.0000000e+00 8.20e+03 1.78e+18  -1.0 1.22e+04    -  1.00e+00 4.99e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1800  0.0000000e+00 8.16e+03 1.79e+18  -1.0 1.20e+04    -  1.00e+00 5.00e-03h  8
    1801  0.0000000e+00 8.12e+03 1.80e+18  -1.0 1.20e+04    -  1.00e+00 5.00e-03h  8
    1802  0.0000000e+00 8.07e+03 1.82e+18  -1.0 1.20e+04    -  1.00e+00 5.00e-03h  8
    1803  0.0000000e+00 8.03e+03 1.83e+18  -1.0 1.19e+04    -  1.00e+00 5.00e-03h  8
    1804  0.0000000e+00 1.76e+04 4.67e+16  -1.0 1.18e+04    -  1.00e+00 6.41e-01w  1
    1805  0.0000000e+00 2.34e+05 4.68e+18  -1.0 1.66e+05    -  7.59e-01 1.29e-01w  1
    1806  0.0000000e+00 1.86e+05 2.33e+18  -1.0 1.62e+04    -  1.00e+00 2.36e-01w  1
    1807  0.0000000e+00 7.99e+03 1.84e+18  -1.0 1.17e+04    -  1.00e+00 5.01e-03h  7
    1808  0.0000000e+00 7.95e+03 1.86e+18  -1.0 1.17e+04    -  1.00e+00 5.01e-03h  8
    1809  0.0000000e+00 7.91e+03 1.87e+18  -1.0 1.17e+04    -  1.00e+00 5.01e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1810  0.0000000e+00 7.87e+03 1.89e+18  -1.0 1.17e+04    -  1.00e+00 5.02e-03h  8
    1811  0.0000000e+00 7.84e+03 1.90e+18  -1.0 1.16e+04    -  1.00e+00 5.02e-03h  8
    1812  0.0000000e+00 7.80e+03 1.91e+18  -1.0 1.16e+04    -  1.00e+00 5.02e-03h  8
    1813  0.0000000e+00 7.76e+03 1.93e+18  -1.0 1.15e+04    -  1.00e+00 5.02e-03h  8
    1814  0.0000000e+00 7.72e+03 1.94e+18  -1.0 1.14e+04    -  1.00e+00 5.03e-03h  8
    1815  0.0000000e+00 7.68e+03 1.96e+18  -1.0 1.14e+04    -  1.00e+00 5.03e-03h  8
    1816  0.0000000e+00 7.64e+03 1.97e+18  -1.0 1.13e+04    -  1.00e+00 5.03e-03h  8
    1817  0.0000000e+00 1.62e+04 4.61e+16  -1.0 1.13e+04    -  1.00e+00 6.45e-01w  1
    1818  0.0000000e+00 2.33e+05 4.63e+18  -1.0 1.67e+05    -  7.55e-01 1.28e-01w  1
    1819  0.0000000e+00 1.85e+05 1.65e+18  -1.0 1.59e+04    -  1.00e+00 2.33e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1820  0.0000000e+00 7.60e+03 1.99e+18  -1.0 1.12e+04    -  1.00e+00 5.04e-03h  7
    1821  0.0000000e+00 7.56e+03 2.00e+18  -1.0 1.12e+04    -  1.00e+00 5.04e-03h  8
    1822  0.0000000e+00 7.53e+03 2.02e+18  -1.0 1.12e+04    -  1.00e+00 5.04e-03h  8
    1823  0.0000000e+00 7.49e+03 2.04e+18  -1.0 1.11e+04    -  1.00e+00 5.05e-03h  8
    1824  0.0000000e+00 7.45e+03 2.05e+18  -1.0 1.11e+04    -  1.00e+00 5.05e-03h  8
    1825  0.0000000e+00 7.41e+03 2.07e+18  -1.0 1.10e+04    -  1.00e+00 5.05e-03h  8
    1826  0.0000000e+00 7.38e+03 2.08e+18  -1.0 1.10e+04    -  1.00e+00 5.06e-03h  8
    1827  0.0000000e+00 7.34e+03 2.10e+18  -1.0 1.09e+04    -  1.00e+00 5.06e-03h  8
    1828  0.0000000e+00 7.30e+03 2.11e+18  -1.0 1.09e+04    -  1.00e+00 5.06e-03h  8
    1829  0.0000000e+00 7.26e+03 2.13e+18  -1.0 1.07e+04    -  1.00e+00 5.06e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1830  0.0000000e+00 1.50e+04 4.49e+16  -1.0 1.08e+04    -  1.00e+00 6.49e-01w  1
    1831  0.0000000e+00 2.33e+05 4.55e+18  -1.0 1.69e+05    -  7.46e-01 1.27e-01w  1
    1832  0.0000000e+00 1.86e+05 9.69e+17  -1.0 1.60e+04    -  1.00e+00 2.30e-01w  1
    1833  0.0000000e+00 7.23e+03 2.15e+18  -1.0 1.10e+04    -  1.00e+00 5.07e-03h  7
    1834  0.0000000e+00 7.19e+03 2.16e+18  -1.0 1.07e+04    -  1.00e+00 5.07e-03h  8
    1835  0.0000000e+00 7.15e+03 2.18e+18  -1.0 1.07e+04    -  1.00e+00 5.07e-03h  8
    1836  0.0000000e+00 7.12e+03 2.20e+18  -1.0 1.06e+04    -  1.00e+00 5.08e-03h  8
    1837  0.0000000e+00 7.08e+03 2.21e+18  -1.0 1.06e+04    -  1.00e+00 5.08e-03h  8
    1838  0.0000000e+00 7.05e+03 2.23e+18  -1.0 1.05e+04    -  1.00e+00 5.08e-03h  8
    1839  0.0000000e+00 7.01e+03 2.25e+18  -1.0 1.04e+04    -  1.00e+00 5.09e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1840  0.0000000e+00 6.97e+03 2.27e+18  -1.0 1.04e+04    -  1.00e+00 5.09e-03h  8
    1841  0.0000000e+00 6.94e+03 2.28e+18  -1.0 1.04e+04    -  1.00e+00 5.09e-03h  8
    1842  0.0000000e+00 6.90e+03 2.30e+18  -1.0 1.03e+04    -  1.00e+00 5.10e-03h  8
    1843  0.0000000e+00 1.40e+04 4.32e+16  -1.0 1.03e+04    -  1.00e+00 6.53e-01w  1
    1844  0.0000000e+00 2.32e+05 4.44e+18  -1.0 1.71e+05    -  7.37e-01 1.25e-01w  1
    1845  0.0000000e+00 1.86e+05 8.84e+16  -1.0 1.55e+04    -  1.00e+00 2.26e-01w  1
    1846  0.0000000e+00 6.87e+03 2.32e+18  -1.0 1.10e+04    -  1.00e+00 5.10e-03h  7
    1847  0.0000000e+00 6.83e+03 2.34e+18  -1.0 1.03e+04    -  1.00e+00 5.10e-03h  8
    1848  0.0000000e+00 6.80e+03 2.36e+18  -1.0 1.02e+04    -  1.00e+00 5.11e-03h  8
    1849  0.0000000e+00 6.76e+03 2.37e+18  -1.0 1.01e+04    -  1.00e+00 5.11e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1850  0.0000000e+00 6.73e+03 2.39e+18  -1.0 1.01e+04    -  1.00e+00 5.11e-03h  8
    1851  0.0000000e+00 6.70e+03 2.41e+18  -1.0 1.00e+04    -  1.00e+00 5.12e-03h  8
    1852  0.0000000e+00 6.66e+03 2.43e+18  -1.0 1.00e+04    -  1.00e+00 5.12e-03h  8
    1853  0.0000000e+00 6.63e+03 2.45e+18  -1.0 9.96e+03    -  1.00e+00 5.12e-03h  8
    1854  0.0000000e+00 6.59e+03 2.47e+18  -1.0 9.91e+03    -  1.00e+00 5.13e-03h  8
    1855  0.0000000e+00 6.56e+03 2.49e+18  -1.0 9.83e+03    -  1.00e+00 5.13e-03h  8
    1856  0.0000000e+00 1.30e+04 4.07e+16  -1.0 9.82e+03    -  1.00e+00 6.57e-01w  1
    1857  0.0000000e+00 2.32e+05 4.27e+18  -1.0 1.74e+05    -  7.24e-01 1.23e-01w  1
    1858  0.0000000e+00 1.86e+05 8.39e+17  -1.0 1.53e+04    -  1.00e+00 2.22e-01w  1
    1859  0.0000000e+00 6.53e+03 2.51e+18  -1.0 1.07e+04    -  1.00e+00 5.14e-03h  7
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1860  0.0000000e+00 6.49e+03 2.53e+18  -1.0 9.70e+03    -  1.00e+00 5.14e-03h  8
    1861  0.0000000e+00 6.46e+03 2.55e+18  -1.0 9.71e+03    -  1.00e+00 5.14e-03h  8
    1862  0.0000000e+00 6.43e+03 2.57e+18  -1.0 9.67e+03    -  1.00e+00 5.15e-03h  8
    1863  0.0000000e+00 6.39e+03 2.59e+18  -1.0 9.63e+03    -  1.00e+00 5.15e-03h  8
    1864  0.0000000e+00 6.36e+03 2.61e+18  -1.0 9.56e+03    -  1.00e+00 5.15e-03h  8
    1865  0.0000000e+00 6.33e+03 2.63e+18  -1.0 9.51e+03    -  1.00e+00 5.16e-03h  8
    1866  0.0000000e+00 6.29e+03 2.65e+18  -1.0 9.49e+03    -  1.00e+00 5.16e-03h  8
    1867  0.0000000e+00 6.26e+03 2.67e+18  -1.0 9.43e+03    -  1.00e+00 5.16e-03h  8
    1868  0.0000000e+00 6.23e+03 2.69e+18  -1.0 9.25e+03    -  1.00e+00 5.17e-03h  8
    1869  0.0000000e+00 1.19e+04 3.74e+16  -1.0 9.34e+03    -  1.00e+00 6.62e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1870  0.0000000e+00 2.32e+05 4.04e+18  -1.0 1.78e+05    -  7.08e-01 1.20e-01w  1
    1871  0.0000000e+00 1.87e+05 1.89e+18  -1.0 1.51e+04    -  1.00e+00 2.18e-01w  1
    1872  0.0000000e+00 6.20e+03 2.71e+18  -1.0 1.06e+04    -  1.00e+00 5.17e-03h  7
    1873  0.0000000e+00 6.17e+03 2.73e+18  -1.0 9.30e+03    -  1.00e+00 5.17e-03h  8
    1874  0.0000000e+00 6.13e+03 2.75e+18  -1.0 9.25e+03    -  1.00e+00 5.18e-03h  8
    1875  0.0000000e+00 6.10e+03 2.78e+18  -1.0 9.10e+03    -  1.00e+00 5.18e-03h  8
    1876  0.0000000e+00 6.07e+03 2.80e+18  -1.0 8.81e+03    -  1.00e+00 5.18e-03h  8
    1877  0.0000000e+00 6.04e+03 2.82e+18  -1.0 9.07e+03    -  1.00e+00 5.19e-03h  8
    1878  0.0000000e+00 6.01e+03 2.84e+18  -1.0 9.08e+03    -  1.00e+00 5.19e-03h  8
    1879  0.0000000e+00 5.98e+03 2.86e+18  -1.0 9.07e+03    -  1.00e+00 5.20e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1880  0.0000000e+00 5.95e+03 2.89e+18  -1.0 9.04e+03    -  1.00e+00 5.20e-03h  8
    1881  0.0000000e+00 5.91e+03 2.91e+18  -1.0 8.99e+03    -  1.00e+00 5.21e-03h  8
    1882  0.0000000e+00 1.12e+04 3.34e+16  -1.0 8.95e+03    -  1.00e+00 6.67e-01w  1
    1883  0.0000000e+00 2.32e+05 3.81e+18  -1.0 1.82e+05    -  6.94e-01 1.17e-01w  1
    1884  0.0000000e+00 1.87e+05 3.03e+18  -1.0 1.49e+04    -  1.00e+00 2.14e-01w  1
    1885  0.0000000e+00 5.88e+03 2.93e+18  -1.0 1.05e+04    -  1.00e+00 5.21e-03h  7
    1886  0.0000000e+00 5.85e+03 2.96e+18  -1.0 8.88e+03    -  1.00e+00 5.21e-03h  8
    1887  0.0000000e+00 5.82e+03 2.98e+18  -1.0 8.83e+03    -  1.00e+00 5.22e-03h  8
    1888  0.0000000e+00 5.79e+03 3.00e+18  -1.0 8.79e+03    -  1.00e+00 5.22e-03h  8
    1889  0.0000000e+00 5.76e+03 3.03e+18  -1.0 8.78e+03    -  1.00e+00 5.23e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1890  0.0000000e+00 5.73e+03 3.05e+18  -1.0 8.68e+03    -  1.00e+00 5.23e-03h  8
    1891  0.0000000e+00 5.70e+03 3.08e+18  -1.0 8.56e+03    -  1.00e+00 5.23e-03h  8
    1892  0.0000000e+00 5.67e+03 3.10e+18  -1.0 8.60e+03    -  1.00e+00 5.24e-03h  8
    1893  0.0000000e+00 5.64e+03 3.12e+18  -1.0 8.62e+03    -  1.00e+00 5.24e-03h  8
    1894  0.0000000e+00 5.61e+03 3.15e+18  -1.0 8.56e+03    -  1.00e+00 5.24e-03h  8
    1895  0.0000000e+00 1.05e+04 2.84e+16  -1.0 8.55e+03    -  1.00e+00 6.72e-01w  1
    1896  0.0000000e+00 2.32e+05 3.49e+18  -1.0 1.87e+05    -  6.77e-01 1.14e-01w  1
    1897  0.0000000e+00 1.88e+05 4.31e+18  -1.0 1.47e+04    -  1.00e+00 2.09e-01w  1
    1898  0.0000000e+00 5.58e+03 3.17e+18  -1.0 1.04e+04    -  1.00e+00 5.25e-03h  7
    1899  0.0000000e+00 5.55e+03 3.20e+18  -1.0 8.51e+03    -  1.00e+00 5.25e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1900  0.0000000e+00 5.52e+03 3.23e+18  -1.0 8.47e+03    -  1.00e+00 5.26e-03h  8
    1901  0.0000000e+00 5.50e+03 3.25e+18  -1.0 8.43e+03    -  1.00e+00 5.26e-03h  8
    1902  0.0000000e+00 5.47e+03 3.28e+18  -1.0 8.37e+03    -  1.00e+00 5.27e-03h  8
    1903  0.0000000e+00 5.44e+03 3.30e+18  -1.0 8.35e+03    -  1.00e+00 5.27e-03h  8
    1904  0.0000000e+00 5.41e+03 3.33e+18  -1.0 8.25e+03    -  1.00e+00 5.27e-03h  8
    1905  0.0000000e+00 5.38e+03 3.36e+18  -1.0 8.25e+03    -  1.00e+00 5.28e-03h  8
    1906  0.0000000e+00 5.35e+03 3.38e+18  -1.0 8.21e+03    -  1.00e+00 5.28e-03h  8
    1907  0.0000000e+00 5.32e+03 3.41e+18  -1.0 8.17e+03    -  1.00e+00 5.29e-03h  8
    1908  0.0000000e+00 9.76e+03 2.22e+16  -1.0 8.15e+03    -  1.00e+00 6.77e-01w  1
    1909  0.0000000e+00 2.32e+05 3.12e+18  -1.0 1.92e+05    -  6.57e-01 1.11e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1910  0.0000000e+00 1.89e+05 5.71e+18  -1.0 1.45e+04    -  1.00e+00 2.04e-01w  1
    1911  0.0000000e+00 5.27e+03 3.47e+18  -1.0 1.02e+04    -  1.00e+00 1.06e-02h  6
    1912  0.0000000e+00 5.21e+03 3.52e+18  -1.0 8.01e+03    -  1.00e+00 1.06e-02h  7
    1913  0.0000000e+00 5.16e+03 3.58e+18  -1.0 7.97e+03    -  1.00e+00 1.06e-02h  7
    1914  0.0000000e+00 5.10e+03 3.64e+18  -1.0 7.86e+03    -  1.00e+00 1.06e-02h  7
    1915  0.0000000e+00 5.05e+03 3.70e+18  -1.0 7.83e+03    -  1.00e+00 1.06e-02h  7
    1916  0.0000000e+00 5.00e+03 3.76e+18  -1.0 7.76e+03    -  1.00e+00 1.07e-02h  7
    1917  0.0000000e+00 4.94e+03 3.82e+18  -1.0 7.70e+03    -  1.00e+00 1.07e-02h  7
    1918  0.0000000e+00 4.89e+03 3.88e+18  -1.0 7.61e+03    -  1.00e+00 1.07e-02h  7
    1919  0.0000000e+00 4.84e+03 3.94e+18  -1.0 7.55e+03    -  1.00e+00 1.07e-02h  7
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1920  0.0000000e+00 4.79e+03 4.01e+18  -1.0 7.49e+03    -  1.00e+00 1.07e-02h  7
    1921  0.0000000e+00 8.46e+03 5.67e+15  -1.0 7.42e+03    -  1.00e+00 6.88e-01w  1
    1922  0.0000000e+00 2.33e+05 2.16e+18  -1.0 2.05e+05    -  6.17e-01 1.04e-01w  1
    1923  0.0000000e+00 1.92e+05 8.97e+18  -1.0 1.41e+04    -  1.00e+00 1.93e-01w  1
    1924  0.0000000e+00 4.73e+03 4.07e+18  -1.0 1.00e+04    -  1.00e+00 1.08e-02h  6
    1925  0.0000000e+00 4.68e+03 4.14e+18  -1.0 7.35e+03    -  1.00e+00 1.08e-02h  7
    1926  0.0000000e+00 4.63e+03 4.21e+18  -1.0 7.26e+03    -  1.00e+00 1.08e-02h  7
    1927  0.0000000e+00 4.58e+03 4.28e+18  -1.0 7.12e+03    -  1.00e+00 1.08e-02h  7
    1928  0.0000000e+00 4.53e+03 4.35e+18  -1.0 7.10e+03    -  1.00e+00 1.08e-02h  7
    1929  0.0000000e+00 4.48e+03 4.42e+18  -1.0 7.08e+03    -  1.00e+00 1.09e-02h  7
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1930  0.0000000e+00 4.44e+03 4.49e+18  -1.0 7.01e+03    -  1.00e+00 1.09e-02h  7
    1931  0.0000000e+00 4.39e+03 4.56e+18  -1.0 6.95e+03    -  1.00e+00 1.09e-02h  7
    1932  0.0000000e+00 4.34e+03 4.64e+18  -1.0 6.87e+03    -  1.00e+00 1.09e-02h  7
    1933  0.0000000e+00 4.29e+03 4.72e+18  -1.0 6.80e+03    -  1.00e+00 1.09e-02h  7
    1934  0.0000000e+00 7.37e+03 1.78e+16  -1.0 6.74e+03    -  1.00e+00 7.01e-01w  1
    1935  0.0000000e+00 2.34e+05 9.46e+17  -1.0 2.21e+05    -  5.73e-01 9.71e-02w  1
    1936  0.0000000e+00 1.95e+05 1.29e+19  -1.0 1.38e+04    -  1.00e+00 1.81e-01w  1
    1937  0.0000000e+00 4.25e+03 4.79e+18  -1.0 9.82e+03    -  1.00e+00 1.09e-02h  6
    1938  0.0000000e+00 4.20e+03 4.87e+18  -1.0 6.67e+03    -  1.00e+00 1.10e-02h  7
    1939  0.0000000e+00 4.15e+03 4.95e+18  -1.0 6.63e+03    -  1.00e+00 1.10e-02h  7
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1940  0.0000000e+00 4.11e+03 5.03e+18  -1.0 6.57e+03    -  1.00e+00 1.10e-02h  7
    1941  0.0000000e+00 4.06e+03 5.12e+18  -1.0 6.52e+03    -  1.00e+00 1.10e-02h  7
    1942  0.0000000e+00 4.02e+03 5.20e+18  -1.0 6.42e+03    -  1.00e+00 1.10e-02h  7
    1943  0.0000000e+00 3.97e+03 5.29e+18  -1.0 6.37e+03    -  1.00e+00 1.11e-02h  7
    1944  0.0000000e+00 3.93e+03 5.37e+18  -1.0 6.33e+03    -  1.00e+00 1.11e-02h  7
    1945  0.0000000e+00 3.89e+03 5.46e+18  -1.0 6.28e+03    -  1.00e+00 1.11e-02h  7
    1946  0.0000000e+00 3.84e+03 5.55e+18  -1.0 6.20e+03    -  1.00e+00 1.11e-02h  7
    1947  0.0000000e+00 6.43e+03 5.01e+16  -1.0 6.16e+03    -  1.00e+00 7.14e-01w  1
    1948  0.0000000e+00 2.35e+05 6.61e+17  -1.0 2.39e+05    -  5.31e-01 9.01e-02w  1
    1949  0.0000000e+00 1.98e+05 1.49e+19  -1.0 1.36e+04    -  8.72e-01 1.67e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1950  0.0000000e+00 3.80e+03 5.65e+18  -1.0 9.61e+03    -  1.00e+00 1.12e-02h  6
    1951  0.0000000e+00 3.76e+03 5.74e+18  -1.0 6.11e+03    -  1.00e+00 1.12e-02h  7
    1952  0.0000000e+00 3.72e+03 5.83e+18  -1.0 6.04e+03    -  1.00e+00 1.12e-02h  7
    1953  0.0000000e+00 3.68e+03 5.93e+18  -1.0 5.95e+03    -  1.00e+00 1.12e-02h  7
    1954  0.0000000e+00 3.63e+03 6.03e+18  -1.0 5.74e+03    -  1.00e+00 1.12e-02h  7
    1955  0.0000000e+00 3.59e+03 6.13e+18  -1.0 5.82e+03    -  1.00e+00 1.13e-02h  7
    1956  0.0000000e+00 3.55e+03 6.23e+18  -1.0 5.78e+03    -  1.00e+00 1.13e-02h  7
    1957  0.0000000e+00 3.51e+03 6.33e+18  -1.0 5.76e+03    -  1.00e+00 1.13e-02h  7
    1958  0.0000000e+00 3.47e+03 6.44e+18  -1.0 5.73e+03    -  1.00e+00 1.13e-02h  7
    1959  0.0000000e+00 3.43e+03 6.55e+18  -1.0 5.67e+03    -  1.00e+00 1.13e-02h  7
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1960  0.0000000e+00 5.64e+03 9.41e+16  -1.0 5.62e+03    -  1.00e+00 7.28e-01w  1
    1961  0.0000000e+00 2.36e+05 2.66e+18  -1.0 2.59e+05    -  4.90e-01 8.30e-02w  1
    1962  0.0000000e+00 2.02e+05 1.66e+19  -1.0 1.34e+04    -  7.42e-01 1.51e-01w  1
    1963  0.0000000e+00 3.40e+03 6.66e+18  -1.0 9.49e+03    -  1.00e+00 1.14e-02h  6
    1964  0.0000000e+00 3.36e+03 6.77e+18  -1.0 5.57e+03    -  1.00e+00 1.14e-02h  7
    1965  0.0000000e+00 3.32e+03 6.88e+18  -1.0 5.49e+03    -  1.00e+00 1.14e-02h  7
    1966  0.0000000e+00 3.28e+03 6.99e+18  -1.0 5.43e+03    -  1.00e+00 1.14e-02h  7
    1967  0.0000000e+00 3.24e+03 7.11e+18  -1.0 5.41e+03    -  1.00e+00 1.15e-02h  7
    1968  0.0000000e+00 3.21e+03 7.23e+18  -1.0 5.37e+03    -  1.00e+00 1.15e-02h  7
    1969  0.0000000e+00 3.17e+03 7.35e+18  -1.0 5.32e+03    -  1.00e+00 1.15e-02h  7
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1970  0.0000000e+00 3.13e+03 7.47e+18  -1.0 5.27e+03    -  1.00e+00 1.15e-02h  7
    1971  0.0000000e+00 3.10e+03 7.60e+18  -1.0 5.21e+03    -  1.00e+00 1.16e-02h  7
    1972  0.0000000e+00 3.06e+03 7.72e+18  -1.0 5.14e+03    -  1.00e+00 1.16e-02h  7
    1973  0.0000000e+00 4.94e+03 1.53e+17  -1.0 5.10e+03    -  1.00e+00 7.43e-01w  1
    1974  0.0000000e+00 2.36e+05 5.08e+18  -1.0 2.84e+05    -  4.48e-01 7.59e-02w  1
    1975  0.0000000e+00 2.06e+05 1.89e+19  -1.0 1.34e+04    -  6.36e-01 1.34e-01w  1
    1976  0.0000000e+00 3.03e+03 7.85e+18  -1.0 9.41e+03    -  1.00e+00 1.16e-02h  6
    1977  0.0000000e+00 2.99e+03 7.98e+18  -1.0 5.07e+03    -  1.00e+00 1.16e-02h  7
    1978  0.0000000e+00 2.96e+03 8.12e+18  -1.0 5.03e+03    -  1.00e+00 1.17e-02h  7
    1979  0.0000000e+00 2.92e+03 8.25e+18  -1.0 4.98e+03    -  1.00e+00 1.17e-02h  7
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1980  0.0000000e+00 2.89e+03 8.39e+18  -1.0 4.93e+03    -  1.00e+00 1.17e-02h  7
    1981  0.0000000e+00 2.85e+03 8.53e+18  -1.0 4.88e+03    -  1.00e+00 1.17e-02h  7
    1982  0.0000000e+00 2.82e+03 8.67e+18  -1.0 4.83e+03    -  1.00e+00 1.17e-02h  7
    1983  0.0000000e+00 2.76e+03 8.97e+18  -1.0 4.78e+03    -  1.00e+00 2.35e-02h  6
    1984  0.0000000e+00 2.69e+03 9.28e+18  -1.0 4.20e+03    -  1.00e+00 2.36e-02h  6
    1985  0.0000000e+00 2.63e+03 9.60e+18  -1.0 4.17e+03    -  1.00e+00 2.37e-02h  6
    1986  0.0000000e+00 3.55e+03 2.64e+17  -1.0 3.98e+03    -  1.00e+00 7.63e-01w  1
    1987  0.0000000e+00 2.51e+05 8.88e+18  -1.0 3.38e+05    -  3.75e-01 6.57e-02w  1
    1988  0.0000000e+00 2.25e+05 2.26e+19  -1.0 1.35e+04    -  5.23e-01 1.09e-01w  1
    1989  0.0000000e+00 2.57e+03 9.93e+18  -1.0 9.07e+03    -  1.00e+00 2.38e-02h  5
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1990  0.0000000e+00 2.50e+03 1.03e+19  -1.0 3.86e+03    -  1.00e+00 2.39e-02h  6
    1991  0.0000000e+00 2.44e+03 1.06e+19  -1.0 3.75e+03    -  1.00e+00 2.40e-02h  6
    1992  0.0000000e+00 2.39e+03 1.10e+19  -1.0 3.59e+03    -  1.00e+00 2.42e-02h  6
    1993  0.0000000e+00 2.33e+03 1.14e+19  -1.0 3.44e+03    -  1.00e+00 2.43e-02h  6
    1994  0.0000000e+00 2.22e+03 1.22e+19  -1.0 3.26e+03    -  1.00e+00 4.87e-02h  5
    1995  0.0000000e+00 2.53e+03 4.49e+17  -1.0 6.81e+02    -  1.00e+00 7.85e-01h  1
    1996  0.0000000e+00 2.60e+05 1.63e+17  -1.0 1.33e+05    -  1.72e-01 1.72e-01s 20
    1997  0.0000000e+00 2.55e+05 1.62e+17  -1.0 3.99e+04    -  2.01e-02 2.01e-02s 20
    1998r 0.0000000e+00 2.55e+05 1.00e+03   1.8 0.00e+00    -  0.00e+00 0.00e+00R  1
    1999r 0.0000000e+00 2.03e+05 4.43e+04   1.8 8.46e+04    -  1.42e-06 2.10e-04f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2000r 0.0000000e+00 1.43e+05 4.56e+04   1.8 6.56e+04    -  8.61e-05 3.08e-04f  1
    2001r 0.0000000e+00 1.43e+05 2.02e+05   1.8 2.43e-01   6.0 2.42e-03 1.09e-01f  1
    2002  0.0000000e+00 1.40e+05 7.93e+01  -1.0 1.21e+05    -  6.49e-02 2.21e-02h  1
    2003  0.0000000e+00 1.40e+05 1.98e+02  -1.0 1.97e+05    -  1.11e-02 3.44e-04h  1
    2004  0.0000000e+00 1.41e+05 4.06e+03  -1.0 2.43e+06    -  1.72e-03 2.10e-02f  1
    2005  0.0000000e+00 1.41e+05 3.76e+03  -1.0 2.99e+06    -  1.13e-02 2.37e-03h  2
    2006  0.0000000e+00 1.42e+05 3.85e+03  -1.0 3.09e+06    -  1.45e-02 4.72e-03h  2
    2007  0.0000000e+00 1.42e+05 3.61e+03  -1.0 3.34e+06    -  1.37e-02 7.25e-04h  4
    2008  0.0000000e+00 1.41e+05 2.57e+03  -1.0 3.38e+06    -  1.22e-02 2.67e-03h  3
    2009  0.0000000e+00 1.41e+05 5.77e+03  -1.0 3.56e+06    -  3.44e-02 2.72e-03h  3
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2010  0.0000000e+00 1.41e+05 1.22e+04  -1.0 3.79e+06    -  1.06e-02 5.07e-05h  8
    2011  0.0000000e+00 1.41e+05 1.32e+05  -1.0 3.79e+06    -  1.48e-01 1.25e-05h 10
    2012r 0.0000000e+00 1.41e+05 9.99e+02   1.5 0.00e+00    -  0.00e+00 3.90e-07R 15
    2013r 0.0000000e+00 5.48e+04 1.02e+03   1.5 4.00e+04    -  1.67e-03 9.97e-04f  1
    2014  0.0000000e+00 5.48e+04 1.95e+03  -1.0 7.68e+05    -  2.20e-01 1.70e-04h  1
    2015  0.0000000e+00 5.48e+04 3.77e+03  -1.0 4.23e+06    -  9.69e-03 3.78e-04h  3
    2016  0.0000000e+00 5.48e+04 2.35e+04  -1.0 4.28e+06    -  2.59e-02 1.33e-04h  4
    2017  0.0000000e+00 5.48e+04 6.08e+05  -1.0 4.30e+06    -  2.38e-02 2.86e-05h  6
    2018r 0.0000000e+00 5.48e+04 9.99e+02   1.1 0.00e+00    -  0.00e+00 4.31e-07R 12
    2019r 0.0000000e+00 5.45e+04 9.99e+02   1.1 1.22e+04    -  1.19e-03 1.15e-03f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020  0.0000000e+00 5.45e+04 2.20e+03  -1.0 9.82e+05    -  1.03e-01 6.47e-05h  1
    2021  0.0000000e+00 5.44e+04 2.20e+03  -1.0 4.68e+06    -  8.14e-04 8.14e-04s 12
    2022r 0.0000000e+00 5.44e+04 9.99e+02   0.5 0.00e+00    -  0.00e+00 0.00e+00R  1
    2023r 0.0000000e+00 5.43e+04 9.97e+02   0.5 4.63e+03    -  7.09e-03 9.91e-04f  1
    2024r 0.0000000e+00 5.43e+04 9.99e+02  -0.3 0.00e+00    -  0.00e+00 3.04e-07R  7
    2025r 0.0000000e+00 5.50e+04 9.88e+02  -0.3 5.23e+03    -  1.14e-02 1.07e-03f  1
    2026  0.0000000e+00 5.49e+04 9.49e+01  -1.0 1.61e+06    -  8.66e-02 5.31e-04h  2
    2027  0.0000000e+00 5.49e+04 9.47e+01  -1.0 4.94e+06    -  1.53e-03 1.53e-03s 13
    2028  0.0000000e+00 5.48e+04 1.56e+02  -1.0 1.19e+07    -  1.39e-03 1.39e-03s 13
    2029r 0.0000000e+00 5.48e+04 9.99e+02   0.9 0.00e+00    -  0.00e+00 0.00e+00R  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2030r 0.0000000e+00 5.47e+04 9.92e+02   0.9 8.42e+03    -  7.19e-03 9.08e-05f  1
    2031r 0.0000000e+00 4.96e+04 9.89e+02   0.9 1.37e+03    -  2.30e-03 8.46e-03f  1
    2032  0.0000000e+00 4.96e+04 3.50e+02  -1.0 2.84e+05    -  1.63e-01 2.96e-04h  2
    2033  0.0000000e+00 4.95e+04 3.49e+02  -1.0 5.23e+06    -  1.57e-03 1.57e-03s 13
    2034  0.0000000e+00 4.94e+04 3.48e+02  -1.0 2.66e+06    -  1.66e-03 1.66e-03s 13
    2035  0.0000000e+00 4.92e+04 3.45e+02  -1.0 2.82e+06    -  5.22e-03 5.22e-03s 13
    2036  0.0000000e+00 5.17e+04 4.91e+02  -1.0 3.23e+06    -  2.95e-03 2.95e-03s 13
    2037r 0.0000000e+00 5.17e+04 9.99e+02   1.1 0.00e+00    -  0.00e+00 0.00e+00R  1
    2038r 0.0000000e+00 4.82e+04 9.98e+02   1.1 1.67e+04    -  2.36e-03 9.92e-04f  1
    2039r 0.0000000e+00 4.63e+04 9.95e+02   1.1 1.54e+03    -  2.32e-03 3.44e-03f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2040r 0.0000000e+00 4.67e+04 9.85e+02   1.1 1.08e+03    -  1.08e-02 7.02e-03f  1
    2041  0.0000000e+00 4.67e+04 1.34e+03  -1.0 5.13e+05    -  5.59e-02 3.43e-04h  2
    2042  0.0000000e+00 4.65e+04 1.34e+03  -1.0 1.58e+06    -  3.58e-03 3.58e-03s 14
    2043r 0.0000000e+00 4.65e+04 9.99e+02  -0.3 0.00e+00    -  0.00e+00 0.00e+00R  1
    2044r 0.0000000e+00 4.67e+04 9.50e+02  -0.3 9.80e+02    -  7.49e-02 9.90e-04f  1
    2045r 0.0000000e+00 4.89e+04 9.04e+02  -0.3 2.89e+02    -  2.17e-02 5.81e-02f  1
    2046  0.0000000e+00 4.89e+04 2.05e+02  -1.0 2.17e+05    -  5.09e-02 5.35e-04h  1
    2047  0.0000000e+00 4.88e+04 2.05e+02  -1.0 1.25e+06    -  1.02e-03 1.02e-03s 12
    2048  0.0000000e+00 4.88e+04 2.05e+02  -1.0 1.13e+06    -  1.16e-03 1.16e-03s 12
    2049  0.0000000e+00 4.87e+04 2.05e+02  -1.0 1.14e+06    -  1.17e-03 1.17e-03s 12
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2050r 0.0000000e+00 4.87e+04 9.99e+02  -1.0 0.00e+00    -  0.00e+00 0.00e+00R  1
    2051r 0.0000000e+00 4.89e+04 9.85e+02  -1.0 8.85e+02    -  7.20e-02 1.29e-03f  1
    2052r 0.0000000e+00 4.89e+04 9.95e+02  -1.0 0.00e+00    -  0.00e+00 4.24e-07R  5
    2053r 0.0000000e+00 4.97e+04 9.78e+02  -1.0 1.59e+03    -  1.24e-01 5.39e-03f  1
    2054  0.0000000e+00 4.96e+04 4.12e+00  -1.0 2.88e+05    -  2.96e-03 2.96e-03s 14
    2055  0.0000000e+00 4.95e+04 5.77e+00  -1.0 4.39e+05    -  2.05e-03 2.05e-03s 14
    2056r 0.0000000e+00 4.95e+04 9.97e+02  -1.0 0.00e+00    -  0.00e+00 0.00e+00R  1
    2057r 0.0000000e+00 4.98e+04 9.87e+02  -1.0 1.44e+03    -  2.24e-01 3.02e-03f  1
    2058r 0.0000000e+00 5.21e+04 9.07e+02  -1.0 1.82e+03    -  4.82e-02 1.15e-01f  1
    2059  0.0000000e+00 5.19e+04 1.76e+00  -1.0 8.90e+04    -  8.73e-03 4.85e-03f  7
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2060  0.0000000e+00 5.09e+04 1.03e+01  -1.0 1.79e+05    -  4.46e-02 3.49e-02f  3
    2061  0.0000000e+00 5.07e+04 1.67e+01  -1.0 2.59e+05    -  2.50e-02 2.61e-03f  5
    2062  0.0000000e+00 5.07e+04 3.72e+01  -1.0 3.81e+05    -  4.62e-02 1.91e-03f  6
    2063  0.0000000e+00 1.06e+05 4.10e+02  -1.0 4.39e+05    -  2.89e-02 5.29e-02f  1
    2064  0.0000000e+00 1.06e+05 3.21e+02  -1.0 1.17e+06    -  2.01e-02 1.31e-02h  2
    2065  0.0000000e+00 1.06e+05 3.17e+02  -1.0 1.56e+06    -  3.65e-02 1.32e-03h  4
    2066  0.0000000e+00 1.06e+05 3.15e+02  -1.0 1.67e+06    -  2.24e-02 5.44e-04h  5
    2067  0.0000000e+00 1.06e+05 3.10e+02  -1.0 1.69e+06    -  6.11e-02 5.05e-04h  5
    2068  0.0000000e+00 1.06e+05 1.21e+03  -1.0 1.50e+06    -  2.24e-02 8.25e-03h  1
    2069  0.0000000e+00 1.06e+05 1.20e+03  -1.0 2.98e+05    -  9.24e-03 2.47e-04h  5
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2070  0.0000000e+00 1.06e+05 1.19e+03  -1.0 3.17e+05    -  1.02e-02 2.65e-04h  6
    2071  0.0000000e+00 1.06e+05 1.19e+03  -1.0 3.49e+05    -  1.33e-02 1.99e-04h  7
    2072  0.0000000e+00 1.06e+05 1.59e+03  -1.0 3.78e+05    -  3.02e-02 2.95e-04h  7
    2073  0.0000000e+00 1.06e+05 2.42e+03  -1.0 4.31e+05    -  1.68e-02 6.21e-05h 10
    2074  0.0000000e+00 1.06e+05 3.23e+03  -1.0 4.45e+05    -  1.47e-02 8.13e-05h 10
    2075  0.0000000e+00 2.56e+05 3.11e+04  -1.0 4.64e+05    -  3.15e-02 2.00e-02f  2
    2076  0.0000000e+00 2.56e+05 3.11e+04  -1.0 3.68e+05    -  3.42e-02 4.21e-04h  7
    2077  0.0000000e+00 2.55e+05 3.11e+04  -1.0 3.54e+05    -  4.03e-02 4.16e-04h  7
    2078  0.0000000e+00 2.56e+05 3.00e+04  -1.0 5.18e+05    -  5.89e-03 5.23e-03h  3
    2079  0.0000000e+00 2.57e+05 1.39e+04  -1.0 5.92e+05    -  2.72e-02 1.29e-02w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2080  0.0000000e+00 2.57e+05 1.07e+14  -1.0 1.22e-05  18.9 9.90e-01 1.00e+00w  1
    2081  0.0000000e+00 2.57e+05 2.07e+14  -1.0 7.15e-05  18.5 9.90e-01 1.00e+00w  1
    2082  0.0000000e+00 2.55e+05 2.93e+04  -1.0 2.88e-04  18.0 2.72e-02 3.23e-03h  2
    2083  0.0000000e+00 2.54e+05 2.79e+04  -1.0 1.51e+06    -  2.36e-03 3.39e-03h  2
    2084  0.0000000e+00 2.54e+05 2.85e+12  -1.0 2.63e-05  17.0 9.90e-01 1.00e+00h  1
    2085  0.0000000e+00 2.54e+05 5.69e+12  -1.0 1.58e-04  16.6 9.90e-01 1.00e+00h  1
    2086r 0.0000000e+00 2.54e+05 9.99e+02   1.8 0.00e+00  16.1 0.00e+00 4.20e-07R  7
    2087r 0.0000000e+00 9.96e+04 6.88e+03   1.8 8.20e+04    -  6.01e-02 1.01e-03f  1
    2088  0.0000000e+00 9.96e+04 3.64e+02  -1.0 6.11e+05    -  3.37e-02 6.06e-05h  9
    2089  0.0000000e+00 9.96e+04 4.53e+02  -1.0 6.12e+05    -  8.45e-02 6.03e-05h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2090  0.0000000e+00 9.96e+04 9.17e+02  -1.0 6.12e+05    -  6.20e-02 6.00e-05h  9
    2091  0.0000000e+00 9.96e+04 1.05e+04  -1.0 6.12e+05    -  2.75e-01 5.98e-05h  9
    2092  0.0000000e+00 9.96e+04 3.25e+04  -1.0 6.11e+05    -  3.69e-02 5.95e-05h  9
    2093  0.0000000e+00 9.96e+04 8.14e+05  -1.0 6.12e+05    -  4.20e-01 5.92e-05h  9
    2094  0.0000000e+00 9.96e+04 2.47e+06  -1.0 6.12e+05    -  3.53e-02 5.89e-05h  9
    2095  0.0000000e+00 9.96e+04 4.64e+07  -1.0 6.13e+05    -  2.76e-01 5.86e-05h  9
    2096  0.0000000e+00 9.96e+04 1.59e+08  -1.0 6.13e+05    -  3.71e-02 5.83e-05h  9
    2097  0.0000000e+00 9.96e+04 2.01e+09  -1.0 6.14e+05    -  1.77e-01 5.79e-05h  9
    2098  0.0000000e+00 2.12e+05 8.29e+10  -1.0 6.14e+05    -  3.99e-02 1.48e-02w  1
    2099  0.0000000e+00 2.17e+05 3.57e+10  -1.0 7.62e+04    -  2.34e-03 2.71e-02w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2100  0.0000000e+00 2.18e+05 2.55e+10  -1.0 8.54e+05    -  1.94e-03 3.53e-04w  1
    2101  0.0000000e+00 9.96e+04 7.30e+09  -1.0 5.11e+06    -  3.99e-02 5.76e-05h  8
    2102  0.0000000e+00 9.96e+04 2.73e+11  -1.0 6.15e+05    -  7.87e-01 5.73e-05h  9
    2103  0.0000000e+00 9.88e+04 4.28e+11  -1.0 3.49e+06    -  8.03e-03 7.31e-03h  2
    2104  0.0000000e+00 9.88e+04 1.05e+12  -1.0 2.62e-04  15.6 9.81e-01 1.00e+00h  1
    2105r 0.0000000e+00 9.88e+04 1.00e+03   1.2 0.00e+00  15.1 0.00e+00 3.12e-07R  6
    2106r 0.0000000e+00 9.94e+04 4.66e+03   1.2 1.86e+04    -  1.57e-01 9.92e-04f  1
    2107  0.0000000e+00 9.94e+04 2.61e+01  -1.0 8.42e+06    -  2.38e-03 9.14e-05h  1
    2108  0.0000000e+00 9.94e+04 2.86e+01  -1.0 3.08e+07    -  2.31e-04 1.32e-04f  5
    2109  0.0000000e+00 9.94e+04 8.31e+01  -1.0 6.17e+07    -  3.30e-04 3.20e-05f  4
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2110r 0.0000000e+00 9.94e+04 9.99e+02   0.6 0.00e+00    -  0.00e+00 4.78e-07R  6
    2111r 0.0000000e+00 9.29e+04 9.90e+02   0.6 7.45e+03    -  8.88e-03 9.91e-04f  1
    2112r 0.0000000e+00 9.29e+04 9.99e+02  -0.1 0.00e+00    -  0.00e+00 3.55e-07R  9
    2113r 0.0000000e+00 9.26e+04 9.90e+02  -0.1 2.99e+03    -  8.78e-03 9.90e-04f  1
    2114  0.0000000e+00 9.23e+04 1.24e+01  -1.0 5.03e+06    -  4.53e-03 2.30e-03f  2
    2115  0.0000000e+00 9.23e+04 1.63e+01  -1.0 7.19e+06    -  6.29e-03 7.06e-04h  1
    2116  0.0000000e+00 9.22e+04 2.32e+01  -1.0 8.21e+06    -  3.24e-03 4.87e-04f  5
    2117  0.0000000e+00 9.22e+04 1.38e+02  -1.0 7.58e+06    -  5.82e-03 1.85e-04h  4
    2118  0.0000000e+00 9.22e+04 3.82e+02  -1.0 6.35e+06    -  5.90e-03 2.41e-04f  5
    2119  0.0000000e+00 9.22e+04 1.17e+03  -1.0 6.27e+06    -  5.73e-03 9.20e-05h  6
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2120  0.0000000e+00 9.22e+04 4.93e+03  -1.0 6.05e+06    -  1.16e-02 2.92e-05h  8
    2121  0.0000000e+00 3.76e+05 2.66e+04  -1.0 4.91e+06    -  2.47e-03 4.03e-03f  2
    2122  0.0000000e+00 3.81e+05 2.61e+04  -1.0 6.65e+05    -  2.81e-02 3.31e-03h  4
    2123  0.0000000e+00 3.79e+05 2.47e+04  -1.0 1.36e+06    -  9.45e-03 5.38e-03h  2
    2124  0.0000000e+00 3.67e+05 4.09e+12  -1.0 2.75e-02  14.2 7.93e-03 1.00e+00h  1
    2125  0.0000000e+00 3.67e+05 4.22e+16  -1.0 2.42e+05  13.7 1.69e-05 1.15e-05h  2
    WARNING: Problem in step computation; switching to emergency mode.
    2126r 0.0000000e+00 3.67e+05 9.99e+02   1.9 0.00e+00  19.5 0.00e+00 0.00e+00R  1
    2127r 0.0000000e+00 5.48e+04 1.13e+03   1.9 1.19e+05    -  5.08e-03 1.02e-03f  1
    2128  0.0000000e+00 5.48e+04 3.38e+02  -1.0 5.50e+05    -  9.62e-02 1.64e-04h  8
    2129  0.0000000e+00 5.48e+04 3.75e+02  -1.0 4.86e+05    -  1.69e-01 1.47e-04h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2130  0.0000000e+00 5.48e+04 4.96e+02  -1.0 4.88e+05    -  2.19e-01 7.37e-05h 10
    2131  0.0000000e+00 5.48e+04 3.19e+03  -1.0 4.90e+05    -  3.36e-01 1.27e-04h  9
    2132  0.0000000e+00 5.56e+04 1.40e+04  -1.0 5.61e+06    -  5.34e-03 1.55e-02f  2
    2133  0.0000000e+00 5.56e+04 2.52e+04  -1.0 4.96e+05    -  2.62e-01 1.19e-02h  3
    2134  0.0000000e+00 5.54e+04 5.59e+04  -1.0 1.74e+05    -  1.99e-01 4.07e-03h  2
    2135  0.0000000e+00 5.51e+04 1.60e+05  -1.0 1.76e+05    -  2.01e-01 6.61e-03h  2
    2136  0.0000000e+00 5.50e+04 4.78e+05  -1.0 1.17e+05    -  1.97e-01 1.14e-03h  5
    2137  0.0000000e+00 5.50e+04 1.92e+06  -1.0 1.14e+05    -  2.02e-01 2.27e-04h  8
    2138  0.0000000e+00 5.31e+04 5.17e+07  -1.0 1.13e+05    -  1.96e-01 4.70e-02w  1
    2139  0.0000000e+00 8.13e+04 1.57e+07  -1.0 2.73e+04    -  3.70e-01 1.40e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2140  0.0000000e+00 8.18e+04 1.70e+07  -1.0 2.74e+04    -  6.31e-02 5.75e-02w  1
    2141  0.0000000e+00 5.50e+04 7.48e+06  -1.0 2.59e+04    -  1.96e-01 1.83e-04h  8
    2142  0.0000000e+00 5.35e+04 4.45e+08  -1.0 1.15e+05    -  1.92e-01 5.84e-02h  1
    2143  0.0000000e+00 5.70e+04 4.13e+08  -1.0 3.12e+04    -  5.13e-01 5.11e-02h  3
    2144  0.0000000e+00 5.88e+04 3.85e+08  -1.0 3.05e+04    -  1.66e-02 5.26e-02h  3
    2145  0.0000000e+00 6.00e+04 3.60e+08  -1.0 3.07e+04    -  1.21e-03 5.31e-02f  3
    2146  0.0000000e+00 6.92e+04 3.08e+08  -1.0 2.94e+04    -  7.25e-02 1.05e-01h  2
    2147  0.0000000e+00 7.27e+04 2.76e+08  -1.0 2.72e+04    -  1.59e-01 9.54e-02h  2
    2148  0.0000000e+00 7.60e+04 2.73e+08  -1.0 8.38e+05    -  2.86e-03 2.26e-02f  2
    2149  0.0000000e+00 8.90e+04 2.70e+08  -1.0 2.15e+06    -  9.92e-03 1.01e-02h  2
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2150  0.0000000e+00 8.91e+04 2.69e+08  -1.0 2.11e+06    -  1.96e-02 6.24e-04h  6
    2151  0.0000000e+00 8.91e+04 2.65e+08  -1.0 1.70e+06    -  3.26e-02 7.69e-04h  6
    2152  0.0000000e+00 8.91e+04 2.58e+08  -1.0 1.26e+06    -  7.35e-02 1.04e-03h  6
    2153  0.0000000e+00 1.56e+05 4.24e+08  -1.0 7.60e+05    -  2.51e-01 5.56e-02w  1
    2154  0.0000000e+00 1.55e+05 4.16e+09  -1.0 2.17e+05    -  8.16e-01 7.84e-03w  1
    2155  0.0000000e+00 1.55e+05 3.57e+10  -1.0 1.07e+05    -  9.85e-01 4.16e-04w  1
    2156  0.0000000e+00 8.96e+04 5.04e+08  -1.0 8.07e+04    -  2.51e-01 6.94e-03h  3
    2157  0.0000000e+00 8.98e+04 3.15e+09  -1.0 2.95e+05    -  7.04e-01 4.60e-03h  6
    2158  0.0000000e+00 8.96e+04 7.10e+09  -1.0 2.69e+05    -  2.11e-01 2.33e-03h  7
    2159  0.0000000e+00 8.96e+04 2.16e+10  -1.0 4.11e+05    -  3.58e-01 3.95e-04h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2160  0.0000000e+00 8.95e+04 3.45e+10  -1.0 5.73e+05    -  1.13e-01 2.78e-04h  9
    2161  0.0000000e+00 2.12e+05 3.23e+10  -1.0 6.20e+05    -  6.53e-02 6.53e-02s 18
    2162  0.0000000e+00 2.12e+05 3.22e+10  -1.0 4.64e+05    -  2.34e-03 2.34e-03s 18
    2163  0.0000000e+00 2.12e+05 3.22e+10  -1.0 1.62e+05    -  1.07e-04 1.07e-04s 18
    2164  0.0000000e+00 2.12e+05 3.22e+10  -1.0 7.57e+04    -  1.68e-04 1.68e-04s 18
    2165r 0.0000000e+00 2.12e+05 1.00e+03   1.7 0.00e+00    -  0.00e+00 0.00e+00R  1
    2166r 0.0000000e+00 6.39e+04 1.12e+03   1.7 5.52e+04    -  3.11e-03 1.00e-03f  1
    2167  0.0000000e+00 5.73e+04 4.66e+02  -1.0 1.37e+05    -  3.98e-01 1.19e-01h  1
    2168  0.0000000e+00 5.72e+04 9.90e+04  -1.0 4.75e+04    -  8.34e-01 2.22e-03h  1
    2169  0.0000000e+00 5.73e+04 3.02e+05  -1.0 3.48e+06    -  2.23e-02 1.52e-03h  3
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2170  0.0000000e+00 5.76e+04 1.60e+06  -1.0 3.54e+06    -  3.36e-02 2.43e-03h  2
    2171  0.0000000e+00 5.77e+04 1.60e+07  -1.0 3.63e+06    -  3.34e-02 1.34e-03h  2
    2172  0.0000000e+00 5.78e+04 9.39e+08  -1.0 3.67e+06    -  9.19e-02 1.42e-03h  1
    2173  0.0000000e+00 5.78e+04 1.35e+12  -1.0 3.47e+06    -  2.42e-02 1.61e-05h  1
    WARNING: Problem in step computation; switching to emergency mode.
    2174r 0.0000000e+00 5.78e+04 1.00e+03   1.1 0.00e+00    -  0.00e+00 0.00e+00R  1
    2175r 0.0000000e+00 5.76e+04 1.05e+03   1.1 5.40e+04    -  1.67e-03 9.60e-04f  1
    2176  0.0000000e+00 5.74e+04 1.47e+03  -1.0 3.19e+05    -  1.67e-01 4.02e-03h  1
    2177  0.0000000e+00 5.73e+04 4.22e+04  -1.0 3.54e+05    -  8.83e-03 3.79e-04h  1
    2178  0.0000000e+00 5.71e+04 3.16e+06  -1.0 4.37e+05    -  1.96e-01 4.07e-03h  1
    2179  0.0000000e+00 5.71e+04 8.60e+08  -1.0 3.89e+05    -  1.29e-02 4.78e-05h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2180r 0.0000000e+00 5.71e+04 9.99e+02   0.9 0.00e+00    -  0.00e+00 4.90e-07R  4
    2181r 0.0000000e+00 5.70e+04 1.00e+03   0.9 1.00e+04    -  7.72e-04 1.21e-03f  1
    2182  0.0000000e+00 5.67e+04 4.78e+02  -1.0 2.86e+05    -  4.62e-02 4.27e-03h  1
    2183  0.0000000e+00 5.67e+04 8.99e+05  -1.0 4.29e+05    -  2.06e-01 1.61e-04h  1
    2184  0.0000000e+00 5.66e+04 1.54e+07  -1.0 6.94e+05    -  4.33e-02 2.37e-03h  1
    2185  0.0000000e+00 5.66e+04 1.39e+11  -1.0 6.92e+05    -  2.11e-01 2.43e-05h  1
    WARNING: Problem in step computation; switching to emergency mode.
    2186r 0.0000000e+00 5.66e+04 9.99e+02   0.4 0.00e+00    -  0.00e+00 0.00e+00R  1
    2187r 0.0000000e+00 5.65e+04 9.96e+02   0.4 4.44e+04    -  5.30e-03 1.13e-03f  1
    2188  0.0000000e+00 5.64e+04 7.34e+02  -1.0 4.71e+05    -  2.10e-02 8.14e-04h  1
    2189  0.0000000e+00 5.64e+04 6.27e+03  -1.0 4.06e+05    -  1.06e-03 2.35e-04h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2190  0.0000000e+00 5.64e+04 7.33e+04  -1.0 4.41e+05    -  5.35e-03 1.00e-03f  1
    2191  0.0000000e+00 5.64e+04 1.03e+08  -1.0 2.61e+05    -  2.12e-02 1.84e-05h  1
    2192r 0.0000000e+00 5.64e+04 9.99e+02   0.2 0.00e+00    -  0.00e+00 3.85e-07R  2
    2193r 0.0000000e+00 5.63e+04 9.96e+02   0.2 1.42e+04    -  5.41e-03 1.01e-03f  1
    2194  0.0000000e+00 5.63e+04 8.13e+02  -1.0 5.78e+05    -  1.20e-02 3.55e-04f  1
    2195  0.0000000e+00 5.63e+04 1.89e+02  -1.0 4.17e+05    -  3.52e-04 4.24e-04h  1
    2196  0.0000000e+00 5.63e+04 6.03e+04  -1.0 4.31e+05    -  1.99e-03 2.51e-04h  1
    2197r 0.0000000e+00 5.63e+04 9.99e+02  -0.0 0.00e+00    -  0.00e+00 4.13e-07R  5
    2198r 0.0000000e+00 5.60e+04 1.00e+03  -0.0 1.42e+04    -  3.49e-04 1.80e-03f  1
    2199  0.0000000e+00 5.60e+04 1.20e+03  -1.0 5.59e+05    -  1.92e-02 3.23e-04h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2200  0.0000000e+00 5.60e+04 9.49e+05  -1.0 4.00e+05    -  8.60e-02 2.42e-04h  1
    2201  0.0000000e+00 5.60e+04 1.91e+08  -1.0 4.13e+05    -  8.00e-02 3.90e-04h  1
    2202r 0.0000000e+00 5.60e+04 9.99e+02  -0.3 0.00e+00    -  0.00e+00 2.51e-07R  5
    2203r 0.0000000e+00 5.53e+04 9.93e+02  -0.3 2.16e+04    -  2.57e-02 4.00e-03f  1
    2204  0.0000000e+00 5.53e+04 1.02e+03  -1.0 5.13e+05    -  2.07e-02 3.47e-04h  1
    2205  0.0000000e+00 5.53e+04 5.15e+05  -1.0 3.76e+05    -  6.55e-02 3.23e-04h  1
    2206  0.0000000e+00 5.52e+04 1.08e+08  -1.0 3.86e+05    -  8.01e-02 3.77e-04h  1
    2207r 0.0000000e+00 5.52e+04 9.99e+02  -0.5 0.00e+00    -  0.00e+00 2.69e-07R  5
    2208r 0.0000000e+00 5.53e+04 9.91e+02  -0.5 2.69e+04    -  2.15e-02 5.45e-03f  1
    2209  0.0000000e+00 5.53e+04 1.19e+03  -1.0 5.74e+05    -  1.82e-02 2.55e-04h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2210  0.0000000e+00 5.53e+04 1.03e+06  -1.0 3.94e+05    -  9.51e-02 2.57e-04h  1
    2211  0.0000000e+00 5.53e+04 8.49e+08  -1.0 4.04e+05    -  2.36e-01 2.83e-04h  1
    2212  0.0000000e+00 5.53e+04 1.61e+12  -1.0 3.46e+05    -  3.59e-02 1.88e-05h  1
    2213  0.0000000e+00 5.51e+04 5.33e+14  -1.0 3.36e+05    -  9.75e-01 2.90e-03h  1
    2214  0.0000000e+00 5.51e+04 5.22e+14  -1.0 3.37e+05    -  2.90e-02 8.03e-04h 11
    WARNING: Problem in step computation; switching to emergency mode.
    2215r 0.0000000e+00 5.51e+04 1.00e+03  -1.0 0.00e+00    -  0.00e+00 0.00e+00R  1
    2216r 0.0000000e+00 5.53e+04 9.84e+02  -1.0 4.92e+03    -  1.02e-01 7.34e-03f  1
    2217r 0.0000000e+00 5.56e+04 8.72e+02  -1.0 5.07e+03    -  9.09e-02 1.14e-01f  1
    2218r 0.0000000e+00 5.53e+04 1.13e+03  -1.0 4.06e+03    -  5.70e-01 1.54e-01f  1
    2219r 0.0000000e+00 5.29e+04 1.02e+03  -1.0 3.93e+03    -  2.27e-01 1.50e-01f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2220r 0.0000000e+00 5.24e+04 4.34e+02  -1.0 2.96e+03    -  2.53e-01 8.53e-01f  1
    2221  0.0000000e+00 4.05e+04 1.35e+03  -1.0 1.50e+04    -  6.05e-01 1.94e-01f  2
    2222  0.0000000e+00 3.62e+04 2.29e+03  -1.0 1.04e+05    -  4.54e-01 1.07e-01f  3
    2223  0.0000000e+00 3.52e+04 4.19e+03  -1.0 2.73e+04    -  5.70e-01 2.72e-02h  5
    2224  0.0000000e+00 3.47e+04 6.28e+03  -1.0 3.64e+04    -  6.21e-01 1.37e-02h  6
    2225  0.0000000e+00 3.38e+04 9.45e+03  -1.0 5.64e+04    -  9.89e-01 2.74e-02h  5
    2226  0.0000000e+00 3.38e+04 1.22e+04  -1.0 5.51e+04    -  8.26e-01 1.08e-04h 13
    2227  0.0000000e+00 3.37e+04 1.53e+04  -1.0 6.29e+04    -  9.96e-01 8.61e-04h 10
    2228  0.0000000e+00 3.37e+04 3.11e+04  -1.0 6.28e+04    -  9.34e-01 5.38e-05h 14
    2229  0.0000000e+00 3.37e+04 5.19e+04  -1.0 6.60e+04    -  5.26e-01 4.30e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2230  0.0000000e+00 3.37e+04 8.51e+04  -1.0 6.43e+04    -  5.09e-01 1.08e-04h 13
    2231  0.0000000e+00 5.57e+05 2.03e+05  -1.0 6.32e+04    -  5.13e-01 4.41e-01w  1
    2232  0.0000000e+00 6.09e+05 1.18e+05  -1.0 6.95e+04    -  5.75e-01 3.38e-01w  1
    2233  0.0000000e+00 3.82e+05 8.03e+04  -1.0 1.33e+04    -  1.00e+00 4.97e-01w  1
    2234r 0.0000000e+00 3.37e+04 1.00e+03  -1.0 0.00e+00    -  0.00e+00 8.21e-10R 29
    2235r 0.0000000e+00 3.37e+04 1.01e+03  -1.0 4.62e+03    -  1.03e-01 1.15e-03f  1
    2236  0.0000000e+00 2.97e+04 2.70e+03  -1.0 1.14e+04    -  7.53e-01 1.10e-01h  3
    2237  0.0000000e+00 2.88e+04 6.86e+03  -1.0 3.40e+04    -  9.89e-01 2.80e-02h  5
    2238  0.0000000e+00 2.88e+04 1.07e+04  -1.0 4.39e+04    -  9.12e-01 1.10e-04h 13
    2239  0.0000000e+00 2.88e+04 1.47e+04  -1.0 5.42e+04    -  9.92e-01 8.79e-04h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2240  0.0000000e+00 2.88e+04 1.85e+04  -1.0 5.35e+04    -  1.00e+00 2.75e-05h 15
    2241  0.0000000e+00 2.88e+04 3.19e+04  -1.0 5.74e+04    -  1.00e+00 4.39e-04h 11
    2242  0.0000000e+00 2.88e+04 7.08e+04  -1.0 5.71e+04    -  1.00e+00 2.75e-05h 15
    2243  0.0000000e+00 2.88e+04 1.56e+05  -1.0 5.70e+04    -  1.00e+00 1.10e-04h 13
    2244  0.0000000e+00 2.88e+04 3.45e+05  -1.0 5.29e+04    -  1.00e+00 5.49e-05h 14
    2245  0.0000000e+00 2.88e+04 7.59e+05  -1.0 4.47e+04    -  1.00e+00 1.10e-04h 13
    2246  0.0000000e+00 1.93e+05 1.14e+06  -1.0 2.88e+04    -  1.00e+00 4.50e-01w  1
    2247  0.0000000e+00 6.38e+05 9.34e+05  -1.0 5.05e+05    -  1.33e-01 6.43e-02w  1
    2248  0.0000000e+00 4.51e+05 8.27e+05  -1.0 1.49e+04    -  1.00e+00 4.93e-01w  1
    2249  0.0000000e+00 2.88e+04 1.67e+06  -1.0 6.24e+03    -  1.00e+00 1.10e-04h 12
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2250  0.0000000e+00 2.88e+04 3.67e+06  -1.0 8.48e+03    -  1.00e+00 8.80e-04h 10
    2251  0.0000000e+00 2.80e+04 7.83e+06  -1.0 1.45e+04    -  1.00e+00 2.73e-02h  5
    2252  0.0000000e+00 2.72e+04 1.70e+07  -1.0 2.54e+04    -  1.00e+00 2.52e-02h  5
    2253  0.0000000e+00 2.69e+04 3.74e+07  -1.0 1.98e+04    -  1.00e+00 1.30e-02h  6
    2254  0.0000000e+00 2.68e+04 8.26e+07  -1.0 1.36e+04    -  1.00e+00 3.41e-03h  8
    2255  0.0000000e+00 2.67e+04 1.82e+08  -1.0 1.02e+04    -  1.00e+00 1.75e-03h  9
    2256  0.0000000e+00 2.67e+04 3.98e+08  -1.0 9.38e+03    -  1.00e+00 1.77e-03h  9
    2257  0.0000000e+00 2.66e+04 8.68e+08  -1.0 9.07e+03    -  1.00e+00 1.77e-03h  9
    2258  0.0000000e+00 2.66e+04 1.89e+09  -1.0 8.80e+03    -  1.00e+00 1.77e-03h  9
    2259  0.0000000e+00 1.65e+04 2.39e+09  -1.0 8.59e+03    -  1.00e+00 4.54e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2260  0.0000000e+00 1.56e+04 3.16e+09  -1.0 1.98e+04    -  1.00e+00 7.05e-02w  1
    2261  0.0000000e+00 1.68e+04 2.70e+09  -1.0 3.86e+03    -  1.00e+00 4.95e-01w  1
    2262  0.0000000e+00 2.65e+04 4.13e+09  -1.0 1.68e+03    -  1.00e+00 1.77e-03h  8
    2263  0.0000000e+00 2.65e+04 8.99e+09  -1.0 8.43e+03    -  1.00e+00 1.77e-03h  9
    2264  0.0000000e+00 2.64e+04 1.96e+10  -1.0 8.32e+03    -  1.00e+00 1.78e-03h  9
    2265  0.0000000e+00 2.64e+04 4.27e+10  -1.0 8.23e+03    -  1.00e+00 1.78e-03h  9
    2266  0.0000000e+00 2.63e+04 9.29e+10  -1.0 8.17e+03    -  1.00e+00 1.78e-03h  9
    2267  0.0000000e+00 2.63e+04 2.02e+11  -1.0 8.13e+03    -  1.00e+00 1.78e-03h  9
    2268  0.0000000e+00 2.63e+04 4.40e+11  -1.0 8.09e+03    -  1.00e+00 1.78e-03h  9
    2269  0.0000000e+00 2.62e+04 9.58e+11  -1.0 8.07e+03    -  1.00e+00 1.78e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2270  0.0000000e+00 2.62e+04 2.09e+12  -1.0 8.05e+03    -  1.00e+00 1.78e-03h  9
    2271  0.0000000e+00 2.61e+04 4.54e+12  -1.0 8.04e+03    -  1.00e+00 1.78e-03h  9
    2272  0.0000000e+00 1.66e+04 5.12e+12  -1.0 8.03e+03    -  1.00e+00 4.55e-01w  1
    2273  0.0000000e+00 1.59e+04 6.69e+12  -1.0 2.32e+04    -  1.00e+00 7.35e-02w  1
    2274  0.0000000e+00 1.80e+04 5.71e+12  -1.0 4.08e+03    -  1.00e+00 4.95e-01w  1
    2275  0.0000000e+00 2.61e+04 6.15e+12  -1.0 1.57e+03    -  1.00e+00 1.78e-03h  8
    2276  0.0000000e+00 2.61e+04 6.15e+12  -1.0 2.42e+04    -  1.00e+00 1.11e-04h 13
    2277  0.0000000e+00 2.61e+04 6.15e+12  -1.0 5.02e+04    -  9.08e-01 5.55e-05h 14
    2278  0.0000000e+00 2.61e+04 6.15e+12  -1.0 5.44e+04    -  1.00e+00 5.55e-05h 14
    2279  0.0000000e+00 2.61e+04 6.15e+12  -1.0 5.47e+04    -  1.00e+00 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2280  0.0000000e+00 2.61e+04 6.15e+12  -1.0 5.47e+04    -  1.00e+00 5.55e-05h 14
    2281  0.0000000e+00 2.61e+04 6.15e+12  -1.0 5.47e+04    -  1.00e+00 5.55e-05h 14
    2282  0.0000000e+00 2.61e+04 6.15e+12  -1.0 5.47e+04    -  1.00e+00 5.55e-05h 14
    2283  0.0000000e+00 2.61e+04 6.15e+12  -1.0 5.47e+04    -  1.00e+00 5.55e-05h 14
    2284  0.0000000e+00 2.61e+04 6.15e+12  -1.0 5.47e+04    -  1.00e+00 5.55e-05h 14
    2285  0.0000000e+00 4.28e+05 9.71e+12  -1.0 5.47e+04    -  1.00e+00 4.54e-01w  1
    2286  0.0000000e+00 8.16e+05 9.45e+12  -1.0 3.16e+07    -  1.92e-03 1.04e-03w  1
    2287  0.0000000e+00 8.00e+05 1.33e+13  -1.0 1.44e+05    -  1.00e+00 5.30e-02w  1
    2288  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.48e+04    -  1.00e+00 5.55e-05h 13
    2289  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.47e+04    -  1.00e+00 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2290  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.47e+04    -  1.00e+00 5.55e-05h 14
    2291  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.47e+04    -  8.01e-01 5.55e-05h 14
    2292  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.47e+04    -  6.61e-01 5.55e-05h 14
    2293  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.47e+04    -  6.62e-01 5.55e-05h 14
    2294  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.47e+04    -  6.62e-01 5.55e-05h 14
    2295  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.47e+04    -  6.63e-01 5.55e-05h 14
    2296  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.47e+04    -  6.62e-01 5.55e-05h 14
    2297  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.47e+04    -  6.63e-01 5.55e-05h 14
    2298  0.0000000e+00 4.27e+05 1.22e+13  -1.0 5.47e+04    -  6.63e-01 4.54e-01w  1
    2299  0.0000000e+00 8.14e+05 1.20e+13  -1.0 2.99e+07    -  2.03e-03 1.09e-03w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2300  0.0000000e+00 7.99e+05 9.90e+12  -1.0 1.45e+05    -  1.00e+00 5.25e-02w  1
    2301  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.61e+04    -  6.63e-01 5.55e-05h 13
    2302  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.47e+04    -  6.63e-01 5.55e-05h 14
    2303  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.47e+04    -  6.63e-01 5.55e-05h 14
    2304  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.47e+04    -  6.63e-01 5.55e-05h 14
    2305  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.47e+04    -  6.63e-01 5.55e-05h 14
    2306  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.47e+04    -  6.63e-01 5.55e-05h 14
    2307  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.47e+04    -  6.63e-01 5.55e-05h 14
    2308  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.47e+04    -  6.63e-01 5.55e-05h 14
    2309  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.47e+04    -  6.63e-01 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2310  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.47e+04    -  6.63e-01 5.55e-05h 14
    2311  0.0000000e+00 4.26e+05 1.22e+13  -1.0 5.47e+04    -  6.63e-01 4.54e-01w  1
    2312  0.0000000e+00 8.13e+05 1.20e+13  -1.0 2.66e+07    -  2.28e-03 1.23e-03w  1
    2313  0.0000000e+00 7.98e+05 9.91e+12  -1.0 1.44e+05    -  1.00e+00 5.29e-02w  1
    2314  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.50e+04    -  6.63e-01 5.55e-05h 13
    2315  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.47e+04    -  6.64e-01 5.55e-05h 14
    2316  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.47e+04    -  6.64e-01 5.55e-05h 14
    2317  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.47e+04    -  6.64e-01 5.55e-05h 14
    2318  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.47e+04    -  6.64e-01 5.55e-05h 14
    2319  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.47e+04    -  6.64e-01 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2320  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.64e-01 5.55e-05h 14
    2321  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.47e+04    -  6.64e-01 5.55e-05h 14
    2322  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.64e-01 5.55e-05h 14
    2323  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.64e-01 5.55e-05h 14
    2324  0.0000000e+00 4.26e+05 1.22e+13  -1.0 5.46e+04    -  6.64e-01 4.54e-01w  1
    2325  0.0000000e+00 8.12e+05 1.20e+13  -1.0 2.40e+07    -  2.53e-03 1.36e-03w  1
    2326  0.0000000e+00 7.96e+05 9.92e+12  -1.0 1.44e+05    -  1.00e+00 5.34e-02w  1
    2327  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.40e+04    -  6.64e-01 5.55e-05h 13
    2328  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.64e-01 5.55e-05h 14
    2329  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.64e-01 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2330  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.64e-01 5.55e-05h 14
    2331  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.64e-01 5.55e-05h 14
    2332  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.65e-01 5.55e-05h 14
    2333  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.65e-01 5.55e-05h 14
    2334  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.65e-01 5.55e-05h 14
    2335  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.65e-01 5.55e-05h 14
    2336  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.65e-01 5.55e-05h 14
    2337  0.0000000e+00 4.25e+05 1.22e+13  -1.0 5.46e+04    -  6.65e-01 4.54e-01w  1
    2338  0.0000000e+00 8.11e+05 1.20e+13  -1.0 2.18e+07    -  2.78e-03 1.50e-03w  1
    2339  0.0000000e+00 7.95e+05 9.92e+12  -1.0 1.43e+05    -  1.00e+00 5.38e-02w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2340  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.30e+04    -  6.65e-01 5.55e-05h 13
    2341  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.65e-01 5.55e-05h 14
    2342  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.65e-01 5.55e-05h 14
    2343  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.65e-01 5.55e-05h 14
    2344  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.65e-01 5.55e-05h 14
    2345  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.65e-01 5.55e-05h 14
    2346  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.65e-01 5.55e-05h 14
    2347  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.66e-01 5.55e-05h 14
    2348  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.65e-01 5.55e-05h 14
    2349  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.66e-01 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2350  0.0000000e+00 4.25e+05 1.22e+13  -1.0 5.46e+04    -  6.66e-01 4.55e-01w  1
    2351  0.0000000e+00 8.10e+05 1.19e+13  -1.0 2.00e+07    -  3.02e-03 1.63e-03w  1
    2352  0.0000000e+00 7.94e+05 9.93e+12  -1.0 1.43e+05    -  1.00e+00 5.42e-02w  1
    2353  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.21e+04    -  6.66e-01 5.55e-05h 13
    2354  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.66e-01 5.55e-05h 14
    2355  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.66e-01 5.55e-05h 14
    2356  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.45e+04    -  6.66e-01 5.55e-05h 14
    2357  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.66e-01 5.55e-05h 14
    2358  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.66e-01 5.55e-05h 14
    2359  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.66e-01 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2360  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.66e-01 5.55e-05h 14
    2361  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.66e-01 5.55e-05h 14
    2362  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.66e-01 5.55e-05h 14
    2363  0.0000000e+00 4.24e+05 1.22e+13  -1.0 5.45e+04    -  6.66e-01 4.55e-01w  1
    2364  0.0000000e+00 8.09e+05 1.19e+13  -1.0 1.85e+07    -  3.28e-03 1.77e-03w  1
    2365  0.0000000e+00 7.93e+05 9.94e+12  -1.0 1.42e+05    -  1.00e+00 5.46e-02w  1
    2366  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.12e+04    -  6.66e-01 5.55e-05h 13
    2367  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.46e+04    -  6.67e-01 5.55e-05h 14
    2368  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.45e+04    -  6.66e-01 5.55e-05h 14
    2369  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.45e+04    -  6.67e-01 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2370  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.45e+04    -  6.67e-01 5.55e-05h 14
    2371  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.45e+04    -  6.67e-01 5.55e-05h 14
    2372  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.45e+04    -  6.67e-01 5.55e-05h 14
    2373  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.45e+04    -  6.67e-01 5.55e-05h 14
    2374  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.45e+04    -  6.67e-01 5.55e-05h 14
    2375  0.0000000e+00 2.60e+04 6.15e+12  -1.0 5.45e+04    -  6.67e-01 5.55e-05h 14
    2376  0.0000000e+00 4.24e+05 1.22e+13  -1.0 5.45e+04    -  6.67e-01 4.55e-01w  1
    2377  0.0000000e+00 8.08e+05 1.19e+13  -1.0 1.71e+07    -  3.53e-03 1.91e-03w  1
    2378  0.0000000e+00 7.92e+05 9.94e+12  -1.0 1.42e+05    -  1.00e+00 5.51e-02w  1
    2379  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.03e+04    -  6.67e-01 5.55e-05h 13
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2380  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.45e+04    -  6.67e-01 5.55e-05h 14
    2381  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.45e+04    -  6.67e-01 5.55e-05h 14
    2382  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.45e+04    -  6.67e-01 5.55e-05h 14
    2383  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.45e+04    -  6.67e-01 5.55e-05h 14
    2384  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.45e+04    -  6.68e-01 5.55e-05h 14
    2385  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.45e+04    -  6.68e-01 5.55e-05h 14
    2386  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.45e+04    -  6.68e-01 5.55e-05h 14
    2387  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.45e+04    -  6.68e-01 5.55e-05h 14
    2388  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.45e+04    -  6.68e-01 5.55e-05h 14
    2389  0.0000000e+00 4.23e+05 1.22e+13  -1.0 5.45e+04    -  6.68e-01 4.55e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2390  0.0000000e+00 8.07e+05 1.19e+13  -1.0 1.59e+07    -  3.80e-03 2.05e-03w  1
    2391  0.0000000e+00 7.91e+05 9.95e+12  -1.0 1.41e+05    -  1.00e+00 5.56e-02w  1
    2392  0.0000000e+00 2.59e+04 6.15e+12  -1.0 4.94e+04    -  6.68e-01 5.55e-05h 13
    2393  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.45e+04    -  6.68e-01 5.55e-05h 14
    2394  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.45e+04    -  6.68e-01 5.55e-05h 14
    2395  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.45e+04    -  6.68e-01 5.55e-05h 14
    2396  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.44e+04    -  6.68e-01 5.55e-05h 14
    2397  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.45e+04    -  6.68e-01 5.55e-05h 14
    2398  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.44e+04    -  6.68e-01 5.55e-05h 14
    2399  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.45e+04    -  6.68e-01 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2400  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.45e+04    -  6.68e-01 5.55e-05h 14
    2401  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.45e+04    -  6.69e-01 5.55e-05h 14
    2402  0.0000000e+00 4.23e+05 1.22e+13  -1.0 5.45e+04    -  6.69e-01 4.55e-01w  1
    2403  0.0000000e+00 8.06e+05 1.19e+13  -1.0 1.47e+07    -  4.11e-03 2.22e-03w  1
    2404  0.0000000e+00 7.89e+05 9.95e+12  -1.0 1.40e+05    -  1.00e+00 5.62e-02w  1
    2405  0.0000000e+00 2.59e+04 6.15e+12  -1.0 4.84e+04    -  6.69e-01 5.55e-05h 13
    2406  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.44e+04    -  6.69e-01 5.55e-05h 14
    2407  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.45e+04    -  6.69e-01 5.55e-05h 14
    2408  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.44e+04    -  6.69e-01 5.55e-05h 14
    2409  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.44e+04    -  6.69e-01 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2410  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.44e+04    -  6.69e-01 5.55e-05h 14
    2411  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.44e+04    -  6.69e-01 5.55e-05h 14
    2412  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.44e+04    -  6.69e-01 5.55e-05h 14
    2413  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.44e+04    -  6.69e-01 5.55e-05h 14
    2414  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.44e+04    -  6.69e-01 5.55e-05h 14
    2415  0.0000000e+00 4.22e+05 1.22e+13  -1.0 5.44e+04    -  6.69e-01 4.55e-01w  1
    2416  0.0000000e+00 8.05e+05 1.19e+13  -1.0 1.33e+07    -  4.53e-03 2.44e-03w  1
    2417  0.0000000e+00 7.88e+05 9.95e+12  -1.0 1.39e+05    -  1.00e+00 5.70e-02w  1
    2418  0.0000000e+00 2.59e+04 6.15e+12  -1.0 4.69e+04    -  6.69e-01 5.55e-05h 13
    2419  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.44e+04    -  6.70e-01 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2420  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.44e+04    -  6.70e-01 5.55e-05h 14
    2421  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.44e+04    -  6.70e-01 5.55e-05h 14
    2422  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.44e+04    -  6.70e-01 5.55e-05h 14
    2423  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.44e+04    -  6.70e-01 5.55e-05h 14
    2424  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.44e+04    -  6.70e-01 5.55e-05h 14
    2425  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.44e+04    -  6.70e-01 5.55e-05h 14
    2426  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.44e+04    -  6.70e-01 5.55e-05h 14
    2427  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.44e+04    -  6.70e-01 5.55e-05h 14
    2428  0.0000000e+00 4.21e+05 1.22e+13  -1.0 5.44e+04    -  6.70e-01 4.55e-01w  1
    2429  0.0000000e+00 8.03e+05 1.19e+13  -1.0 1.15e+07    -  5.25e-03 2.83e-03w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2430  0.0000000e+00 7.86e+05 9.95e+12  -1.0 1.37e+05    -  1.00e+00 5.86e-02w  1
    2431  0.0000000e+00 2.59e+04 6.15e+12  -1.0 4.45e+04    -  6.70e-01 5.55e-05h 13
    2432  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.44e+04    -  6.71e-01 5.55e-05h 14
    2433  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.43e+04    -  6.71e-01 5.55e-05h 14
    2434  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.43e+04    -  6.71e-01 5.55e-05h 14
    2435  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.43e+04    -  6.71e-01 5.55e-05h 14
    2436  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.43e+04    -  6.71e-01 5.55e-05h 14
    2437  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.43e+04    -  6.71e-01 5.55e-05h 14
    2438  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.43e+04    -  6.71e-01 5.55e-05h 14
    2439  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.42e+04    -  6.71e-01 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2440  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.41e+04    -  6.72e-01 5.55e-05h 14
    2441  0.0000000e+00 4.20e+05 1.22e+13  -1.0 5.42e+04    -  6.72e-01 4.55e-01w  1
    2442  0.0000000e+00 8.01e+05 1.19e+13  -1.0 8.87e+06    -  6.81e-03 3.67e-03w  1
    2443  0.0000000e+00 7.82e+05 9.91e+12  -1.0 1.33e+05    -  1.00e+00 6.26e-02w  1
    2444  0.0000000e+00 2.59e+04 6.15e+12  -1.0 4.00e+04    -  6.72e-01 5.55e-05h 13
    2445  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.42e+04    -  6.72e-01 5.55e-05h 14
    2446  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.42e+04    -  6.72e-01 5.55e-05h 14
    2447  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.42e+04    -  6.72e-01 5.55e-05h 14
    2448  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.41e+04    -  6.72e-01 5.55e-05h 14
    2449  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.41e+04    -  6.73e-01 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2450  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.41e+04    -  6.73e-01 5.55e-05h 14
    2451  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.41e+04    -  6.73e-01 5.55e-05h 14
    2452  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.40e+04    -  6.73e-01 5.55e-05h 14
    2453  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.40e+04    -  6.74e-01 5.55e-05h 14
    2454  0.0000000e+00 4.16e+05 1.22e+13  -1.0 5.39e+04    -  6.74e-01 4.55e-01w  1
    2455  0.0000000e+00 7.96e+05 1.18e+13  -1.0 5.65e+06    -  1.07e-02 5.75e-03w  1
    2456  0.0000000e+00 7.74e+05 9.78e+12  -1.0 1.22e+05    -  1.00e+00 7.40e-02w  1
    2457  0.0000000e+00 2.59e+04 6.15e+12  -1.0 3.28e+04    -  6.74e-01 5.55e-05h 13
    2458  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.39e+04    -  6.74e-01 5.55e-05h 14
    2459  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.39e+04    -  6.74e-01 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2460  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.37e+04    -  6.75e-01 5.55e-05h 14
    2461  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.37e+04    -  6.75e-01 5.55e-05h 14
    2462  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.37e+04    -  6.76e-01 5.55e-05h 14
    2463  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.36e+04    -  6.76e-01 5.55e-05h 14
    2464  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.35e+04    -  6.76e-01 5.55e-05h 14
    2465  0.0000000e+00 2.59e+04 6.15e+12  -1.0 5.35e+04    -  6.77e-01 5.55e-05h 14
    2466  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.34e+04    -  6.77e-01 5.55e-05h 14
    2467  0.0000000e+00 4.10e+05 1.21e+13  -1.0 5.33e+04    -  6.78e-01 4.55e-01w  1
    2468  0.0000000e+00 7.85e+05 1.17e+13  -1.0 2.88e+06    -  2.10e-02 1.13e-02w  1
    2469  0.0000000e+00 7.53e+05 9.23e+12  -1.0 9.98e+04    -  1.00e+00 1.23e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2470  0.0000000e+00 2.58e+04 6.15e+12  -1.0 2.28e+04    -  6.78e-01 5.55e-05h 13
    2471  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.32e+04    -  6.78e-01 5.55e-05h 14
    2472  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.30e+04    -  6.79e-01 5.55e-05h 14
    2473  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.30e+04    -  6.80e-01 5.55e-05h 14
    2474  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.28e+04    -  6.80e-01 5.55e-05h 14
    2475  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.27e+04    -  6.81e-01 5.55e-05h 14
    2476  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.26e+04    -  6.82e-01 5.55e-05h 14
    2477  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.25e+04    -  6.83e-01 5.55e-05h 14
    2478  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.23e+04    -  6.84e-01 5.55e-05h 14
    2479  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.21e+04    -  6.85e-01 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2480  0.0000000e+00 3.96e+05 1.21e+13  -1.0 5.20e+04    -  6.85e-01 4.55e-01w  1
    2481  0.0000000e+00 7.58e+05 1.14e+13  -1.0 1.23e+06    -  4.90e-02 2.61e-02w  1
    2482  0.0000000e+00 6.61e+05 8.34e+12  -1.0 6.16e+04    -  1.00e+00 2.04e-01w  1
    2483  0.0000000e+00 2.58e+04 6.15e+12  -1.0 1.54e+04    -  6.85e-01 5.55e-05h 13
    2484  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.18e+04    -  6.86e-01 5.55e-05h 14
    2485  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.17e+04    -  6.87e-01 5.55e-05h 14
    2486  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.15e+04    -  6.88e-01 5.55e-05h 14
    2487  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.13e+04    -  6.89e-01 5.55e-05h 14
    2488  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.11e+04    -  6.91e-01 5.55e-05h 14
    2489  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.09e+04    -  6.92e-01 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2490  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.08e+04    -  6.93e-01 5.55e-05h 14
    2491  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.06e+04    -  6.94e-01 5.55e-05h 14
    2492  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.04e+04    -  6.95e-01 5.55e-05h 14
    2493  0.0000000e+00 3.78e+05 1.20e+13  -1.0 5.02e+04    -  6.96e-01 4.55e-01w  1
    2494  0.0000000e+00 6.99e+05 1.05e+13  -1.0 4.68e+05    -  1.26e-01 6.59e-02w  1
    2495  0.0000000e+00 4.83e+05 7.17e+12  -1.0 2.03e+04    -  1.00e+00 3.14e-01w  1
    2496  0.0000000e+00 2.58e+04 6.15e+12  -1.0 8.80e+03    -  6.96e-01 5.55e-05h 13
    2497  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.00e+04    -  6.97e-01 5.55e-05h 14
    2498  0.0000000e+00 2.58e+04 6.15e+12  -1.0 4.98e+04    -  6.98e-01 5.55e-05h 14
    2499  0.0000000e+00 2.58e+04 6.15e+12  -1.0 4.96e+04    -  6.99e-01 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2500  0.0000000e+00 2.58e+04 6.15e+12  -1.0 4.94e+04    -  7.00e-01 5.55e-05h 14
    2501  0.0000000e+00 2.58e+04 6.15e+12  -1.0 4.93e+04    -  7.02e-01 5.55e-05h 14
    2502  0.0000000e+00 2.58e+04 6.15e+12  -1.0 4.90e+04    -  7.02e-01 5.55e-05h 14
    2503  0.0000000e+00 2.58e+04 6.15e+12  -1.0 4.88e+04    -  7.04e-01 5.55e-05h 14
    2504  0.0000000e+00 2.58e+04 6.15e+12  -1.0 4.88e+04    -  7.04e-01 5.55e-05h 14
    2505  0.0000000e+00 2.58e+04 6.15e+12  -1.0 4.86e+04    -  7.05e-01 5.55e-05h 14
    2506  0.0000000e+00 3.63e+05 1.21e+13  -1.0 4.85e+04    -  7.06e-01 4.55e-01w  1
    2507  0.0000000e+00 5.65e+05 8.45e+12  -1.0 1.63e+05    -  3.34e-01 1.69e-01w  1
    2508  0.0000000e+00 3.13e+05 5.27e+12  -1.0 7.18e+03    -  1.00e+00 4.93e-01w  1
    2509  0.0000000e+00 2.58e+04 6.15e+12  -1.0 2.07e+03    -  7.06e-01 5.55e-05h 13
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2510  0.0000000e+00 2.58e+04 6.15e+12  -1.0 4.83e+04    -  7.07e-01 5.55e-05h 14
    2511  0.0000000e+00 2.58e+04 6.15e+12  -1.0 4.82e+04    -  7.08e-01 5.55e-05h 14
    2512  0.0000000e+00 2.58e+04 6.15e+12  -1.0 4.81e+04    -  7.09e-01 5.55e-05h 14
    2513  0.0000000e+00 2.58e+04 6.15e+12  -1.0 4.80e+04    -  7.09e-01 5.55e-05h 14
    2514  0.0000000e+00 2.58e+04 6.15e+12  -1.0 4.86e+04    -  7.05e-01 5.55e-05h 14
    2515  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.22e+04    -  6.81e-01 5.55e-05h 14
    2516  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.35e+04    -  6.78e-01 5.55e-05h 14
    2517  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.39e+04    -  6.76e-01 5.55e-05h 14
    2518  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.40e+04    -  6.75e-01 5.55e-05h 14
    2519  0.0000000e+00 4.20e+05 1.25e+13  -1.0 5.42e+04    -  6.75e-01 4.55e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2520  0.0000000e+00 5.00e+05 7.45e+12  -1.0 9.43e+04    -  4.97e-01 2.58e-01w  1
    2521  0.0000000e+00 2.76e+05 4.88e+12  -1.0 6.77e+03    -  1.00e+00 4.94e-01w  1
    2522  0.0000000e+00 2.58e+04 6.15e+12  -1.0 2.01e+03    -  6.75e-01 5.55e-05h 13
    2523  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.42e+04    -  6.75e-01 5.55e-05h 14
    2524  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.41e+04    -  6.75e-01 5.55e-05h 14
    2525  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.41e+04    -  6.75e-01 5.55e-05h 14
    2526  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.42e+04    -  6.75e-01 5.55e-05h 14
    2527  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.42e+04    -  6.75e-01 5.55e-05h 14
    2528  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.41e+04    -  6.75e-01 5.55e-05h 14
    2529  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.41e+04    -  6.76e-01 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2530  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.42e+04    -  6.76e-01 5.55e-05h 14
    2531  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.41e+04    -  6.76e-01 5.55e-05h 14
    2532  0.0000000e+00 4.18e+05 1.25e+13  -1.0 5.42e+04    -  6.76e-01 4.55e-01w  1
    2533  0.0000000e+00 4.98e+05 7.51e+12  -1.0 9.42e+04    -  4.92e-01 2.58e-01w  1
    2534  0.0000000e+00 2.74e+05 4.87e+12  -1.0 6.67e+03    -  1.00e+00 4.94e-01w  1
    2535  0.0000000e+00 2.58e+04 6.15e+12  -1.0 2.04e+03    -  6.76e-01 5.55e-05h 13
    2536  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.41e+04    -  6.76e-01 5.55e-05h 14
    2537  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.42e+04    -  6.76e-01 5.55e-05h 14
    2538  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.38e+04    -  6.76e-01 5.55e-05h 14
    2539  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.40e+04    -  6.76e-01 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2540  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.41e+04    -  6.76e-01 5.55e-05h 14
    2541  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.42e+04    -  6.76e-01 5.55e-05h 14
    2542  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.40e+04    -  6.76e-01 5.55e-05h 14
    2543  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.41e+04    -  6.76e-01 5.55e-05h 14
    2544  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.40e+04    -  6.77e-01 5.55e-05h 14
    2545  0.0000000e+00 4.17e+05 1.25e+13  -1.0 5.41e+04    -  6.77e-01 4.55e-01w  1
    2546  0.0000000e+00 4.97e+05 7.50e+12  -1.0 9.40e+04    -  4.93e-01 2.58e-01w  1
    2547  0.0000000e+00 2.74e+05 4.88e+12  -1.0 6.66e+03    -  1.00e+00 4.94e-01w  1
    2548  0.0000000e+00 2.58e+04 6.15e+12  -1.0 2.04e+03    -  6.77e-01 5.55e-05h 13
    2549  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.40e+04    -  6.77e-01 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2550  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.41e+04    -  6.77e-01 5.55e-05h 14
    2551  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.41e+04    -  6.77e-01 5.55e-05h 14
    2552  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.40e+04    -  6.77e-01 5.55e-05h 14
    2553  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.41e+04    -  6.77e-01 5.55e-05h 14
    2554  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.41e+04    -  6.77e-01 5.55e-05h 14
    2555  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.41e+04    -  6.77e-01 5.55e-05h 14
    2556  0.0000000e+00 2.58e+04 6.15e+12  -1.0 5.41e+04    -  6.77e-01 5.55e-05h 14
    2557  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.41e+04    -  6.77e-01 5.55e-05h 14
    2558  0.0000000e+00 4.17e+05 1.25e+13  -1.0 5.41e+04    -  6.77e-01 4.55e-01w  1
    2559  0.0000000e+00 4.96e+05 7.50e+12  -1.0 9.39e+04    -  4.93e-01 2.58e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2560  0.0000000e+00 2.73e+05 4.88e+12  -1.0 6.64e+03    -  1.00e+00 4.94e-01w  1
    2561  0.0000000e+00 2.57e+04 6.15e+12  -1.0 2.03e+03    -  6.77e-01 5.55e-05h 13
    2562  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.41e+04    -  6.77e-01 5.55e-05h 14
    2563  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.41e+04    -  6.77e-01 5.55e-05h 14
    2564  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.41e+04    -  6.78e-01 5.55e-05h 14
    2565  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.41e+04    -  6.78e-01 5.55e-05h 14
    2566  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.41e+04    -  6.78e-01 5.55e-05h 14
    2567  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.40e+04    -  6.78e-01 5.55e-05h 14
    2568  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.41e+04    -  6.78e-01 5.55e-05h 14
    2569  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.41e+04    -  6.78e-01 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2570  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.41e+04    -  6.78e-01 5.55e-05h 14
    2571  0.0000000e+00 4.16e+05 1.25e+13  -1.0 5.41e+04    -  6.78e-01 4.55e-01w  1
    2572  0.0000000e+00 4.96e+05 7.49e+12  -1.0 9.37e+04    -  4.93e-01 2.59e-01w  1
    2573  0.0000000e+00 2.73e+05 4.89e+12  -1.0 6.63e+03    -  1.00e+00 4.94e-01w  1
    2574  0.0000000e+00 2.57e+04 6.15e+12  -1.0 2.03e+03    -  6.78e-01 5.55e-05h 13
    2575  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.40e+04    -  6.78e-01 5.55e-05h 14
    2576  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.41e+04    -  6.78e-01 5.55e-05h 14
    2577  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.41e+04    -  6.78e-01 5.55e-05h 14
    2578  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.41e+04    -  6.78e-01 5.55e-05h 14
    2579  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.41e+04    -  6.78e-01 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2580  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.39e+04    -  6.78e-01 5.55e-05h 14
    2581  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.41e+04    -  6.79e-01 5.55e-05h 14
    2582  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.41e+04    -  6.79e-01 5.55e-05h 14
    2583  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.41e+04    -  6.79e-01 5.55e-05h 14
    2584  0.0000000e+00 4.16e+05 1.25e+13  -1.0 5.41e+04    -  6.79e-01 4.55e-01w  1
    2585  0.0000000e+00 4.95e+05 7.49e+12  -1.0 9.36e+04    -  4.94e-01 2.59e-01w  1
    2586  0.0000000e+00 2.73e+05 4.89e+12  -1.0 6.62e+03    -  1.00e+00 4.94e-01w  1
    2587  0.0000000e+00 2.57e+04 6.15e+12  -1.0 2.02e+03    -  6.79e-01 5.55e-05h 13
    2588  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.41e+04    -  6.79e-01 5.55e-05h 14
    2589  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.41e+04    -  6.79e-01 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2590  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.41e+04    -  6.79e-01 5.55e-05h 14
    2591  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.41e+04    -  6.79e-01 5.55e-05h 14
    2592  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.41e+04    -  6.79e-01 5.55e-05h 14
    2593  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.41e+04    -  6.79e-01 5.55e-05h 14
    2594  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.41e+04    -  6.79e-01 5.55e-05h 14
    2595  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.40e+04    -  6.79e-01 5.55e-05h 14
    2596  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.41e+04    -  6.79e-01 5.55e-05h 14
    2597  0.0000000e+00 4.16e+05 1.25e+13  -1.0 5.41e+04    -  6.79e-01 4.55e-01w  1
    2598  0.0000000e+00 4.94e+05 7.48e+12  -1.0 9.34e+04    -  4.95e-01 2.59e-01w  1
    2599  0.0000000e+00 2.72e+05 4.90e+12  -1.0 6.61e+03    -  1.00e+00 4.94e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2600  0.0000000e+00 2.57e+04 6.15e+12  -1.0 2.02e+03    -  6.79e-01 5.55e-05h 13
    2601  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.39e+04    -  6.80e-01 5.55e-05h 14
    2602  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.39e+04    -  6.80e-01 5.55e-05h 14
    2603  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.40e+04    -  6.80e-01 5.55e-05h 14
    2604  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.40e+04    -  6.80e-01 5.55e-05h 14
    2605  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.39e+04    -  6.80e-01 5.55e-05h 14
    2606  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.40e+04    -  6.80e-01 5.55e-05h 14
    2607  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.40e+04    -  6.80e-01 5.55e-05h 14
    2608  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.39e+04    -  6.80e-01 5.55e-05h 14
    2609  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.40e+04    -  6.80e-01 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2610  0.0000000e+00 4.15e+05 1.25e+13  -1.0 5.40e+04    -  6.80e-01 4.55e-01w  1
    2611  0.0000000e+00 4.93e+05 7.48e+12  -1.0 9.33e+04    -  4.95e-01 2.60e-01w  1
    2612  0.0000000e+00 2.71e+05 4.90e+12  -1.0 6.60e+03    -  1.00e+00 4.94e-01w  1
    2613  0.0000000e+00 2.57e+04 6.15e+12  -1.0 2.01e+03    -  6.80e-01 5.55e-05h 13
    2614  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.40e+04    -  6.80e-01 5.55e-05h 14
    2615  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.40e+04    -  6.80e-01 5.55e-05h 14
    2616  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.40e+04    -  6.80e-01 5.55e-05h 14
    2617  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.39e+04    -  6.80e-01 5.55e-05h 14
    2618  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.39e+04    -  6.81e-01 5.55e-05h 14
    2619  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.40e+04    -  6.81e-01 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2620  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.40e+04    -  6.81e-01 5.55e-05h 14
    2621  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.40e+04    -  6.81e-01 5.55e-05h 14
    2622  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.40e+04    -  6.81e-01 5.55e-05h 14
    2623  0.0000000e+00 4.15e+05 1.25e+13  -1.0 5.40e+04    -  6.81e-01 4.55e-01w  1
    2624  0.0000000e+00 4.92e+05 7.47e+12  -1.0 9.31e+04    -  4.96e-01 2.60e-01w  1
    2625  0.0000000e+00 2.71e+05 4.91e+12  -1.0 6.59e+03    -  1.00e+00 4.94e-01w  1
    2626  0.0000000e+00 2.57e+04 6.15e+12  -1.0 2.01e+03    -  6.81e-01 5.55e-05h 13
    2627  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.39e+04    -  6.81e-01 5.55e-05h 14
    2628  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.39e+04    -  6.81e-01 5.55e-05h 14
    2629  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.40e+04    -  6.81e-01 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2630  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.39e+04    -  6.81e-01 5.55e-05h 14
    2631  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.40e+04    -  6.81e-01 5.55e-05h 14
    2632  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.40e+04    -  6.81e-01 5.55e-05h 14
    2633  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.39e+04    -  6.81e-01 5.55e-05h 14
    2634  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.39e+04    -  6.81e-01 5.55e-05h 14
    2635  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.40e+04    -  6.82e-01 5.55e-05h 14
    2636  0.0000000e+00 4.14e+05 1.25e+13  -1.0 5.39e+04    -  6.82e-01 4.55e-01w  1
    2637  0.0000000e+00 4.91e+05 7.46e+12  -1.0 9.29e+04    -  4.96e-01 2.60e-01w  1
    2638  0.0000000e+00 2.70e+05 4.91e+12  -1.0 6.57e+03    -  1.00e+00 4.94e-01w  1
    2639  0.0000000e+00 2.57e+04 6.15e+12  -1.0 2.00e+03    -  6.82e-01 5.55e-05h 13
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2640  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.39e+04    -  6.82e-01 5.55e-05h 14
    2641  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.39e+04    -  6.82e-01 5.55e-05h 14
    2642  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.39e+04    -  6.82e-01 5.55e-05h 14
    2643  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.37e+04    -  6.82e-01 5.56e-05h 14
    2644  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.39e+04    -  6.82e-01 5.55e-05h 14
    2645  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.39e+04    -  6.82e-01 5.55e-05h 14
    2646  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.39e+04    -  6.82e-01 5.55e-05h 14
    2647  0.0000000e+00 2.57e+04 6.15e+12  -1.0 5.39e+04    -  6.82e-01 5.55e-05h 14
    2648  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.39e+04    -  6.82e-01 5.55e-05h 14
    2649  0.0000000e+00 4.14e+05 1.25e+13  -1.0 5.37e+04    -  6.82e-01 4.55e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2650  0.0000000e+00 4.91e+05 7.45e+12  -1.0 9.26e+04    -  4.97e-01 2.60e-01w  1
    2651  0.0000000e+00 2.70e+05 4.92e+12  -1.0 6.56e+03    -  1.00e+00 4.94e-01w  1
    2652  0.0000000e+00 2.56e+04 6.15e+12  -1.0 1.99e+03    -  6.82e-01 5.56e-05h 13
    2653  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.39e+04    -  6.83e-01 5.55e-05h 14
    2654  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.39e+04    -  6.83e-01 5.55e-05h 14
    2655  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.39e+04    -  6.83e-01 5.55e-05h 14
    2656  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.39e+04    -  6.83e-01 5.55e-05h 14
    2657  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.39e+04    -  6.83e-01 5.55e-05h 14
    2658  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.39e+04    -  6.83e-01 5.55e-05h 14
    2659  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.38e+04    -  6.83e-01 5.55e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2660  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.39e+04    -  6.83e-01 5.55e-05h 14
    2661  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.38e+04    -  6.83e-01 5.56e-05h 14
    2662  0.0000000e+00 4.13e+05 1.25e+13  -1.0 5.39e+04    -  6.83e-01 4.55e-01w  1
    2663  0.0000000e+00 4.89e+05 7.45e+12  -1.0 9.24e+04    -  4.98e-01 2.61e-01w  1
    2664  0.0000000e+00 2.69e+05 4.92e+12  -1.0 6.55e+03    -  1.00e+00 4.94e-01w  1
    2665  0.0000000e+00 2.56e+04 6.15e+12  -1.0 1.99e+03    -  6.83e-01 5.55e-05h 13
    2666  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.39e+04    -  6.83e-01 5.55e-05h 14
    2667  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.39e+04    -  6.83e-01 5.55e-05h 14
    2668  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.38e+04    -  6.83e-01 5.56e-05h 14
    2669  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.38e+04    -  6.83e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2670  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.36e+04    -  6.84e-01 5.56e-05h 14
    2671  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.39e+04    -  6.84e-01 5.56e-05h 14
    2672  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.38e+04    -  6.84e-01 5.56e-05h 14
    2673  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.39e+04    -  6.84e-01 5.56e-05h 14
    2674  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.37e+04    -  6.84e-01 5.56e-05h 14
    2675  0.0000000e+00 4.12e+05 1.25e+13  -1.0 5.38e+04    -  6.84e-01 4.55e-01w  1
    2676  0.0000000e+00 4.88e+05 7.44e+12  -1.0 9.21e+04    -  5.00e-01 2.62e-01w  1
    2677  0.0000000e+00 2.68e+05 4.93e+12  -1.0 6.53e+03    -  1.00e+00 4.94e-01w  1
    2678  0.0000000e+00 2.56e+04 6.15e+12  -1.0 1.98e+03    -  6.84e-01 5.56e-05h 13
    2679  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.37e+04    -  6.84e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2680  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.38e+04    -  6.84e-01 5.56e-05h 14
    2681  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.38e+04    -  6.84e-01 5.56e-05h 14
    2682  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.38e+04    -  6.84e-01 5.56e-05h 14
    2683  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.39e+04    -  6.84e-01 5.56e-05h 14
    2684  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.38e+04    -  6.84e-01 5.56e-05h 14
    2685  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.37e+04    -  6.85e-01 5.56e-05h 14
    2686  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.38e+04    -  6.84e-01 5.56e-05h 14
    2687  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.38e+04    -  6.85e-01 5.56e-05h 14
    2688  0.0000000e+00 4.12e+05 1.25e+13  -1.0 5.38e+04    -  6.85e-01 4.55e-01w  1
    2689  0.0000000e+00 4.87e+05 7.42e+12  -1.0 9.16e+04    -  5.01e-01 2.62e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2690  0.0000000e+00 2.68e+05 4.93e+12  -1.0 6.51e+03    -  1.00e+00 4.94e-01w  1
    2691  0.0000000e+00 2.56e+04 6.15e+12  -1.0 1.97e+03    -  6.85e-01 5.56e-05h 13
    2692  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.38e+04    -  6.85e-01 5.56e-05h 14
    2693  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.38e+04    -  6.85e-01 5.56e-05h 14
    2694  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.38e+04    -  6.85e-01 5.56e-05h 14
    2695  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.38e+04    -  6.85e-01 5.56e-05h 14
    2696  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.38e+04    -  6.85e-01 5.56e-05h 14
    2697  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.38e+04    -  6.85e-01 5.56e-05h 14
    2698  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.38e+04    -  6.85e-01 5.56e-05h 14
    2699  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.38e+04    -  6.85e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2700  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.37e+04    -  6.85e-01 5.56e-05h 14
    2701  0.0000000e+00 4.11e+05 1.25e+13  -1.0 5.38e+04    -  6.85e-01 4.55e-01w  1
    2702  0.0000000e+00 4.85e+05 7.40e+12  -1.0 9.10e+04    -  5.04e-01 2.64e-01w  1
    2703  0.0000000e+00 2.67e+05 4.93e+12  -1.0 6.49e+03    -  1.00e+00 4.94e-01w  1
    2704  0.0000000e+00 2.56e+04 6.15e+12  -1.0 1.96e+03    -  6.85e-01 5.56e-05h 13
    2705  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.38e+04    -  6.86e-01 5.56e-05h 14
    2706  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.35e+04    -  6.86e-01 5.56e-05h 14
    2707  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.37e+04    -  6.86e-01 5.56e-05h 14
    2708  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.38e+04    -  6.86e-01 5.56e-05h 14
    2709  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.38e+04    -  6.86e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2710  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.37e+04    -  6.86e-01 5.56e-05h 14
    2711  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.37e+04    -  6.86e-01 5.56e-05h 14
    2712  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.37e+04    -  6.86e-01 5.56e-05h 14
    2713  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.37e+04    -  6.86e-01 5.56e-05h 14
    2714  0.0000000e+00 4.11e+05 1.25e+13  -1.0 5.37e+04    -  6.86e-01 4.55e-01w  1
    2715  0.0000000e+00 4.82e+05 7.37e+12  -1.0 9.02e+04    -  5.08e-01 2.66e-01w  1
    2716  0.0000000e+00 2.65e+05 4.94e+12  -1.0 6.46e+03    -  1.00e+00 4.95e-01w  1
    2717  0.0000000e+00 2.56e+04 6.15e+12  -1.0 1.94e+03    -  6.86e-01 5.56e-05h 13
    2718  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.37e+04    -  6.86e-01 5.56e-05h 14
    2719  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.37e+04    -  6.86e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2720  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.36e+04    -  6.87e-01 5.56e-05h 14
    2721  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.37e+04    -  6.87e-01 5.56e-05h 14
    2722  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.35e+04    -  6.87e-01 5.56e-05h 14
    2723  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.37e+04    -  6.87e-01 5.56e-05h 14
    2724  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.37e+04    -  6.87e-01 5.56e-05h 14
    2725  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.37e+04    -  6.87e-01 5.56e-05h 14
    2726  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.37e+04    -  6.87e-01 5.56e-05h 14
    2727  0.0000000e+00 4.10e+05 1.25e+13  -1.0 5.37e+04    -  6.87e-01 4.55e-01w  1
    2728  0.0000000e+00 4.79e+05 7.33e+12  -1.0 8.88e+04    -  5.14e-01 2.68e-01w  1
    2729  0.0000000e+00 2.63e+05 4.93e+12  -1.0 6.42e+03    -  1.00e+00 4.95e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2730  0.0000000e+00 2.56e+04 6.15e+12  -1.0 1.92e+03    -  6.87e-01 5.56e-05h 13
    2731  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.37e+04    -  6.87e-01 5.56e-05h 14
    2732  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.37e+04    -  6.87e-01 5.56e-05h 14
    2733  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.37e+04    -  6.87e-01 5.56e-05h 14
    2734  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.35e+04    -  6.87e-01 5.56e-05h 14
    2735  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.36e+04    -  6.88e-01 5.56e-05h 14
    2736  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.36e+04    -  6.88e-01 5.56e-05h 14
    2737  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.36e+04    -  6.88e-01 5.56e-05h 14
    2738  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.36e+04    -  6.88e-01 5.56e-05h 14
    2739  0.0000000e+00 2.56e+04 6.15e+12  -1.0 5.36e+04    -  6.88e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2740  0.0000000e+00 4.09e+05 1.25e+13  -1.0 5.34e+04    -  6.88e-01 4.55e-01w  1
    2741  0.0000000e+00 4.74e+05 7.26e+12  -1.0 8.68e+04    -  5.23e-01 2.73e-01w  1
    2742  0.0000000e+00 2.60e+05 4.93e+12  -1.0 6.35e+03    -  1.00e+00 4.95e-01w  1
    2743  0.0000000e+00 2.55e+04 6.15e+12  -1.0 1.88e+03    -  6.88e-01 5.56e-05h 13
    2744  0.0000000e+00 2.55e+04 6.15e+12  -1.0 5.35e+04    -  6.88e-01 5.56e-05h 14
    2745  0.0000000e+00 2.55e+04 6.15e+12  -1.0 5.36e+04    -  6.88e-01 5.56e-05h 14
    2746  0.0000000e+00 2.55e+04 6.15e+12  -1.0 5.36e+04    -  6.89e-01 5.56e-05h 14
    2747  0.0000000e+00 2.55e+04 6.15e+12  -1.0 5.34e+04    -  6.89e-01 5.56e-05h 14
    2748  0.0000000e+00 2.55e+04 6.15e+12  -1.0 5.36e+04    -  6.89e-01 5.56e-05h 14
    2749  0.0000000e+00 2.55e+04 6.15e+12  -1.0 5.36e+04    -  6.89e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2750  0.0000000e+00 2.55e+04 6.15e+12  -1.0 5.33e+04    -  6.89e-01 5.56e-05h 14
    2751  0.0000000e+00 2.55e+04 6.15e+12  -1.0 5.34e+04    -  6.89e-01 5.56e-05h 14
    2752  0.0000000e+00 2.55e+04 6.15e+12  -1.0 5.34e+04    -  6.89e-01 5.56e-05h 14
    2753  0.0000000e+00 4.08e+05 1.25e+13  -1.0 5.35e+04    -  6.89e-01 4.55e-01w  1
    2754  0.0000000e+00 4.65e+05 7.16e+12  -1.0 8.39e+04    -  5.37e-01 2.80e-01w  1
    2755  0.0000000e+00 2.55e+05 4.91e+12  -1.0 6.24e+03    -  1.00e+00 4.95e-01w  1
    2756  0.0000000e+00 2.55e+04 6.15e+12  -1.0 1.83e+03    -  6.89e-01 5.56e-05h 13
    2757  0.0000000e+00 2.55e+04 6.15e+12  -1.0 5.35e+04    -  6.89e-01 5.56e-05h 14
    2758  0.0000000e+00 2.55e+04 6.15e+12  -1.0 5.34e+04    -  6.89e-01 5.56e-05h 14
    2759  0.0000000e+00 2.55e+04 6.15e+12  -1.0 5.35e+04    -  6.90e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2760  0.0000000e+00 2.55e+04 6.15e+12  -1.0 5.34e+04    -  6.90e-01 5.56e-05h 14
    2761  0.0000000e+00 2.55e+04 6.15e+12  -1.0 5.35e+04    -  6.90e-01 5.56e-05h 14
    2762  0.0000000e+00 2.55e+04 6.15e+12  -1.0 5.34e+04    -  6.90e-01 5.56e-05h 14
    2763  0.0000000e+00 2.55e+04 6.15e+12  -1.0 5.34e+04    -  6.90e-01 5.56e-05h 14
    2764  0.0000000e+00 2.55e+04 6.15e+12  -1.0 5.34e+04    -  6.90e-01 5.56e-05h 14
    2765  0.0000000e+00 2.55e+04 6.15e+12  -1.0 5.33e+04    -  6.90e-01 5.56e-05h 14
    2766  0.0000000e+00 4.06e+05 1.25e+13  -1.0 5.33e+04    -  6.91e-01 4.55e-01w  1
    2767  0.0000000e+00 4.54e+05 7.01e+12  -1.0 7.95e+04    -  5.59e-01 2.91e-01w  1
    2768  0.0000000e+00 2.48e+05 4.88e+12  -1.0 6.05e+03    -  1.00e+00 4.95e-01w  1
    2769  0.0000000e+00 2.55e+04 6.15e+12  -1.0 1.73e+03    -  6.91e-01 5.56e-05h 13
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2770  0.0000000e+00 2.55e+04 6.15e+12  -1.0 5.34e+04    -  6.91e-01 5.56e-05h 14
    2771  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.33e+04    -  6.91e-01 5.56e-05h 14
    2772  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.33e+04    -  6.91e-01 5.56e-05h 14
    2773  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.33e+04    -  6.91e-01 5.56e-05h 14
    2774  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.32e+04    -  6.91e-01 5.56e-05h 14
    2775  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.33e+04    -  6.91e-01 5.56e-05h 14
    2776  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.33e+04    -  6.92e-01 5.56e-05h 14
    2777  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.32e+04    -  6.92e-01 5.56e-05h 14
    2778  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.31e+04    -  6.92e-01 5.56e-05h 14
    2779  0.0000000e+00 4.04e+05 1.25e+13  -1.0 5.32e+04    -  6.92e-01 4.55e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2780  0.0000000e+00 4.35e+05 6.78e+12  -1.0 7.33e+04    -  5.94e-01 3.08e-01w  1
    2781  0.0000000e+00 2.36e+05 4.81e+12  -1.0 5.67e+03    -  1.00e+00 4.95e-01w  1
    2782  0.0000000e+00 2.55e+04 6.16e+12  -1.0 1.55e+03    -  6.92e-01 5.56e-05h 13
    2783  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.32e+04    -  6.92e-01 5.56e-05h 14
    2784  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.32e+04    -  6.92e-01 5.56e-05h 14
    2785  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.31e+04    -  6.93e-01 5.56e-05h 14
    2786  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.30e+04    -  6.93e-01 5.56e-05h 14
    2787  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.31e+04    -  6.93e-01 5.56e-05h 14
    2788  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.30e+04    -  6.93e-01 5.56e-05h 14
    2789  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.28e+04    -  6.93e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2790  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.30e+04    -  6.94e-01 5.56e-05h 14
    2791  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.28e+04    -  6.94e-01 5.56e-05h 14
    2792  0.0000000e+00 4.02e+05 1.25e+13  -1.0 5.28e+04    -  6.94e-01 4.55e-01w  1
    2793  0.0000000e+00 4.08e+05 6.45e+12  -1.0 6.47e+04    -  6.51e-01 3.35e-01w  1
    2794  0.0000000e+00 2.18e+05 4.68e+12  -1.0 4.92e+03    -  1.00e+00 4.95e-01w  1
    2795  0.0000000e+00 2.55e+04 6.16e+12  -1.0 1.22e+03    -  6.94e-01 5.56e-05h 13
    2796  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.29e+04    -  6.94e-01 5.56e-05h 14
    2797  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.29e+04    -  6.94e-01 5.56e-05h 14
    2798  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.26e+04    -  6.95e-01 5.56e-05h 14
    2799  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.26e+04    -  6.95e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2800  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.28e+04    -  6.95e-01 5.56e-05h 14
    2801  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.26e+04    -  6.95e-01 5.56e-05h 14
    2802  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.25e+04    -  6.96e-01 5.56e-05h 14
    2803  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.26e+04    -  6.96e-01 5.56e-05h 14
    2804  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.32e+04    -  6.93e-01 5.56e-05h 14
    2805  0.0000000e+00 4.08e+05 1.26e+13  -1.0 5.34e+04    -  6.92e-01 4.55e-01w  1
    2806  0.0000000e+00 3.80e+05 6.16e+12  -1.0 5.58e+04    -  7.19e-01 3.69e-01w  1
    2807  0.0000000e+00 1.98e+05 4.43e+12  -1.0 3.73e+03    -  1.00e+00 4.95e-01w  1
    2808  0.0000000e+00 2.55e+04 6.16e+12  -1.0 1.81e+03    -  6.92e-01 5.56e-05h 13
    2809  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.33e+04    -  6.91e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2810  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.34e+04    -  6.91e-01 5.56e-05h 14
    2811  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.36e+04    -  6.91e-01 5.56e-05h 14
    2812  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.34e+04    -  6.91e-01 5.56e-05h 14
    2813  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.36e+04    -  6.92e-01 5.56e-05h 14
    2814  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.36e+04    -  6.91e-01 5.56e-05h 14
    2815  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.36e+04    -  6.92e-01 5.56e-05h 14
    2816  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.34e+04    -  6.92e-01 5.56e-05h 14
    2817  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.35e+04    -  6.92e-01 5.56e-05h 14
    2818  0.0000000e+00 4.07e+05 1.26e+13  -1.0 5.33e+04    -  6.92e-01 4.55e-01w  1
    2819  0.0000000e+00 3.79e+05 6.21e+12  -1.0 5.56e+04    -  7.15e-01 3.69e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2820  0.0000000e+00 1.98e+05 4.42e+12  -1.0 3.69e+03    -  1.00e+00 4.95e-01w  1
    2821  0.0000000e+00 2.55e+04 6.16e+12  -1.0 1.82e+03    -  6.92e-01 5.56e-05h 13
    2822  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.35e+04    -  6.92e-01 5.56e-05h 14
    2823  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.36e+04    -  6.92e-01 5.56e-05h 14
    2824  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.35e+04    -  6.92e-01 5.56e-05h 14
    2825  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.34e+04    -  6.92e-01 5.56e-05h 14
    2826  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.35e+04    -  6.92e-01 5.56e-05h 14
    2827  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.35e+04    -  6.92e-01 5.56e-05h 14
    2828  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.34e+04    -  6.92e-01 5.56e-05h 14
    2829  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.35e+04    -  6.92e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2830  0.0000000e+00 2.55e+04 6.16e+12  -1.0 5.35e+04    -  6.93e-01 5.56e-05h 14
    2831  0.0000000e+00 4.07e+05 1.26e+13  -1.0 5.35e+04    -  6.92e-01 4.55e-01w  1
    2832  0.0000000e+00 3.78e+05 6.22e+12  -1.0 5.56e+04    -  7.16e-01 3.70e-01w  1
    2833  0.0000000e+00 1.97e+05 4.42e+12  -1.0 3.67e+03    -  1.00e+00 4.95e-01w  1
    2834  0.0000000e+00 2.54e+04 6.16e+12  -1.0 1.86e+03    -  6.92e-01 5.56e-05h 13
    2835  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.35e+04    -  6.93e-01 5.56e-05h 14
    2836  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.35e+04    -  6.93e-01 5.56e-05h 14
    2837  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.35e+04    -  6.93e-01 5.56e-05h 14
    2838  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.34e+04    -  6.93e-01 5.56e-05h 14
    2839  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.34e+04    -  6.93e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2840  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.34e+04    -  6.93e-01 5.56e-05h 14
    2841  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.35e+04    -  6.93e-01 5.56e-05h 14
    2842  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.33e+04    -  6.93e-01 5.56e-05h 14
    2843  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.34e+04    -  6.93e-01 5.56e-05h 14
    2844  0.0000000e+00 4.06e+05 1.26e+13  -1.0 5.35e+04    -  6.93e-01 4.55e-01w  1
    2845  0.0000000e+00 3.78e+05 6.22e+12  -1.0 5.55e+04    -  7.16e-01 3.70e-01w  1
    2846  0.0000000e+00 1.97e+05 4.42e+12  -1.0 3.66e+03    -  1.00e+00 4.95e-01w  1
    2847  0.0000000e+00 2.54e+04 6.16e+12  -1.0 1.87e+03    -  6.93e-01 5.56e-05h 13
    2848  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.35e+04    -  6.93e-01 5.56e-05h 14
    2849  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.35e+04    -  6.93e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2850  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.35e+04    -  6.94e-01 5.56e-05h 14
    2851  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.35e+04    -  6.93e-01 5.56e-05h 14
    2852  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.35e+04    -  6.94e-01 5.56e-05h 14
    2853  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.35e+04    -  6.94e-01 5.56e-05h 14
    2854  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.35e+04    -  6.94e-01 5.56e-05h 14
    2855  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.34e+04    -  6.94e-01 5.56e-05h 14
    2856  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.31e+04    -  6.94e-01 5.56e-05h 14
    2857  0.0000000e+00 4.06e+05 1.26e+13  -1.0 5.33e+04    -  6.94e-01 4.55e-01w  1
    2858  0.0000000e+00 3.77e+05 6.21e+12  -1.0 5.54e+04    -  7.17e-01 3.70e-01w  1
    2859  0.0000000e+00 1.96e+05 4.43e+12  -1.0 3.64e+03    -  1.00e+00 4.95e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2860  0.0000000e+00 2.54e+04 6.16e+12  -1.0 1.86e+03    -  6.94e-01 5.56e-05h 13
    2861  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.34e+04    -  6.94e-01 5.56e-05h 14
    2862  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.35e+04    -  6.94e-01 5.56e-05h 14
    2863  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.34e+04    -  6.94e-01 5.56e-05h 14
    2864  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.34e+04    -  6.94e-01 5.56e-05h 14
    2865  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.35e+04    -  6.94e-01 5.56e-05h 14
    2866  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.35e+04    -  6.94e-01 5.56e-05h 14
    2867  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.33e+04    -  6.95e-01 5.56e-05h 14
    2868  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.34e+04    -  6.95e-01 5.56e-05h 14
    2869  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.35e+04    -  6.95e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2870  0.0000000e+00 4.05e+05 1.26e+13  -1.0 5.35e+04    -  6.95e-01 4.55e-01w  1
    2871  0.0000000e+00 3.76e+05 6.21e+12  -1.0 5.54e+04    -  7.17e-01 3.70e-01w  1
    2872  0.0000000e+00 1.96e+05 4.43e+12  -1.0 3.63e+03    -  1.00e+00 4.95e-01w  1
    2873  0.0000000e+00 2.54e+04 6.16e+12  -1.0 1.91e+03    -  6.95e-01 5.56e-05h 13
    2874  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.35e+04    -  6.95e-01 5.56e-05h 14
    2875  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.33e+04    -  6.95e-01 5.56e-05h 14
    2876  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.32e+04    -  6.95e-01 5.56e-05h 14
    2877  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.34e+04    -  6.95e-01 5.56e-05h 14
    2878  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.34e+04    -  6.95e-01 5.56e-05h 14
    2879  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.34e+04    -  6.95e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2880  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.34e+04    -  6.95e-01 5.56e-05h 14
    2881  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.34e+04    -  6.95e-01 5.56e-05h 14
    2882  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.33e+04    -  6.95e-01 5.56e-05h 14
    2883  0.0000000e+00 4.05e+05 1.26e+13  -1.0 5.34e+04    -  6.95e-01 4.55e-01w  1
    2884  0.0000000e+00 3.76e+05 6.21e+12  -1.0 5.53e+04    -  7.18e-01 3.71e-01w  1
    2885  0.0000000e+00 1.96e+05 4.43e+12  -1.0 3.61e+03    -  1.00e+00 4.95e-01w  1
    2886  0.0000000e+00 2.54e+04 6.16e+12  -1.0 1.92e+03    -  6.95e-01 5.56e-05h 13
    2887  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.34e+04    -  6.96e-01 5.56e-05h 14
    2888  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.31e+04    -  6.96e-01 5.56e-05h 14
    2889  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.32e+04    -  6.96e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2890  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.34e+04    -  6.96e-01 5.56e-05h 14
    2891  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.34e+04    -  6.96e-01 5.56e-05h 14
    2892  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.33e+04    -  6.96e-01 5.56e-05h 14
    2893  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.34e+04    -  6.96e-01 5.56e-05h 14
    2894  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.34e+04    -  6.96e-01 5.56e-05h 14
    2895  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.34e+04    -  6.96e-01 5.56e-05h 14
    2896  0.0000000e+00 4.04e+05 1.26e+13  -1.0 5.34e+04    -  6.96e-01 4.55e-01w  1
    2897  0.0000000e+00 3.75e+05 6.21e+12  -1.0 5.52e+04    -  7.19e-01 3.71e-01w  1
    2898  0.0000000e+00 1.95e+05 4.44e+12  -1.0 3.60e+03    -  1.00e+00 4.95e-01w  1
    2899  0.0000000e+00 2.54e+04 6.16e+12  -1.0 1.93e+03    -  6.96e-01 5.56e-05h 13
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2900  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.33e+04    -  6.96e-01 5.56e-05h 14
    2901  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.33e+04    -  6.96e-01 5.56e-05h 14
    2902  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.34e+04    -  6.96e-01 5.56e-05h 14
    2903  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.34e+04    -  6.96e-01 5.56e-05h 14
    2904  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.33e+04    -  6.97e-01 5.56e-05h 14
    2905  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.34e+04    -  6.97e-01 5.56e-05h 14
    2906  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.34e+04    -  6.97e-01 5.56e-05h 14
    2907  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.34e+04    -  6.97e-01 5.56e-05h 14
    2908  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.33e+04    -  6.97e-01 5.56e-05h 14
    2909  0.0000000e+00 4.04e+05 1.26e+13  -1.0 5.33e+04    -  6.97e-01 4.55e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2910  0.0000000e+00 3.74e+05 6.21e+12  -1.0 5.51e+04    -  7.19e-01 3.71e-01w  1
    2911  0.0000000e+00 1.95e+05 4.44e+12  -1.0 3.59e+03    -  1.00e+00 4.95e-01w  1
    2912  0.0000000e+00 2.54e+04 6.16e+12  -1.0 1.93e+03    -  6.97e-01 5.56e-05h 13
    2913  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.33e+04    -  6.97e-01 5.56e-05h 14
    2914  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.33e+04    -  6.97e-01 5.56e-05h 14
    2915  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.30e+04    -  6.97e-01 5.56e-05h 14
    2916  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.33e+04    -  6.97e-01 5.56e-05h 14
    2917  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.33e+04    -  6.97e-01 5.56e-05h 14
    2918  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.32e+04    -  6.97e-01 5.56e-05h 14
    2919  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.33e+04    -  6.97e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2920  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.30e+04    -  6.97e-01 5.56e-05h 14
    2921  0.0000000e+00 2.54e+04 6.16e+12  -1.0 5.33e+04    -  6.98e-01 5.56e-05h 14
    2922  0.0000000e+00 4.03e+05 1.26e+13  -1.0 5.33e+04    -  6.98e-01 4.56e-01w  1
    2923  0.0000000e+00 3.74e+05 6.20e+12  -1.0 5.51e+04    -  7.20e-01 3.71e-01w  1
    2924  0.0000000e+00 1.95e+05 4.44e+12  -1.0 3.57e+03    -  1.00e+00 4.95e-01w  1
    2925  0.0000000e+00 2.54e+04 6.16e+12  -1.0 1.95e+03    -  6.98e-01 5.56e-05h 13
    2926  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  6.98e-01 5.56e-05h 14
    2927  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.33e+04    -  6.98e-01 5.56e-05h 14
    2928  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.33e+04    -  6.98e-01 5.56e-05h 14
    2929  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.33e+04    -  6.98e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2930  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  6.98e-01 5.56e-05h 14
    2931  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.33e+04    -  6.98e-01 5.56e-05h 14
    2932  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.33e+04    -  6.98e-01 5.56e-05h 14
    2933  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.33e+04    -  6.98e-01 5.56e-05h 14
    2934  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  6.98e-01 5.56e-05h 14
    2935  0.0000000e+00 4.03e+05 1.26e+13  -1.0 5.33e+04    -  6.98e-01 4.56e-01w  1
    2936  0.0000000e+00 3.73e+05 6.20e+12  -1.0 5.50e+04    -  7.20e-01 3.72e-01w  1
    2937  0.0000000e+00 1.94e+05 4.45e+12  -1.0 3.56e+03    -  1.00e+00 4.95e-01w  1
    2938  0.0000000e+00 2.53e+04 6.16e+12  -1.0 1.97e+03    -  6.98e-01 5.56e-05h 13
    2939  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.33e+04    -  6.98e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2940  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.33e+04    -  6.98e-01 5.56e-05h 14
    2941  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.33e+04    -  6.99e-01 5.56e-05h 14
    2942  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.33e+04    -  6.99e-01 5.56e-05h 14
    2943  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.33e+04    -  6.99e-01 5.56e-05h 14
    2944  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  6.99e-01 5.56e-05h 14
    2945  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  6.99e-01 5.56e-05h 14
    2946  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.33e+04    -  6.99e-01 5.56e-05h 14
    2947  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  6.99e-01 5.56e-05h 14
    2948  0.0000000e+00 4.02e+05 1.26e+13  -1.0 5.32e+04    -  6.99e-01 4.56e-01w  1
    2949  0.0000000e+00 3.72e+05 6.20e+12  -1.0 5.49e+04    -  7.21e-01 3.72e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2950  0.0000000e+00 1.94e+05 4.45e+12  -1.0 3.54e+03    -  1.00e+00 4.95e-01w  1
    2951  0.0000000e+00 2.53e+04 6.16e+12  -1.0 1.97e+03    -  6.99e-01 5.56e-05h 13
    2952  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.33e+04    -  6.99e-01 5.56e-05h 14
    2953  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.33e+04    -  6.99e-01 5.56e-05h 14
    2954  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.31e+04    -  6.99e-01 5.56e-05h 14
    2955  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.33e+04    -  6.99e-01 5.56e-05h 14
    2956  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.33e+04    -  6.99e-01 5.56e-05h 14
    2957  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.30e+04    -  6.99e-01 5.56e-05h 14
    2958  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  7.00e-01 5.56e-05h 14
    2959  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  7.00e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2960  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.31e+04    -  7.00e-01 5.56e-05h 14
    2961  0.0000000e+00 4.02e+05 1.26e+13  -1.0 5.29e+04    -  7.00e-01 4.56e-01w  1
    2962  0.0000000e+00 3.72e+05 6.20e+12  -1.0 5.47e+04    -  7.21e-01 3.72e-01w  1
    2963  0.0000000e+00 1.94e+05 4.46e+12  -1.0 3.53e+03    -  1.00e+00 4.95e-01w  1
    2964  0.0000000e+00 2.53e+04 6.16e+12  -1.0 1.95e+03    -  7.00e-01 5.56e-05h 13
    2965  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.30e+04    -  7.00e-01 5.56e-05h 14
    2966  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.31e+04    -  7.00e-01 5.56e-05h 14
    2967  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  7.00e-01 5.56e-05h 14
    2968  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  7.00e-01 5.56e-05h 14
    2969  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  7.00e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2970  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  7.00e-01 5.56e-05h 14
    2971  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  7.00e-01 5.56e-05h 14
    2972  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  7.00e-01 5.56e-05h 14
    2973  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  7.00e-01 5.56e-05h 14
    2974  0.0000000e+00 4.02e+05 1.26e+13  -1.0 5.31e+04    -  7.00e-01 4.56e-01w  1
    2975  0.0000000e+00 3.71e+05 6.20e+12  -1.0 5.48e+04    -  7.22e-01 3.72e-01w  1
    2976  0.0000000e+00 1.93e+05 4.46e+12  -1.0 3.58e+03    -  1.00e+00 4.95e-01w  1
    2977  0.0000000e+00 2.53e+04 6.16e+12  -1.0 2.00e+03    -  7.00e-01 5.56e-05h 13
    2978  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  7.01e-01 5.56e-05h 14
    2979  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  7.01e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2980  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  7.01e-01 5.56e-05h 14
    2981  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.31e+04    -  7.01e-01 5.56e-05h 14
    2982  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  7.01e-01 5.56e-05h 14
    2983  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  7.01e-01 5.56e-05h 14
    2984  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.31e+04    -  7.01e-01 5.56e-05h 14
    2985  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  7.01e-01 5.56e-05h 14
    2986  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  7.01e-01 5.56e-05h 14
    2987  0.0000000e+00 4.01e+05 1.26e+13  -1.0 5.32e+04    -  7.01e-01 4.56e-01w  1
    2988  0.0000000e+00 3.70e+05 6.20e+12  -1.0 5.47e+04    -  7.23e-01 3.73e-01w  1
    2989  0.0000000e+00 1.93e+05 4.46e+12  -1.0 3.65e+03    -  1.00e+00 4.95e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2990  0.0000000e+00 2.53e+04 6.16e+12  -1.0 2.02e+03    -  7.01e-01 5.56e-05h 13
    2991  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  7.01e-01 5.56e-05h 14
    2992  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.30e+04    -  7.01e-01 5.56e-05h 14
    2993  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  7.01e-01 5.56e-05h 14
    2994  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  7.01e-01 5.56e-05h 14
    2995  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.31e+04    -  7.02e-01 5.56e-05h 14
    2996  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  7.02e-01 5.56e-05h 14
    2997  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.30e+04    -  7.02e-01 5.56e-05h 14
    2998  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  7.02e-01 5.56e-05h 14
    2999  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.31e+04    -  7.02e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3000  0.0000000e+00 4.01e+05 1.26e+13  -1.0 5.32e+04    -  7.02e-01 4.56e-01w  1
    3001  0.0000000e+00 3.70e+05 6.20e+12  -1.0 5.46e+04    -  7.23e-01 3.73e-01w  1
    3002  0.0000000e+00 1.92e+05 4.47e+12  -1.0 3.67e+03    -  1.00e+00 4.95e-01w  1
    3003  0.0000000e+00 2.53e+04 6.16e+12  -1.0 2.03e+03    -  7.02e-01 5.56e-05h 13
    3004  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  7.02e-01 5.56e-05h 14
    3005  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  7.02e-01 5.56e-05h 14
    3006  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  7.02e-01 5.56e-05h 14
    3007  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.32e+04    -  7.02e-01 5.56e-05h 14
    3008  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.29e+04    -  7.02e-01 5.56e-05h 14
    3009  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.31e+04    -  7.02e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3010  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.31e+04    -  7.02e-01 5.56e-05h 14
    3011  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.31e+04    -  7.02e-01 5.56e-05h 14
    3012  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.31e+04    -  7.03e-01 5.56e-05h 14
    3013  0.0000000e+00 4.00e+05 1.26e+13  -1.0 5.31e+04    -  7.03e-01 4.56e-01w  1
    3014  0.0000000e+00 3.69e+05 6.19e+12  -1.0 5.46e+04    -  7.24e-01 3.73e-01w  1
    3015  0.0000000e+00 1.92e+05 4.47e+12  -1.0 3.69e+03    -  1.00e+00 4.95e-01w  1
    3016  0.0000000e+00 2.53e+04 6.16e+12  -1.0 2.05e+03    -  7.03e-01 5.56e-05h 13
    3017  0.0000000e+00 2.53e+04 6.16e+12  -1.0 5.31e+04    -  7.03e-01 5.56e-05h 14
    3018  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.31e+04    -  7.03e-01 5.56e-05h 14
    3019  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.31e+04    -  7.03e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3020  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.31e+04    -  7.03e-01 5.56e-05h 14
    3021  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.31e+04    -  7.03e-01 5.56e-05h 14
    3022  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.31e+04    -  7.03e-01 5.56e-05h 14
    3023  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.30e+04    -  7.03e-01 5.56e-05h 14
    3024  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.31e+04    -  7.03e-01 5.56e-05h 14
    3025  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.30e+04    -  7.03e-01 5.56e-05h 14
    3026  0.0000000e+00 4.00e+05 1.26e+13  -1.0 5.31e+04    -  7.03e-01 4.56e-01w  1
    3027  0.0000000e+00 3.69e+05 6.19e+12  -1.0 5.45e+04    -  7.24e-01 3.73e-01w  1
    3028  0.0000000e+00 1.92e+05 4.47e+12  -1.0 3.72e+03    -  1.00e+00 4.95e-01w  1
    3029  0.0000000e+00 2.52e+04 6.16e+12  -1.0 2.06e+03    -  7.03e-01 5.56e-05h 13
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3030  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.31e+04    -  7.03e-01 5.56e-05h 14
    3031  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.31e+04    -  7.03e-01 5.56e-05h 14
    3032  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.31e+04    -  7.04e-01 5.56e-05h 14
    3033  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.31e+04    -  7.04e-01 5.56e-05h 14
    3034  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.31e+04    -  7.04e-01 5.56e-05h 14
    3035  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.28e+04    -  7.04e-01 5.56e-05h 14
    3036  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.29e+04    -  7.04e-01 5.56e-05h 14
    3037  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.31e+04    -  7.04e-01 5.56e-05h 14
    3038  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.31e+04    -  7.04e-01 5.56e-05h 14
    3039  0.0000000e+00 3.99e+05 1.26e+13  -1.0 5.30e+04    -  7.04e-01 4.56e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3040  0.0000000e+00 3.68e+05 6.19e+12  -1.0 5.44e+04    -  7.25e-01 3.74e-01w  1
    3041  0.0000000e+00 1.91e+05 4.48e+12  -1.0 3.71e+03    -  1.00e+00 4.95e-01w  1
    3042  0.0000000e+00 2.52e+04 6.16e+12  -1.0 2.06e+03    -  7.04e-01 5.56e-05h 13
    3043  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.29e+04    -  7.04e-01 5.56e-05h 14
    3044  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.31e+04    -  7.04e-01 5.56e-05h 14
    3045  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.30e+04    -  7.04e-01 5.56e-05h 14
    3046  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.31e+04    -  7.04e-01 5.56e-05h 14
    3047  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.31e+04    -  7.04e-01 5.56e-05h 14
    3048  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.31e+04    -  7.04e-01 5.56e-05h 14
    3049  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.31e+04    -  7.05e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3050  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.29e+04    -  7.05e-01 5.56e-05h 14
    3051  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.30e+04    -  7.05e-01 5.56e-05h 14
    3052  0.0000000e+00 3.99e+05 1.26e+13  -1.0 5.31e+04    -  7.05e-01 4.56e-01w  1
    3053  0.0000000e+00 3.67e+05 6.19e+12  -1.0 5.44e+04    -  7.26e-01 3.74e-01w  1
    3054  0.0000000e+00 1.91e+05 4.48e+12  -1.0 3.77e+03    -  1.00e+00 4.95e-01w  1
    3055  0.0000000e+00 2.52e+04 6.16e+12  -1.0 2.09e+03    -  7.05e-01 5.56e-05h 13
    3056  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.30e+04    -  7.05e-01 5.56e-05h 14
    3057  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.30e+04    -  7.05e-01 5.56e-05h 14
    3058  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.31e+04    -  7.05e-01 5.56e-05h 14
    3059  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.30e+04    -  7.05e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3060  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.29e+04    -  7.05e-01 5.56e-05h 14
    3061  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.30e+04    -  7.05e-01 5.56e-05h 14
    3062  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.30e+04    -  7.05e-01 5.56e-05h 14
    3063  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.30e+04    -  7.05e-01 5.56e-05h 14
    3064  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.27e+04    -  7.05e-01 5.56e-05h 14
    3065  0.0000000e+00 3.98e+05 1.26e+13  -1.0 5.29e+04    -  7.06e-01 4.56e-01w  1
    3066  0.0000000e+00 3.67e+05 6.19e+12  -1.0 5.42e+04    -  7.26e-01 3.74e-01w  1
    3067  0.0000000e+00 1.90e+05 4.48e+12  -1.0 3.74e+03    -  1.00e+00 4.95e-01w  1
    3068  0.0000000e+00 2.52e+04 6.16e+12  -1.0 2.08e+03    -  7.06e-01 5.56e-05h 13
    3069  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.30e+04    -  7.06e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3070  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.28e+04    -  7.06e-01 5.56e-05h 14
    3071  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.30e+04    -  7.06e-01 5.56e-05h 14
    3072  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.29e+04    -  7.06e-01 5.56e-05h 14
    3073  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.29e+04    -  7.06e-01 5.56e-05h 14
    3074  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.30e+04    -  7.06e-01 5.56e-05h 14
    3075  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.29e+04    -  7.06e-01 5.56e-05h 14
    3076  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.30e+04    -  7.06e-01 5.56e-05h 14
    3077  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.29e+04    -  7.06e-01 5.56e-05h 14
    3078  0.0000000e+00 3.98e+05 1.26e+13  -1.0 5.30e+04    -  7.06e-01 4.56e-01w  1
    3079  0.0000000e+00 3.66e+05 6.19e+12  -1.0 5.42e+04    -  7.27e-01 3.74e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3080  0.0000000e+00 1.90e+05 4.49e+12  -1.0 3.81e+03    -  1.00e+00 4.95e-01w  1
    3081  0.0000000e+00 2.52e+04 6.16e+12  -1.0 2.11e+03    -  7.06e-01 5.56e-05h 13
    3082  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.29e+04    -  7.06e-01 5.56e-05h 14
    3083  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.29e+04    -  7.06e-01 5.56e-05h 14
    3084  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.29e+04    -  7.07e-01 5.56e-05h 14
    3085  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.30e+04    -  7.06e-01 5.56e-05h 14
    3086  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.28e+04    -  7.07e-01 5.56e-05h 14
    3087  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.30e+04    -  7.07e-01 5.56e-05h 14
    3088  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.29e+04    -  7.07e-01 5.56e-05h 14
    3089  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.30e+04    -  7.07e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3090  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.30e+04    -  7.07e-01 5.56e-05h 14
    3091  0.0000000e+00 3.97e+05 1.26e+13  -1.0 5.29e+04    -  7.07e-01 4.56e-01w  1
    3092  0.0000000e+00 3.65e+05 6.19e+12  -1.0 5.41e+04    -  7.27e-01 3.75e-01w  1
    3093  0.0000000e+00 1.90e+05 4.49e+12  -1.0 3.83e+03    -  1.00e+00 4.95e-01w  1
    3094  0.0000000e+00 2.52e+04 6.16e+12  -1.0 2.12e+03    -  7.07e-01 5.56e-05h 13
    3095  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.26e+04    -  7.07e-01 5.56e-05h 14
    3096  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.29e+04    -  7.07e-01 5.56e-05h 14
    3097  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.29e+04    -  7.07e-01 5.56e-05h 14
    3098  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.28e+04    -  7.07e-01 5.56e-05h 14
    3099  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.29e+04    -  7.07e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3100  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.29e+04    -  7.07e-01 5.56e-05h 14
    3101  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.30e+04    -  7.07e-01 5.56e-05h 14
    3102  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.29e+04    -  7.07e-01 5.56e-05h 14
    3103  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.29e+04    -  7.08e-01 5.56e-05h 14
    3104  0.0000000e+00 3.97e+05 1.26e+13  -1.0 5.29e+04    -  7.08e-01 4.56e-01w  1
    3105  0.0000000e+00 3.65e+05 6.18e+12  -1.0 5.40e+04    -  7.28e-01 3.75e-01w  1
    3106  0.0000000e+00 1.89e+05 4.50e+12  -1.0 3.86e+03    -  1.00e+00 4.95e-01w  1
    3107  0.0000000e+00 2.52e+04 6.16e+12  -1.0 2.13e+03    -  7.08e-01 5.56e-05h 13
    3108  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.29e+04    -  7.08e-01 5.56e-05h 14
    3109  0.0000000e+00 2.52e+04 6.16e+12  -1.0 5.29e+04    -  7.08e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3110  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.29e+04    -  7.08e-01 5.56e-05h 14
    3111  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.29e+04    -  7.08e-01 5.56e-05h 14
    3112  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.24e+04    -  7.08e-01 5.57e-05h 14
    3113  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.29e+04    -  7.08e-01 5.56e-05h 14
    3114  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.29e+04    -  7.08e-01 5.56e-05h 14
    3115  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.29e+04    -  7.08e-01 5.56e-05h 14
    3116  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.29e+04    -  7.08e-01 5.56e-05h 14
    3117  0.0000000e+00 3.96e+05 1.26e+13  -1.0 5.29e+04    -  7.08e-01 4.56e-01w  1
    3118  0.0000000e+00 3.64e+05 6.18e+12  -1.0 5.40e+04    -  7.29e-01 3.75e-01w  1
    3119  0.0000000e+00 1.89e+05 4.50e+12  -1.0 3.89e+03    -  1.00e+00 4.95e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3120  0.0000000e+00 2.51e+04 6.16e+12  -1.0 2.15e+03    -  7.08e-01 5.56e-05h 13
    3121  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.27e+04    -  7.08e-01 5.56e-05h 14
    3122  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.09e-01 5.56e-05h 14
    3123  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.29e+04    -  7.09e-01 5.56e-05h 14
    3124  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.09e-01 5.56e-05h 14
    3125  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.29e+04    -  7.09e-01 5.56e-05h 14
    3126  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.09e-01 5.56e-05h 14
    3127  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.29e+04    -  7.09e-01 5.56e-05h 14
    3128  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.29e+04    -  7.09e-01 5.56e-05h 14
    3129  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.09e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3130  0.0000000e+00 3.96e+05 1.26e+13  -1.0 5.29e+04    -  7.09e-01 4.56e-01w  1
    3131  0.0000000e+00 3.63e+05 6.18e+12  -1.0 5.39e+04    -  7.29e-01 3.75e-01w  1
    3132  0.0000000e+00 1.89e+05 4.50e+12  -1.0 3.91e+03    -  1.00e+00 4.95e-01w  1
    3133  0.0000000e+00 2.51e+04 6.16e+12  -1.0 2.16e+03    -  7.09e-01 5.56e-05h 13
    3134  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.29e+04    -  7.09e-01 5.56e-05h 14
    3135  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.09e-01 5.56e-05h 14
    3136  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.09e-01 5.56e-05h 14
    3137  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.29e+04    -  7.09e-01 5.56e-05h 14
    3138  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.29e+04    -  7.09e-01 5.56e-05h 14
    3139  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.09e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3140  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.10e-01 5.56e-05h 14
    3141  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.29e+04    -  7.10e-01 5.56e-05h 14
    3142  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.10e-01 5.56e-05h 14
    3143  0.0000000e+00 3.95e+05 1.26e+13  -1.0 5.28e+04    -  7.10e-01 4.56e-01w  1
    3144  0.0000000e+00 3.63e+05 6.18e+12  -1.0 5.38e+04    -  7.30e-01 3.76e-01w  1
    3145  0.0000000e+00 1.88e+05 4.51e+12  -1.0 3.92e+03    -  1.00e+00 4.95e-01w  1
    3146  0.0000000e+00 2.51e+04 6.16e+12  -1.0 2.17e+03    -  7.10e-01 5.56e-05h 13
    3147  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.29e+04    -  7.10e-01 5.56e-05h 14
    3148  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.29e+04    -  7.10e-01 5.56e-05h 14
    3149  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.29e+04    -  7.10e-01 5.56e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3150  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.29e+04    -  7.10e-01 5.56e-05h 14
    3151  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.10e-01 5.56e-05h 14
    3152  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.10e-01 5.57e-05h 14
    3153  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.10e-01 5.56e-05h 14
    3154  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.10e-01 5.57e-05h 14
    3155  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.26e+04    -  7.11e-01 5.57e-05h 14
    3156  0.0000000e+00 3.95e+05 1.26e+13  -1.0 5.28e+04    -  7.11e-01 4.56e-01w  1
    3157  0.0000000e+00 3.62e+05 6.18e+12  -1.0 5.38e+04    -  7.30e-01 3.76e-01w  1
    3158  0.0000000e+00 1.88e+05 4.51e+12  -1.0 3.96e+03    -  1.00e+00 4.95e-01w  1
    3159  0.0000000e+00 2.51e+04 6.16e+12  -1.0 2.19e+03    -  7.11e-01 5.57e-05h 13
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3160  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.11e-01 5.57e-05h 14
    3161  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.11e-01 5.57e-05h 14
    3162  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.11e-01 5.57e-05h 14
    3163  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.27e+04    -  7.11e-01 5.57e-05h 14
    3164  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.11e-01 5.57e-05h 14
    3165  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.27e+04    -  7.11e-01 5.57e-05h 14
    3166  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.26e+04    -  7.11e-01 5.57e-05h 14
    3167  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.11e-01 5.57e-05h 14
    3168  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.27e+04    -  7.11e-01 5.57e-05h 14
    3169  0.0000000e+00 3.94e+05 1.26e+13  -1.0 5.28e+04    -  7.11e-01 4.56e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3170  0.0000000e+00 3.62e+05 6.18e+12  -1.0 5.37e+04    -  7.31e-01 3.76e-01w  1
    3171  0.0000000e+00 1.88e+05 4.51e+12  -1.0 3.99e+03    -  1.00e+00 4.95e-01w  1
    3172  0.0000000e+00 2.51e+04 6.16e+12  -1.0 2.20e+03    -  7.11e-01 5.57e-05h 13
    3173  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.27e+04    -  7.11e-01 5.57e-05h 14
    3174  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.11e-01 5.57e-05h 14
    3175  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.11e-01 5.57e-05h 14
    3176  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.11e-01 5.57e-05h 14
    3177  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.27e+04    -  7.12e-01 5.57e-05h 14
    3178  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.12e-01 5.57e-05h 14
    3179  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.12e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3180  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.26e+04    -  7.12e-01 5.57e-05h 14
    3181  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.12e-01 5.57e-05h 14
    3182  0.0000000e+00 3.94e+05 1.26e+13  -1.0 5.25e+04    -  7.12e-01 4.56e-01w  1
    3183  0.0000000e+00 3.61e+05 6.17e+12  -1.0 5.35e+04    -  7.31e-01 3.76e-01w  1
    3184  0.0000000e+00 1.87e+05 4.52e+12  -1.0 3.95e+03    -  1.00e+00 4.95e-01w  1
    3185  0.0000000e+00 2.51e+04 6.16e+12  -1.0 2.18e+03    -  7.12e-01 5.57e-05h 13
    3186  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.27e+04    -  7.12e-01 5.57e-05h 14
    3187  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.26e+04    -  7.12e-01 5.57e-05h 14
    3188  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.12e-01 5.57e-05h 14
    3189  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.27e+04    -  7.12e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3190  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.12e-01 5.57e-05h 14
    3191  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.12e-01 5.57e-05h 14
    3192  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.27e+04    -  7.12e-01 5.57e-05h 14
    3193  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.28e+04    -  7.12e-01 5.57e-05h 14
    3194  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.25e+04    -  7.13e-01 5.57e-05h 14
    3195  0.0000000e+00 3.93e+05 1.26e+13  -1.0 5.27e+04    -  7.13e-01 4.56e-01w  1
    3196  0.0000000e+00 3.60e+05 6.18e+12  -1.0 5.35e+04    -  7.32e-01 3.77e-01w  1
    3197  0.0000000e+00 1.87e+05 4.52e+12  -1.0 4.02e+03    -  1.00e+00 4.95e-01w  1
    3198  0.0000000e+00 2.51e+04 6.16e+12  -1.0 2.22e+03    -  7.13e-01 5.57e-05h 13
    3199  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.27e+04    -  7.13e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3200  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.27e+04    -  7.13e-01 5.57e-05h 14
    3201  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.27e+04    -  7.13e-01 5.57e-05h 14
    3202  0.0000000e+00 2.51e+04 6.16e+12  -1.0 5.25e+04    -  7.13e-01 5.57e-05h 14
    3203  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.27e+04    -  7.13e-01 5.57e-05h 14
    3204  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.27e+04    -  7.13e-01 5.57e-05h 14
    3205  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.27e+04    -  7.13e-01 5.57e-05h 14
    3206  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.26e+04    -  7.13e-01 5.57e-05h 14
    3207  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.25e+04    -  7.13e-01 5.57e-05h 14
    3208  0.0000000e+00 3.93e+05 1.26e+13  -1.0 5.27e+04    -  7.13e-01 4.56e-01w  1
    3209  0.0000000e+00 3.60e+05 6.17e+12  -1.0 5.35e+04    -  7.33e-01 3.77e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3210  0.0000000e+00 1.87e+05 4.53e+12  -1.0 4.06e+03    -  1.00e+00 4.95e-01w  1
    3211  0.0000000e+00 2.50e+04 6.16e+12  -1.0 2.24e+03    -  7.13e-01 5.57e-05h 13
    3212  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.27e+04    -  7.13e-01 5.57e-05h 14
    3213  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.27e+04    -  7.13e-01 5.57e-05h 14
    3214  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.27e+04    -  7.14e-01 5.57e-05h 14
    3215  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.27e+04    -  7.14e-01 5.57e-05h 14
    3216  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.27e+04    -  7.14e-01 5.57e-05h 14
    3217  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.26e+04    -  7.14e-01 5.57e-05h 14
    3218  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.27e+04    -  7.14e-01 5.57e-05h 14
    3219  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.27e+04    -  7.14e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3220  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.27e+04    -  7.14e-01 5.57e-05h 14
    3221  0.0000000e+00 3.93e+05 1.26e+13  -1.0 5.27e+04    -  7.14e-01 4.56e-01w  1
    3222  0.0000000e+00 3.59e+05 6.17e+12  -1.0 5.34e+04    -  7.33e-01 3.77e-01w  1
    3223  0.0000000e+00 1.86e+05 4.53e+12  -1.0 4.08e+03    -  1.00e+00 4.95e-01w  1
    3224  0.0000000e+00 2.50e+04 6.16e+12  -1.0 2.25e+03    -  7.14e-01 5.57e-05h 13
    3225  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.27e+04    -  7.14e-01 5.57e-05h 14
    3226  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.27e+04    -  7.14e-01 5.57e-05h 14
    3227  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.27e+04    -  7.14e-01 5.57e-05h 14
    3228  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.24e+04    -  7.14e-01 5.57e-05h 14
    3229  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.26e+04    -  7.15e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3230  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.27e+04    -  7.14e-01 5.57e-05h 14
    3231  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.26e+04    -  7.15e-01 5.57e-05h 14
    3232  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.26e+04    -  7.15e-01 5.57e-05h 14
    3233  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.27e+04    -  7.15e-01 5.57e-05h 14
    3234  0.0000000e+00 3.92e+05 1.26e+13  -1.0 5.26e+04    -  7.15e-01 4.56e-01w  1
    3235  0.0000000e+00 3.59e+05 6.17e+12  -1.0 5.33e+04    -  7.34e-01 3.77e-01w  1
    3236  0.0000000e+00 1.86e+05 4.53e+12  -1.0 4.10e+03    -  1.00e+00 4.95e-01w  1
    3237  0.0000000e+00 2.50e+04 6.16e+12  -1.0 2.25e+03    -  7.15e-01 5.57e-05h 13
    3238  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.27e+04    -  7.15e-01 5.57e-05h 14
    3239  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.26e+04    -  7.15e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3240  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.27e+04    -  7.15e-01 5.57e-05h 14
    3241  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.26e+04    -  7.15e-01 5.57e-05h 14
    3242  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.25e+04    -  7.15e-01 5.57e-05h 14
    3243  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.26e+04    -  7.15e-01 5.57e-05h 14
    3244  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.26e+04    -  7.15e-01 5.57e-05h 14
    3245  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.26e+04    -  7.15e-01 5.57e-05h 14
    3246  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.26e+04    -  7.15e-01 5.57e-05h 14
    3247  0.0000000e+00 3.92e+05 1.26e+13  -1.0 5.26e+04    -  7.15e-01 4.56e-01w  1
    3248  0.0000000e+00 3.58e+05 6.17e+12  -1.0 5.33e+04    -  7.34e-01 3.77e-01w  1
    3249  0.0000000e+00 1.85e+05 4.54e+12  -1.0 4.13e+03    -  1.00e+00 4.95e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3250  0.0000000e+00 2.50e+04 6.16e+12  -1.0 2.27e+03    -  7.15e-01 5.57e-05h 13
    3251  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.26e+04    -  7.16e-01 5.57e-05h 14
    3252  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.26e+04    -  7.16e-01 5.57e-05h 14
    3253  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.26e+04    -  7.16e-01 5.57e-05h 14
    3254  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.25e+04    -  7.16e-01 5.57e-05h 14
    3255  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.26e+04    -  7.16e-01 5.57e-05h 14
    3256  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.26e+04    -  7.16e-01 5.57e-05h 14
    3257  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.26e+04    -  7.16e-01 5.57e-05h 14
    3258  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.25e+04    -  7.16e-01 5.57e-05h 14
    3259  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.23e+04    -  7.16e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3260  0.0000000e+00 3.91e+05 1.26e+13  -1.0 5.26e+04    -  7.16e-01 4.56e-01w  1
    3261  0.0000000e+00 3.57e+05 6.17e+12  -1.0 5.32e+04    -  7.35e-01 3.78e-01w  1
    3262  0.0000000e+00 1.85e+05 4.54e+12  -1.0 4.16e+03    -  1.00e+00 4.95e-01w  1
    3263  0.0000000e+00 2.50e+04 6.16e+12  -1.0 2.29e+03    -  7.16e-01 5.57e-05h 13
    3264  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.26e+04    -  7.16e-01 5.57e-05h 14
    3265  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.26e+04    -  7.16e-01 5.57e-05h 14
    3266  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.25e+04    -  7.16e-01 5.57e-05h 14
    3267  0.0000000e+00 2.50e+04 6.16e+12  -1.0 5.26e+04    -  7.16e-01 5.57e-05h 14
    3268  0.0000000e+00 2.50e+04 6.17e+12  -1.0 5.24e+04    -  7.17e-01 5.57e-05h 14
    3269  0.0000000e+00 2.50e+04 6.17e+12  -1.0 5.25e+04    -  7.17e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3270  0.0000000e+00 2.50e+04 6.17e+12  -1.0 5.25e+04    -  7.17e-01 5.57e-05h 14
    3271  0.0000000e+00 2.50e+04 6.17e+12  -1.0 5.25e+04    -  7.17e-01 5.57e-05h 14
    3272  0.0000000e+00 2.50e+04 6.17e+12  -1.0 5.26e+04    -  7.17e-01 5.57e-05h 14
    3273  0.0000000e+00 3.91e+05 1.26e+13  -1.0 5.26e+04    -  7.17e-01 4.56e-01w  1
    3274  0.0000000e+00 3.57e+05 6.17e+12  -1.0 5.31e+04    -  7.35e-01 3.78e-01w  1
    3275  0.0000000e+00 1.85e+05 4.54e+12  -1.0 4.18e+03    -  1.00e+00 4.95e-01w  1
    3276  0.0000000e+00 2.50e+04 6.17e+12  -1.0 2.30e+03    -  7.17e-01 5.57e-05h 13
    3277  0.0000000e+00 2.50e+04 6.17e+12  -1.0 5.26e+04    -  7.17e-01 5.57e-05h 14
    3278  0.0000000e+00 2.50e+04 6.17e+12  -1.0 5.26e+04    -  7.17e-01 5.57e-05h 14
    3279  0.0000000e+00 2.50e+04 6.17e+12  -1.0 5.25e+04    -  7.17e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3280  0.0000000e+00 2.50e+04 6.17e+12  -1.0 5.26e+04    -  7.17e-01 5.57e-05h 14
    3281  0.0000000e+00 2.50e+04 6.17e+12  -1.0 5.25e+04    -  7.17e-01 5.57e-05h 14
    3282  0.0000000e+00 2.50e+04 6.17e+12  -1.0 5.26e+04    -  7.17e-01 5.57e-05h 14
    3283  0.0000000e+00 2.50e+04 6.17e+12  -1.0 5.26e+04    -  7.17e-01 5.57e-05h 14
    3284  0.0000000e+00 2.50e+04 6.17e+12  -1.0 5.24e+04    -  7.17e-01 5.57e-05h 14
    3285  0.0000000e+00 2.50e+04 6.17e+12  -1.0 5.25e+04    -  7.18e-01 5.57e-05h 14
    3286  0.0000000e+00 3.90e+05 1.26e+13  -1.0 5.24e+04    -  7.18e-01 4.56e-01w  1
    3287  0.0000000e+00 3.56e+05 6.17e+12  -1.0 5.30e+04    -  7.36e-01 3.78e-01w  1
    3288  0.0000000e+00 1.84e+05 4.55e+12  -1.0 4.15e+03    -  1.00e+00 4.95e-01w  1
    3289  0.0000000e+00 2.50e+04 6.17e+12  -1.0 2.28e+03    -  7.18e-01 5.57e-05h 13
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3290  0.0000000e+00 2.50e+04 6.17e+12  -1.0 5.25e+04    -  7.18e-01 5.57e-05h 14
    3291  0.0000000e+00 2.50e+04 6.17e+12  -1.0 5.26e+04    -  7.18e-01 5.57e-05h 14
    3292  0.0000000e+00 2.50e+04 6.17e+12  -1.0 5.25e+04    -  7.18e-01 5.57e-05h 14
    3293  0.0000000e+00 2.50e+04 6.17e+12  -1.0 5.25e+04    -  7.18e-01 5.57e-05h 14
    3294  0.0000000e+00 2.50e+04 6.17e+12  -1.0 5.25e+04    -  7.18e-01 5.57e-05h 14
    3295  0.0000000e+00 2.50e+04 6.17e+12  -1.0 5.25e+04    -  7.18e-01 5.57e-05h 14
    3296  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.23e+04    -  7.18e-01 5.57e-05h 14
    3297  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.18e-01 5.57e-05h 14
    3298  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.25e+04    -  7.18e-01 5.57e-05h 14
    3299  0.0000000e+00 3.90e+05 1.26e+13  -1.0 5.25e+04    -  7.18e-01 4.56e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3300  0.0000000e+00 3.56e+05 6.17e+12  -1.0 5.30e+04    -  7.37e-01 3.78e-01w  1
    3301  0.0000000e+00 1.84e+05 4.55e+12  -1.0 4.23e+03    -  1.00e+00 4.95e-01w  1
    3302  0.0000000e+00 2.49e+04 6.17e+12  -1.0 2.32e+03    -  7.18e-01 5.57e-05h 13
    3303  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.25e+04    -  7.18e-01 5.57e-05h 14
    3304  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.25e+04    -  7.18e-01 5.57e-05h 14
    3305  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.19e-01 5.57e-05h 14
    3306  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.23e+04    -  7.19e-01 5.57e-05h 14
    3307  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.25e+04    -  7.19e-01 5.57e-05h 14
    3308  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.25e+04    -  7.19e-01 5.57e-05h 14
    3309  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.25e+04    -  7.19e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3310  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.19e-01 5.57e-05h 14
    3311  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.19e-01 5.57e-05h 14
    3312  0.0000000e+00 3.89e+05 1.26e+13  -1.0 5.25e+04    -  7.19e-01 4.56e-01w  1
    3313  0.0000000e+00 3.55e+05 6.17e+12  -1.0 5.29e+04    -  7.37e-01 3.79e-01w  1
    3314  0.0000000e+00 1.84e+05 4.56e+12  -1.0 4.25e+03    -  1.00e+00 4.95e-01w  1
    3315  0.0000000e+00 2.49e+04 6.17e+12  -1.0 2.34e+03    -  7.19e-01 5.57e-05h 13
    3316  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.25e+04    -  7.19e-01 5.57e-05h 14
    3317  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.25e+04    -  7.19e-01 5.57e-05h 14
    3318  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.25e+04    -  7.19e-01 5.57e-05h 14
    3319  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.25e+04    -  7.19e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3320  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.25e+04    -  7.19e-01 5.57e-05h 14
    3321  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.25e+04    -  7.19e-01 5.57e-05h 14
    3322  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.25e+04    -  7.20e-01 5.57e-05h 14
    3323  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.25e+04    -  7.20e-01 5.57e-05h 14
    3324  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.23e+04    -  7.20e-01 5.57e-05h 14
    3325  0.0000000e+00 3.89e+05 1.26e+13  -1.0 5.25e+04    -  7.20e-01 4.56e-01w  1
    3326  0.0000000e+00 3.54e+05 6.17e+12  -1.0 5.29e+04    -  7.38e-01 3.79e-01w  1
    3327  0.0000000e+00 1.83e+05 4.56e+12  -1.0 4.28e+03    -  1.00e+00 4.95e-01w  1
    3328  0.0000000e+00 2.49e+04 6.17e+12  -1.0 2.35e+03    -  7.20e-01 5.57e-05h 13
    3329  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.25e+04    -  7.20e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3330  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.25e+04    -  7.20e-01 5.57e-05h 14
    3331  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.25e+04    -  7.20e-01 5.57e-05h 14
    3332  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.25e+04    -  7.20e-01 5.57e-05h 14
    3333  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.25e+04    -  7.20e-01 5.57e-05h 14
    3334  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.20e-01 5.57e-05h 14
    3335  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.20e-01 5.57e-05h 14
    3336  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.20e-01 5.57e-05h 14
    3337  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.25e+04    -  7.20e-01 5.57e-05h 14
    3338  0.0000000e+00 3.89e+05 1.26e+13  -1.0 5.22e+04    -  7.20e-01 4.56e-01w  1
    3339  0.0000000e+00 3.54e+05 6.16e+12  -1.0 5.27e+04    -  7.38e-01 3.79e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3340  0.0000000e+00 1.83e+05 4.56e+12  -1.0 4.23e+03    -  1.00e+00 4.95e-01w  1
    3341  0.0000000e+00 2.49e+04 6.17e+12  -1.0 2.32e+03    -  7.20e-01 5.57e-05h 13
    3342  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.21e-01 5.57e-05h 14
    3343  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.23e+04    -  7.21e-01 5.57e-05h 14
    3344  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.21e-01 5.57e-05h 14
    3345  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.21e-01 5.57e-05h 14
    3346  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.23e+04    -  7.21e-01 5.57e-05h 14
    3347  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.21e-01 5.57e-05h 14
    3348  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.21e-01 5.57e-05h 14
    3349  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.21e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3350  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.23e+04    -  7.21e-01 5.57e-05h 14
    3351  0.0000000e+00 3.88e+05 1.26e+13  -1.0 5.24e+04    -  7.21e-01 4.56e-01w  1
    3352  0.0000000e+00 3.53e+05 6.16e+12  -1.0 5.27e+04    -  7.39e-01 3.79e-01w  1
    3353  0.0000000e+00 1.83e+05 4.57e+12  -1.0 4.32e+03    -  1.00e+00 4.95e-01w  1
    3354  0.0000000e+00 2.49e+04 6.17e+12  -1.0 2.37e+03    -  7.21e-01 5.57e-05h 13
    3355  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.21e-01 5.57e-05h 14
    3356  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.21e-01 5.57e-05h 14
    3357  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.21e-01 5.57e-05h 14
    3358  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.21e-01 5.57e-05h 14
    3359  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.22e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3360  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.23e+04    -  7.22e-01 5.57e-05h 14
    3361  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.22e+04    -  7.22e-01 5.57e-05h 14
    3362  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.22e-01 5.57e-05h 14
    3363  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.22e-01 5.57e-05h 14
    3364  0.0000000e+00 3.88e+05 1.26e+13  -1.0 5.24e+04    -  7.22e-01 4.56e-01w  1
    3365  0.0000000e+00 3.53e+05 6.16e+12  -1.0 5.26e+04    -  7.39e-01 3.80e-01w  1
    3366  0.0000000e+00 1.82e+05 4.57e+12  -1.0 4.34e+03    -  1.00e+00 4.95e-01w  1
    3367  0.0000000e+00 2.49e+04 6.17e+12  -1.0 2.38e+03    -  7.22e-01 5.57e-05h 13
    3368  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.22e-01 5.57e-05h 14
    3369  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.23e+04    -  7.22e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3370  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.22e-01 5.57e-05h 14
    3371  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.23e+04    -  7.22e-01 5.57e-05h 14
    3372  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.23e+04    -  7.22e-01 5.57e-05h 14
    3373  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.22e-01 5.57e-05h 14
    3374  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.23e+04    -  7.22e-01 5.57e-05h 14
    3375  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.23e+04    -  7.22e-01 5.57e-05h 14
    3376  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.23e-01 5.57e-05h 14
    3377  0.0000000e+00 3.87e+05 1.26e+13  -1.0 5.23e+04    -  7.23e-01 4.56e-01w  1
    3378  0.0000000e+00 3.52e+05 6.16e+12  -1.0 5.26e+04    -  7.40e-01 3.80e-01w  1
    3379  0.0000000e+00 1.82e+05 4.57e+12  -1.0 4.37e+03    -  1.00e+00 4.95e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3380  0.0000000e+00 2.49e+04 6.17e+12  -1.0 2.39e+03    -  7.23e-01 5.57e-05h 13
    3381  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.23e-01 5.57e-05h 14
    3382  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.24e+04    -  7.23e-01 5.57e-05h 14
    3383  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.23e+04    -  7.23e-01 5.57e-05h 14
    3384  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.22e+04    -  7.23e-01 5.57e-05h 14
    3385  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.23e+04    -  7.23e-01 5.57e-05h 14
    3386  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.23e+04    -  7.23e-01 5.57e-05h 14
    3387  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.21e+04    -  7.23e-01 5.57e-05h 14
    3388  0.0000000e+00 2.49e+04 6.17e+12  -1.0 5.23e+04    -  7.23e-01 5.57e-05h 14
    3389  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.23e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3390  0.0000000e+00 3.87e+05 1.26e+13  -1.0 5.23e+04    -  7.23e-01 4.56e-01w  1
    3391  0.0000000e+00 3.51e+05 6.16e+12  -1.0 5.25e+04    -  7.40e-01 3.80e-01w  1
    3392  0.0000000e+00 1.82e+05 4.58e+12  -1.0 4.38e+03    -  1.00e+00 4.95e-01w  1
    3393  0.0000000e+00 2.48e+04 6.17e+12  -1.0 2.40e+03    -  7.23e-01 5.57e-05h 13
    3394  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.23e+04    -  7.23e-01 5.57e-05h 14
    3395  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.23e+04    -  7.23e-01 5.57e-05h 14
    3396  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.24e-01 5.57e-05h 14
    3397  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.23e+04    -  7.24e-01 5.57e-05h 14
    3398  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.23e+04    -  7.24e-01 5.57e-05h 14
    3399  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.23e+04    -  7.24e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3400  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.23e+04    -  7.24e-01 5.57e-05h 14
    3401  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.23e+04    -  7.24e-01 5.57e-05h 14
    3402  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.24e-01 5.57e-05h 14
    3403  0.0000000e+00 3.86e+05 1.26e+13  -1.0 5.23e+04    -  7.24e-01 4.56e-01w  1
    3404  0.0000000e+00 3.51e+05 6.16e+12  -1.0 5.25e+04    -  7.41e-01 3.80e-01w  1
    3405  0.0000000e+00 1.81e+05 4.58e+12  -1.0 4.42e+03    -  1.00e+00 4.95e-01w  1
    3406  0.0000000e+00 2.48e+04 6.17e+12  -1.0 2.42e+03    -  7.24e-01 5.57e-05h 13
    3407  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.23e+04    -  7.24e-01 5.57e-05h 14
    3408  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.24e-01 5.57e-05h 14
    3409  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.23e+04    -  7.24e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3410  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.23e+04    -  7.24e-01 5.57e-05h 14
    3411  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.20e+04    -  7.24e-01 5.57e-05h 14
    3412  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.24e-01 5.57e-05h 14
    3413  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.23e+04    -  7.25e-01 5.57e-05h 14
    3414  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.23e+04    -  7.25e-01 5.57e-05h 14
    3415  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.25e-01 5.57e-05h 14
    3416  0.0000000e+00 3.86e+05 1.26e+13  -1.0 5.21e+04    -  7.25e-01 4.56e-01w  1
    3417  0.0000000e+00 3.50e+05 6.16e+12  -1.0 5.23e+04    -  7.41e-01 3.80e-01w  1
    3418  0.0000000e+00 1.81e+05 4.59e+12  -1.0 4.38e+03    -  1.00e+00 4.95e-01w  1
    3419  0.0000000e+00 2.48e+04 6.17e+12  -1.0 2.40e+03    -  7.25e-01 5.57e-05h 13
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3420  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.23e+04    -  7.25e-01 5.57e-05h 14
    3421  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.23e+04    -  7.25e-01 5.57e-05h 14
    3422  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.23e+04    -  7.25e-01 5.57e-05h 14
    3423  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.23e+04    -  7.25e-01 5.57e-05h 14
    3424  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.25e-01 5.57e-05h 14
    3425  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.23e+04    -  7.25e-01 5.57e-05h 14
    3426  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.25e-01 5.57e-05h 14
    3427  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.25e-01 5.57e-05h 14
    3428  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.21e+04    -  7.25e-01 5.57e-05h 14
    3429  0.0000000e+00 3.85e+05 1.26e+13  -1.0 5.23e+04    -  7.25e-01 4.56e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3430  0.0000000e+00 3.50e+05 6.16e+12  -1.0 5.23e+04    -  7.42e-01 3.81e-01w  1
    3431  0.0000000e+00 1.81e+05 4.59e+12  -1.0 4.47e+03    -  1.00e+00 4.95e-01w  1
    3432  0.0000000e+00 2.48e+04 6.17e+12  -1.0 2.45e+03    -  7.25e-01 5.57e-05h 13
    3433  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.26e-01 5.57e-05h 14
    3434  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.26e-01 5.57e-05h 14
    3435  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.26e-01 5.57e-05h 14
    3436  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.26e-01 5.57e-05h 14
    3437  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.26e-01 5.57e-05h 14
    3438  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.26e-01 5.57e-05h 14
    3439  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.26e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3440  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.21e+04    -  7.26e-01 5.57e-05h 14
    3441  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.26e-01 5.57e-05h 14
    3442  0.0000000e+00 3.85e+05 1.26e+13  -1.0 5.22e+04    -  7.26e-01 4.56e-01w  1
    3443  0.0000000e+00 3.49e+05 6.16e+12  -1.0 5.22e+04    -  7.43e-01 3.81e-01w  1
    3444  0.0000000e+00 1.80e+05 4.59e+12  -1.0 4.48e+03    -  1.00e+00 4.95e-01w  1
    3445  0.0000000e+00 2.48e+04 6.17e+12  -1.0 2.45e+03    -  7.26e-01 5.57e-05h 13
    3446  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.21e+04    -  7.26e-01 5.57e-05h 14
    3447  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.21e+04    -  7.26e-01 5.57e-05h 14
    3448  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.26e-01 5.57e-05h 14
    3449  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.26e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3450  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.27e-01 5.57e-05h 14
    3451  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.27e-01 5.57e-05h 14
    3452  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.27e-01 5.57e-05h 14
    3453  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.21e+04    -  7.27e-01 5.57e-05h 14
    3454  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.21e+04    -  7.27e-01 5.57e-05h 14
    3455  0.0000000e+00 3.84e+05 1.26e+13  -1.0 5.22e+04    -  7.27e-01 4.56e-01w  1
    3456  0.0000000e+00 3.48e+05 6.16e+12  -1.0 5.22e+04    -  7.43e-01 3.81e-01w  1
    3457  0.0000000e+00 1.80e+05 4.60e+12  -1.0 4.51e+03    -  1.00e+00 4.95e-01w  1
    3458  0.0000000e+00 2.48e+04 6.17e+12  -1.0 2.47e+03    -  7.27e-01 5.57e-05h 13
    3459  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.21e+04    -  7.27e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3460  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.21e+04    -  7.27e-01 5.57e-05h 14
    3461  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.27e-01 5.57e-05h 14
    3462  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.21e+04    -  7.27e-01 5.57e-05h 14
    3463  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.21e+04    -  7.27e-01 5.57e-05h 14
    3464  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.21e+04    -  7.27e-01 5.57e-05h 14
    3465  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.21e+04    -  7.27e-01 5.57e-05h 14
    3466  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.27e-01 5.57e-05h 14
    3467  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.21e+04    -  7.27e-01 5.57e-05h 14
    3468  0.0000000e+00 3.84e+05 1.26e+13  -1.0 5.21e+04    -  7.28e-01 4.56e-01w  1
    3469  0.0000000e+00 3.48e+05 6.16e+12  -1.0 5.21e+04    -  7.44e-01 3.81e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3470  0.0000000e+00 1.80e+05 4.60e+12  -1.0 4.53e+03    -  1.00e+00 4.95e-01w  1
    3471  0.0000000e+00 2.48e+04 6.17e+12  -1.0 2.48e+03    -  7.28e-01 5.57e-05h 13
    3472  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.28e-01 5.57e-05h 14
    3473  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.28e-01 5.57e-05h 14
    3474  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.22e+04    -  7.28e-01 5.57e-05h 14
    3475  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.21e+04    -  7.28e-01 5.57e-05h 14
    3476  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.20e+04    -  7.28e-01 5.57e-05h 14
    3477  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.21e+04    -  7.28e-01 5.57e-05h 14
    3478  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.20e+04    -  7.28e-01 5.57e-05h 14
    3479  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.21e+04    -  7.28e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3480  0.0000000e+00 2.48e+04 6.17e+12  -1.0 5.21e+04    -  7.28e-01 5.57e-05h 14
    3481  0.0000000e+00 3.84e+05 1.26e+13  -1.0 5.21e+04    -  7.28e-01 4.56e-01w  1
    3482  0.0000000e+00 3.47e+05 6.16e+12  -1.0 5.21e+04    -  7.44e-01 3.82e-01w  1
    3483  0.0000000e+00 1.79e+05 4.60e+12  -1.0 4.55e+03    -  1.00e+00 4.95e-01w  1
    3484  0.0000000e+00 2.48e+04 6.17e+12  -1.0 2.49e+03    -  7.28e-01 5.57e-05h 13
    3485  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.28e-01 5.57e-05h 14
    3486  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.28e-01 5.57e-05h 14
    3487  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.21e+04    -  7.29e-01 5.57e-05h 14
    3488  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.19e+04    -  7.29e-01 5.57e-05h 14
    3489  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.18e+04    -  7.29e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3490  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.29e-01 5.57e-05h 14
    3491  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.29e-01 5.57e-05h 14
    3492  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.18e+04    -  7.29e-01 5.57e-05h 14
    3493  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.19e+04    -  7.29e-01 5.57e-05h 14
    3494  0.0000000e+00 3.83e+05 1.26e+13  -1.0 5.20e+04    -  7.29e-01 4.56e-01w  1
    3495  0.0000000e+00 3.47e+05 6.15e+12  -1.0 5.20e+04    -  7.45e-01 3.82e-01w  1
    3496  0.0000000e+00 1.79e+05 4.61e+12  -1.0 4.56e+03    -  1.00e+00 4.95e-01w  1
    3497  0.0000000e+00 2.47e+04 6.17e+12  -1.0 2.49e+03    -  7.29e-01 5.57e-05h 13
    3498  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.21e+04    -  7.29e-01 5.57e-05h 14
    3499  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.21e+04    -  7.29e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3500  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.21e+04    -  7.29e-01 5.57e-05h 14
    3501  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.21e+04    -  7.29e-01 5.57e-05h 14
    3502  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.21e+04    -  7.29e-01 5.57e-05h 14
    3503  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.21e+04    -  7.29e-01 5.57e-05h 14
    3504  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.16e+04    -  7.30e-01 5.57e-05h 14
    3505  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.30e-01 5.57e-05h 14
    3506  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.30e-01 5.57e-05h 14
    3507  0.0000000e+00 3.83e+05 1.26e+13  -1.0 5.21e+04    -  7.30e-01 4.56e-01w  1
    3508  0.0000000e+00 3.46e+05 6.15e+12  -1.0 5.19e+04    -  7.45e-01 3.82e-01w  1
    3509  0.0000000e+00 1.79e+05 4.61e+12  -1.0 4.60e+03    -  1.00e+00 4.95e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3510  0.0000000e+00 2.47e+04 6.17e+12  -1.0 2.51e+03    -  7.30e-01 5.57e-05h 13
    3511  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.30e-01 5.57e-05h 14
    3512  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.21e+04    -  7.30e-01 5.57e-05h 14
    3513  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.17e+04    -  7.30e-01 5.57e-05h 14
    3514  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.30e-01 5.57e-05h 14
    3515  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.30e-01 5.57e-05h 14
    3516  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.21e+04    -  7.30e-01 5.57e-05h 14
    3517  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.21e+04    -  7.30e-01 5.57e-05h 14
    3518  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.19e+04    -  7.30e-01 5.57e-05h 14
    3519  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.21e+04    -  7.30e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3520  0.0000000e+00 3.82e+05 1.26e+13  -1.0 5.21e+04    -  7.30e-01 4.56e-01w  1
    3521  0.0000000e+00 3.46e+05 6.15e+12  -1.0 5.19e+04    -  7.46e-01 3.82e-01w  1
    3522  0.0000000e+00 1.78e+05 4.62e+12  -1.0 4.63e+03    -  1.00e+00 4.95e-01w  1
    3523  0.0000000e+00 2.47e+04 6.17e+12  -1.0 2.53e+03    -  7.30e-01 5.57e-05h 13
    3524  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.21e+04    -  7.30e-01 5.57e-05h 14
    3525  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.19e+04    -  7.30e-01 5.57e-05h 14
    3526  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.31e-01 5.57e-05h 14
    3527  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.31e-01 5.57e-05h 14
    3528  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.31e-01 5.57e-05h 14
    3529  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.31e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3530  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.17e+04    -  7.31e-01 5.57e-05h 14
    3531  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.31e-01 5.57e-05h 14
    3532  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.31e-01 5.57e-05h 14
    3533  0.0000000e+00 3.82e+05 1.26e+13  -1.0 5.20e+04    -  7.31e-01 4.57e-01w  1
    3534  0.0000000e+00 3.45e+05 6.15e+12  -1.0 5.18e+04    -  7.47e-01 3.82e-01w  1
    3535  0.0000000e+00 1.78e+05 4.62e+12  -1.0 4.63e+03    -  1.00e+00 4.95e-01w  1
    3536  0.0000000e+00 2.47e+04 6.17e+12  -1.0 2.53e+03    -  7.31e-01 5.57e-05h 13
    3537  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.19e+04    -  7.31e-01 5.57e-05h 14
    3538  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.31e-01 5.57e-05h 14
    3539  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.31e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3540  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.31e-01 5.57e-05h 14
    3541  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.19e+04    -  7.31e-01 5.57e-05h 14
    3542  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.31e-01 5.57e-05h 14
    3543  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.32e-01 5.57e-05h 14
    3544  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.19e+04    -  7.32e-01 5.57e-05h 14
    3545  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.32e-01 5.57e-05h 14
    3546  0.0000000e+00 3.81e+05 1.26e+13  -1.0 5.20e+04    -  7.32e-01 4.57e-01w  1
    3547  0.0000000e+00 3.44e+05 6.15e+12  -1.0 5.17e+04    -  7.47e-01 3.83e-01w  1
    3548  0.0000000e+00 1.78e+05 4.62e+12  -1.0 4.67e+03    -  1.00e+00 4.95e-01w  1
    3549  0.0000000e+00 2.47e+04 6.17e+12  -1.0 2.55e+03    -  7.32e-01 5.57e-05h 13
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3550  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.19e+04    -  7.32e-01 5.57e-05h 14
    3551  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.32e-01 5.57e-05h 14
    3552  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.32e-01 5.57e-05h 14
    3553  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.32e-01 5.57e-05h 14
    3554  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.32e-01 5.57e-05h 14
    3555  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.32e-01 5.57e-05h 14
    3556  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.32e-01 5.57e-05h 14
    3557  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.18e+04    -  7.32e-01 5.57e-05h 14
    3558  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.19e+04    -  7.32e-01 5.57e-05h 14
    3559  0.0000000e+00 3.81e+05 1.26e+13  -1.0 5.20e+04    -  7.32e-01 4.57e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3560  0.0000000e+00 3.44e+05 6.15e+12  -1.0 5.17e+04    -  7.48e-01 3.83e-01w  1
    3561  0.0000000e+00 1.77e+05 4.63e+12  -1.0 4.70e+03    -  1.00e+00 4.95e-01w  1
    3562  0.0000000e+00 2.47e+04 6.17e+12  -1.0 2.56e+03    -  7.32e-01 5.57e-05h 13
    3563  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.33e-01 5.57e-05h 14
    3564  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.33e-01 5.57e-05h 14
    3565  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.19e+04    -  7.33e-01 5.57e-05h 14
    3566  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.33e-01 5.57e-05h 14
    3567  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.19e+04    -  7.33e-01 5.57e-05h 14
    3568  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.19e+04    -  7.33e-01 5.57e-05h 14
    3569  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.20e+04    -  7.33e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3570  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.19e+04    -  7.33e-01 5.57e-05h 14
    3571  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.19e+04    -  7.33e-01 5.57e-05h 14
    3572  0.0000000e+00 3.80e+05 1.26e+13  -1.0 5.14e+04    -  7.33e-01 4.57e-01w  1
    3573  0.0000000e+00 3.43e+05 6.15e+12  -1.0 5.14e+04    -  7.48e-01 3.83e-01w  1
    3574  0.0000000e+00 1.77e+05 4.63e+12  -1.0 4.57e+03    -  1.00e+00 4.95e-01w  1
    3575  0.0000000e+00 2.47e+04 6.17e+12  -1.0 2.50e+03    -  7.33e-01 5.57e-05h 13
    3576  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.19e+04    -  7.33e-01 5.57e-05h 14
    3577  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.19e+04    -  7.33e-01 5.57e-05h 14
    3578  0.0000000e+00 2.47e+04 6.17e+12  -1.0 5.19e+04    -  7.33e-01 5.57e-05h 14
    3579  0.0000000e+00 2.46e+04 6.17e+12  -1.0 5.19e+04    -  7.33e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3580  0.0000000e+00 2.46e+04 6.17e+12  -1.0 5.19e+04    -  7.34e-01 5.57e-05h 14
    3581  0.0000000e+00 2.46e+04 6.17e+12  -1.0 5.18e+04    -  7.34e-01 5.57e-05h 14
    3582  0.0000000e+00 2.46e+04 6.17e+12  -1.0 5.19e+04    -  7.34e-01 5.57e-05h 14
    3583  0.0000000e+00 2.46e+04 6.17e+12  -1.0 5.18e+04    -  7.34e-01 5.57e-05h 14
    3584  0.0000000e+00 2.46e+04 6.17e+12  -1.0 5.17e+04    -  7.34e-01 5.57e-05h 14
    3585  0.0000000e+00 3.80e+05 1.26e+13  -1.0 5.19e+04    -  7.34e-01 4.57e-01w  1
    3586  0.0000000e+00 3.43e+05 6.15e+12  -1.0 5.15e+04    -  7.49e-01 3.83e-01w  1
    3587  0.0000000e+00 1.77e+05 4.64e+12  -1.0 4.73e+03    -  1.00e+00 4.95e-01w  1
    3588  0.0000000e+00 2.46e+04 6.17e+12  -1.0 2.58e+03    -  7.34e-01 5.57e-05h 13
    3589  0.0000000e+00 2.46e+04 6.17e+12  -1.0 5.19e+04    -  7.34e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3590  0.0000000e+00 2.46e+04 6.17e+12  -1.0 5.19e+04    -  7.34e-01 5.57e-05h 14
    3591  0.0000000e+00 2.46e+04 6.17e+12  -1.0 5.18e+04    -  7.34e-01 5.57e-05h 14
    3592  0.0000000e+00 2.46e+04 6.17e+12  -1.0 5.17e+04    -  7.34e-01 5.57e-05h 14
    3593  0.0000000e+00 2.46e+04 6.17e+12  -1.0 5.18e+04    -  7.34e-01 5.57e-05h 14
    3594  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.19e+04    -  7.34e-01 5.57e-05h 14
    3595  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.18e+04    -  7.34e-01 5.57e-05h 14
    3596  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.18e+04    -  7.34e-01 5.57e-05h 14
    3597  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.19e+04    -  7.35e-01 5.57e-05h 14
    3598  0.0000000e+00 3.80e+05 1.26e+13  -1.0 5.19e+04    -  7.35e-01 4.57e-01w  1
    3599  0.0000000e+00 3.42e+05 6.15e+12  -1.0 5.15e+04    -  7.49e-01 3.84e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3600  0.0000000e+00 1.76e+05 4.64e+12  -1.0 4.77e+03    -  1.00e+00 4.95e-01w  1
    3601  0.0000000e+00 2.46e+04 6.18e+12  -1.0 2.60e+03    -  7.35e-01 5.57e-05h 13
    3602  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.19e+04    -  7.35e-01 5.57e-05h 14
    3603  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.18e+04    -  7.35e-01 5.57e-05h 14
    3604  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.19e+04    -  7.35e-01 5.57e-05h 14
    3605  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.18e+04    -  7.35e-01 5.57e-05h 14
    3606  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.16e+04    -  7.35e-01 5.57e-05h 14
    3607  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.18e+04    -  7.35e-01 5.57e-05h 14
    3608  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.19e+04    -  7.35e-01 5.57e-05h 14
    3609  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.19e+04    -  7.35e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3610  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.19e+04    -  7.35e-01 5.57e-05h 14
    3611  0.0000000e+00 3.79e+05 1.26e+13  -1.0 5.18e+04    -  7.35e-01 4.57e-01w  1
    3612  0.0000000e+00 3.42e+05 6.15e+12  -1.0 5.14e+04    -  7.50e-01 3.84e-01w  1
    3613  0.0000000e+00 1.76e+05 4.64e+12  -1.0 4.77e+03    -  1.00e+00 4.95e-01w  1
    3614  0.0000000e+00 2.46e+04 6.18e+12  -1.0 2.60e+03    -  7.35e-01 5.57e-05h 13
    3615  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.19e+04    -  7.35e-01 5.57e-05h 14
    3616  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.19e+04    -  7.35e-01 5.57e-05h 14
    3617  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.19e+04    -  7.35e-01 5.57e-05h 14
    3618  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.18e+04    -  7.36e-01 5.57e-05h 14
    3619  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.19e+04    -  7.36e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3620  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.18e+04    -  7.36e-01 5.57e-05h 14
    3621  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.18e+04    -  7.36e-01 5.57e-05h 14
    3622  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.18e+04    -  7.36e-01 5.57e-05h 14
    3623  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.18e+04    -  7.36e-01 5.57e-05h 14
    3624  0.0000000e+00 3.79e+05 1.26e+13  -1.0 5.18e+04    -  7.36e-01 4.57e-01w  1
    3625  0.0000000e+00 3.41e+05 6.15e+12  -1.0 5.14e+04    -  7.50e-01 3.84e-01w  1
    3626  0.0000000e+00 1.76e+05 4.65e+12  -1.0 4.81e+03    -  1.00e+00 4.95e-01w  1
    3627  0.0000000e+00 2.46e+04 6.18e+12  -1.0 2.62e+03    -  7.36e-01 5.57e-05h 13
    3628  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.15e+04    -  7.36e-01 5.58e-05h 14
    3629  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.18e+04    -  7.36e-01 5.57e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3630  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.17e+04    -  7.36e-01 5.57e-05h 14
    3631  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.17e+04    -  7.36e-01 5.57e-05h 14
    3632  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.18e+04    -  7.36e-01 5.57e-05h 14
    3633  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.18e+04    -  7.36e-01 5.57e-05h 14
    3634  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.18e+04    -  7.36e-01 5.57e-05h 14
    3635  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.18e+04    -  7.37e-01 5.57e-05h 14
    3636  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.17e+04    -  7.37e-01 5.57e-05h 14
    3637  0.0000000e+00 3.78e+05 1.26e+13  -1.0 5.18e+04    -  7.37e-01 4.57e-01w  1
    3638  0.0000000e+00 3.40e+05 6.15e+12  -1.0 5.13e+04    -  7.51e-01 3.84e-01w  1
    3639  0.0000000e+00 1.75e+05 4.65e+12  -1.0 4.84e+03    -  1.00e+00 4.95e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3640  0.0000000e+00 2.46e+04 6.18e+12  -1.0 2.63e+03    -  7.37e-01 5.57e-05h 13
    3641  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.17e+04    -  7.37e-01 5.58e-05h 14
    3642  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.18e+04    -  7.37e-01 5.57e-05h 14
    3643  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.18e+04    -  7.37e-01 5.57e-05h 14
    3644  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.18e+04    -  7.37e-01 5.57e-05h 14
    3645  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.18e+04    -  7.37e-01 5.57e-05h 14
    3646  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.18e+04    -  7.37e-01 5.57e-05h 14
    3647  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.17e+04    -  7.37e-01 5.58e-05h 14
    3648  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.16e+04    -  7.37e-01 5.58e-05h 14
    3649  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.18e+04    -  7.37e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3650  0.0000000e+00 3.78e+05 1.26e+13  -1.0 5.17e+04    -  7.37e-01 4.57e-01w  1
    3651  0.0000000e+00 3.40e+05 6.15e+12  -1.0 5.12e+04    -  7.51e-01 3.84e-01w  1
    3652  0.0000000e+00 1.75e+05 4.66e+12  -1.0 4.84e+03    -  1.00e+00 4.95e-01w  1
    3653  0.0000000e+00 2.46e+04 6.18e+12  -1.0 2.64e+03    -  7.37e-01 5.58e-05h 13
    3654  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.17e+04    -  7.37e-01 5.58e-05h 14
    3655  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.18e+04    -  7.38e-01 5.58e-05h 14
    3656  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.17e+04    -  7.38e-01 5.58e-05h 14
    3657  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.18e+04    -  7.38e-01 5.58e-05h 14
    3658  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.17e+04    -  7.38e-01 5.58e-05h 14
    3659  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.18e+04    -  7.38e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3660  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.17e+04    -  7.38e-01 5.58e-05h 14
    3661  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.17e+04    -  7.38e-01 5.58e-05h 14
    3662  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.17e+04    -  7.38e-01 5.58e-05h 14
    3663  0.0000000e+00 3.77e+05 1.26e+13  -1.0 5.18e+04    -  7.38e-01 4.57e-01w  1
    3664  0.0000000e+00 3.39e+05 6.15e+12  -1.0 5.12e+04    -  7.52e-01 3.85e-01w  1
    3665  0.0000000e+00 1.75e+05 4.66e+12  -1.0 4.88e+03    -  1.00e+00 4.95e-01w  1
    3666  0.0000000e+00 2.46e+04 6.18e+12  -1.0 2.66e+03    -  7.38e-01 5.58e-05h 13
    3667  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.17e+04    -  7.38e-01 5.58e-05h 14
    3668  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.17e+04    -  7.38e-01 5.58e-05h 14
    3669  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.16e+04    -  7.38e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3670  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.17e+04    -  7.38e-01 5.58e-05h 14
    3671  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.16e+04    -  7.38e-01 5.58e-05h 14
    3672  0.0000000e+00 2.46e+04 6.18e+12  -1.0 5.17e+04    -  7.39e-01 5.58e-05h 14
    3673  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.17e+04    -  7.39e-01 5.58e-05h 14
    3674  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.17e+04    -  7.39e-01 5.58e-05h 14
    3675  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.17e+04    -  7.39e-01 5.58e-05h 14
    3676  0.0000000e+00 3.77e+05 1.26e+13  -1.0 5.17e+04    -  7.39e-01 4.57e-01w  1
    3677  0.0000000e+00 3.39e+05 6.15e+12  -1.0 5.11e+04    -  7.52e-01 3.85e-01w  1
    3678  0.0000000e+00 1.74e+05 4.66e+12  -1.0 4.89e+03    -  1.00e+00 4.95e-01w  1
    3679  0.0000000e+00 2.45e+04 6.18e+12  -1.0 2.66e+03    -  7.39e-01 5.58e-05h 13
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3680  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.17e+04    -  7.39e-01 5.58e-05h 14
    3681  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.17e+04    -  7.39e-01 5.58e-05h 14
    3682  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.17e+04    -  7.39e-01 5.58e-05h 14
    3683  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.17e+04    -  7.39e-01 5.58e-05h 14
    3684  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.17e+04    -  7.39e-01 5.58e-05h 14
    3685  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.17e+04    -  7.39e-01 5.58e-05h 14
    3686  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.17e+04    -  7.39e-01 5.58e-05h 14
    3687  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.17e+04    -  7.39e-01 5.58e-05h 14
    3688  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.14e+04    -  7.39e-01 5.58e-05h 14
    3689  0.0000000e+00 3.76e+05 1.26e+13  -1.0 5.17e+04    -  7.40e-01 4.57e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3690  0.0000000e+00 3.38e+05 6.15e+12  -1.0 5.11e+04    -  7.53e-01 3.85e-01w  1
    3691  0.0000000e+00 1.74e+05 4.67e+12  -1.0 4.92e+03    -  1.00e+00 4.95e-01w  1
    3692  0.0000000e+00 2.45e+04 6.18e+12  -1.0 2.68e+03    -  7.40e-01 5.58e-05h 13
    3693  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.15e+04    -  7.40e-01 5.58e-05h 14
    3694  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.16e+04    -  7.40e-01 5.58e-05h 14
    3695  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.17e+04    -  7.40e-01 5.58e-05h 14
    3696  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.17e+04    -  7.40e-01 5.58e-05h 14
    3697  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.16e+04    -  7.40e-01 5.58e-05h 14
    3698  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.15e+04    -  7.40e-01 5.58e-05h 14
    3699  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.15e+04    -  7.40e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3700  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.17e+04    -  7.40e-01 5.58e-05h 14
    3701  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.17e+04    -  7.40e-01 5.58e-05h 14
    3702  0.0000000e+00 3.76e+05 1.26e+13  -1.0 5.17e+04    -  7.40e-01 4.57e-01w  1
    3703  0.0000000e+00 3.38e+05 6.15e+12  -1.0 5.10e+04    -  7.54e-01 3.85e-01w  1
    3704  0.0000000e+00 1.74e+05 4.67e+12  -1.0 4.94e+03    -  1.00e+00 4.95e-01w  1
    3705  0.0000000e+00 2.45e+04 6.18e+12  -1.0 2.69e+03    -  7.40e-01 5.58e-05h 13
    3706  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.15e+04    -  7.40e-01 5.58e-05h 14
    3707  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.17e+04    -  7.40e-01 5.58e-05h 14
    3708  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.16e+04    -  7.40e-01 5.58e-05h 14
    3709  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.17e+04    -  7.40e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3710  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.17e+04    -  7.41e-01 5.58e-05h 14
    3711  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.17e+04    -  7.41e-01 5.58e-05h 14
    3712  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.17e+04    -  7.41e-01 5.58e-05h 14
    3713  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.17e+04    -  7.41e-01 5.58e-05h 14
    3714  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.15e+04    -  7.41e-01 5.58e-05h 14
    3715  0.0000000e+00 3.75e+05 1.26e+13  -1.0 5.16e+04    -  7.41e-01 4.57e-01w  1
    3716  0.0000000e+00 3.37e+05 6.15e+12  -1.0 5.09e+04    -  7.54e-01 3.86e-01w  1
    3717  0.0000000e+00 1.73e+05 4.67e+12  -1.0 4.96e+03    -  1.00e+00 4.95e-01w  1
    3718  0.0000000e+00 2.45e+04 6.18e+12  -1.0 2.70e+03    -  7.41e-01 5.58e-05h 13
    3719  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.16e+04    -  7.41e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3720  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.16e+04    -  7.41e-01 5.58e-05h 14
    3721  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.16e+04    -  7.41e-01 5.58e-05h 14
    3722  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.16e+04    -  7.41e-01 5.58e-05h 14
    3723  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.16e+04    -  7.41e-01 5.58e-05h 14
    3724  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.15e+04    -  7.41e-01 5.58e-05h 14
    3725  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.15e+04    -  7.41e-01 5.58e-05h 14
    3726  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.16e+04    -  7.41e-01 5.58e-05h 14
    3727  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.16e+04    -  7.42e-01 5.58e-05h 14
    3728  0.0000000e+00 3.75e+05 1.26e+13  -1.0 5.16e+04    -  7.42e-01 4.57e-01w  1
    3729  0.0000000e+00 3.36e+05 6.15e+12  -1.0 5.09e+04    -  7.55e-01 3.86e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3730  0.0000000e+00 1.73e+05 4.68e+12  -1.0 4.99e+03    -  1.00e+00 4.95e-01w  1
    3731  0.0000000e+00 2.45e+04 6.18e+12  -1.0 2.71e+03    -  7.42e-01 5.58e-05h 13
    3732  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.16e+04    -  7.42e-01 5.58e-05h 14
    3733  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.16e+04    -  7.42e-01 5.58e-05h 14
    3734  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.16e+04    -  7.42e-01 5.58e-05h 14
    3735  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.15e+04    -  7.42e-01 5.58e-05h 14
    3736  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.15e+04    -  7.42e-01 5.58e-05h 14
    3737  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.16e+04    -  7.42e-01 5.58e-05h 14
    3738  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.16e+04    -  7.42e-01 5.58e-05h 14
    3739  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.14e+04    -  7.42e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3740  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.15e+04    -  7.42e-01 5.58e-05h 14
    3741  0.0000000e+00 3.75e+05 1.26e+13  -1.0 5.15e+04    -  7.42e-01 4.57e-01w  1
    3742  0.0000000e+00 3.36e+05 6.14e+12  -1.0 5.08e+04    -  7.55e-01 3.86e-01w  1
    3743  0.0000000e+00 1.73e+05 4.68e+12  -1.0 4.99e+03    -  1.00e+00 4.95e-01w  1
    3744  0.0000000e+00 2.45e+04 6.18e+12  -1.0 2.71e+03    -  7.42e-01 5.58e-05h 13
    3745  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.16e+04    -  7.42e-01 5.58e-05h 14
    3746  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.15e+04    -  7.42e-01 5.58e-05h 14
    3747  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.15e+04    -  7.43e-01 5.58e-05h 14
    3748  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.16e+04    -  7.43e-01 5.58e-05h 14
    3749  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.16e+04    -  7.43e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3750  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.16e+04    -  7.43e-01 5.58e-05h 14
    3751  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.15e+04    -  7.43e-01 5.58e-05h 14
    3752  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.16e+04    -  7.43e-01 5.58e-05h 14
    3753  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.16e+04    -  7.43e-01 5.58e-05h 14
    3754  0.0000000e+00 3.74e+05 1.26e+13  -1.0 5.15e+04    -  7.43e-01 4.57e-01w  1
    3755  0.0000000e+00 3.35e+05 6.14e+12  -1.0 5.07e+04    -  7.56e-01 3.86e-01w  1
    3756  0.0000000e+00 1.73e+05 4.69e+12  -1.0 5.03e+03    -  1.00e+00 4.95e-01w  1
    3757  0.0000000e+00 2.45e+04 6.18e+12  -1.0 2.73e+03    -  7.43e-01 5.58e-05h 13
    3758  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.16e+04    -  7.43e-01 5.58e-05h 14
    3759  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.15e+04    -  7.43e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3760  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.15e+04    -  7.43e-01 5.58e-05h 14
    3761  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.16e+04    -  7.43e-01 5.58e-05h 14
    3762  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.16e+04    -  7.43e-01 5.58e-05h 14
    3763  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.15e+04    -  7.43e-01 5.58e-05h 14
    3764  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.15e+04    -  7.43e-01 5.58e-05h 14
    3765  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.15e+04    -  7.44e-01 5.58e-05h 14
    3766  0.0000000e+00 2.45e+04 6.18e+12  -1.0 5.15e+04    -  7.44e-01 5.58e-05h 14
    3767  0.0000000e+00 3.74e+05 1.26e+13  -1.0 5.15e+04    -  7.44e-01 4.57e-01w  1
    3768  0.0000000e+00 3.35e+05 6.14e+12  -1.0 5.06e+04    -  7.56e-01 3.86e-01w  1
    3769  0.0000000e+00 1.72e+05 4.69e+12  -1.0 5.03e+03    -  1.00e+00 4.95e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3770  0.0000000e+00 2.44e+04 6.18e+12  -1.0 2.73e+03    -  7.44e-01 5.58e-05h 13
    3771  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.15e+04    -  7.44e-01 5.58e-05h 14
    3772  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.12e+04    -  7.44e-01 5.58e-05h 14
    3773  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.15e+04    -  7.44e-01 5.58e-05h 14
    3774  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.15e+04    -  7.44e-01 5.58e-05h 14
    3775  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.15e+04    -  7.44e-01 5.58e-05h 14
    3776  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.15e+04    -  7.44e-01 5.58e-05h 14
    3777  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.15e+04    -  7.44e-01 5.58e-05h 14
    3778  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.44e-01 5.58e-05h 14
    3779  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.15e+04    -  7.44e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3780  0.0000000e+00 3.73e+05 1.26e+13  -1.0 5.15e+04    -  7.44e-01 4.57e-01w  1
    3781  0.0000000e+00 3.34e+05 6.14e+12  -1.0 5.06e+04    -  7.57e-01 3.87e-01w  1
    3782  0.0000000e+00 1.72e+05 4.69e+12  -1.0 5.07e+03    -  1.00e+00 4.95e-01w  1
    3783  0.0000000e+00 2.44e+04 6.18e+12  -1.0 2.75e+03    -  7.44e-01 5.58e-05h 13
    3784  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.15e+04    -  7.44e-01 5.58e-05h 14
    3785  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.15e+04    -  7.44e-01 5.58e-05h 14
    3786  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.15e+04    -  7.45e-01 5.58e-05h 14
    3787  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.15e+04    -  7.45e-01 5.58e-05h 14
    3788  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.45e-01 5.58e-05h 14
    3789  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.45e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3790  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.15e+04    -  7.45e-01 5.58e-05h 14
    3791  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.15e+04    -  7.45e-01 5.58e-05h 14
    3792  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.45e-01 5.58e-05h 14
    3793  0.0000000e+00 3.73e+05 1.26e+13  -1.0 5.15e+04    -  7.45e-01 4.57e-01w  1
    3794  0.0000000e+00 3.34e+05 6.14e+12  -1.0 5.05e+04    -  7.57e-01 3.87e-01w  1
    3795  0.0000000e+00 1.72e+05 4.70e+12  -1.0 5.10e+03    -  1.00e+00 4.95e-01w  1
    3796  0.0000000e+00 2.44e+04 6.18e+12  -1.0 2.77e+03    -  7.45e-01 5.58e-05h 13
    3797  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.15e+04    -  7.45e-01 5.58e-05h 14
    3798  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.45e-01 5.58e-05h 14
    3799  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.45e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3800  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.45e-01 5.58e-05h 14
    3801  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.15e+04    -  7.45e-01 5.58e-05h 14
    3802  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.45e-01 5.58e-05h 14
    3803  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.46e-01 5.58e-05h 14
    3804  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.15e+04    -  7.46e-01 5.58e-05h 14
    3805  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.13e+04    -  7.46e-01 5.58e-05h 14
    3806  0.0000000e+00 3.72e+05 1.26e+13  -1.0 5.14e+04    -  7.46e-01 4.57e-01w  1
    3807  0.0000000e+00 3.33e+05 6.14e+12  -1.0 5.05e+04    -  7.58e-01 3.87e-01w  1
    3808  0.0000000e+00 1.71e+05 4.70e+12  -1.0 5.12e+03    -  1.00e+00 4.95e-01w  1
    3809  0.0000000e+00 2.44e+04 6.18e+12  -1.0 2.78e+03    -  7.46e-01 5.58e-05h 13
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3810  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.46e-01 5.58e-05h 14
    3811  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.46e-01 5.58e-05h 14
    3812  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.46e-01 5.58e-05h 14
    3813  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.46e-01 5.58e-05h 14
    3814  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.46e-01 5.58e-05h 14
    3815  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.13e+04    -  7.46e-01 5.58e-05h 14
    3816  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.46e-01 5.58e-05h 14
    3817  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.46e-01 5.58e-05h 14
    3818  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.46e-01 5.58e-05h 14
    3819  0.0000000e+00 3.72e+05 1.26e+13  -1.0 5.14e+04    -  7.46e-01 4.57e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3820  0.0000000e+00 3.33e+05 6.14e+12  -1.0 5.04e+04    -  7.58e-01 3.87e-01w  1
    3821  0.0000000e+00 1.71e+05 4.71e+12  -1.0 5.14e+03    -  1.00e+00 4.95e-01w  1
    3822  0.0000000e+00 2.44e+04 6.18e+12  -1.0 2.79e+03    -  7.46e-01 5.58e-05h 13
    3823  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.13e+04    -  7.47e-01 5.58e-05h 14
    3824  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.47e-01 5.58e-05h 14
    3825  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.47e-01 5.58e-05h 14
    3826  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.47e-01 5.58e-05h 14
    3827  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.13e+04    -  7.47e-01 5.58e-05h 14
    3828  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.13e+04    -  7.47e-01 5.58e-05h 14
    3829  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.47e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3830  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.47e-01 5.58e-05h 14
    3831  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.13e+04    -  7.47e-01 5.58e-05h 14
    3832  0.0000000e+00 3.72e+05 1.26e+13  -1.0 5.14e+04    -  7.47e-01 4.57e-01w  1
    3833  0.0000000e+00 3.32e+05 6.14e+12  -1.0 5.04e+04    -  7.59e-01 3.87e-01w  1
    3834  0.0000000e+00 1.71e+05 4.71e+12  -1.0 5.16e+03    -  1.00e+00 4.95e-01w  1
    3835  0.0000000e+00 2.44e+04 6.18e+12  -1.0 2.80e+03    -  7.47e-01 5.58e-05h 13
    3836  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.47e-01 5.58e-05h 14
    3837  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.12e+04    -  7.47e-01 5.58e-05h 14
    3838  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.47e-01 5.58e-05h 14
    3839  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.12e+04    -  7.47e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3840  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.13e+04    -  7.48e-01 5.58e-05h 14
    3841  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.12e+04    -  7.48e-01 5.58e-05h 14
    3842  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.48e-01 5.58e-05h 14
    3843  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.48e-01 5.58e-05h 14
    3844  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.48e-01 5.58e-05h 14
    3845  0.0000000e+00 3.71e+05 1.26e+13  -1.0 5.12e+04    -  7.48e-01 4.57e-01w  1
    3846  0.0000000e+00 3.31e+05 6.14e+12  -1.0 5.03e+04    -  7.59e-01 3.88e-01w  1
    3847  0.0000000e+00 1.70e+05 4.71e+12  -1.0 5.14e+03    -  1.00e+00 4.95e-01w  1
    3848  0.0000000e+00 2.44e+04 6.18e+12  -1.0 2.79e+03    -  7.48e-01 5.58e-05h 13
    3849  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.14e+04    -  7.48e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3850  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.13e+04    -  7.48e-01 5.58e-05h 14
    3851  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.13e+04    -  7.48e-01 5.58e-05h 14
    3852  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.13e+04    -  7.48e-01 5.58e-05h 14
    3853  0.0000000e+00 2.44e+04 6.18e+12  -1.0 5.13e+04    -  7.48e-01 5.58e-05h 14
    3854  0.0000000e+00 2.44e+04 6.19e+12  -1.0 5.13e+04    -  7.48e-01 5.58e-05h 14
    3855  0.0000000e+00 2.44e+04 6.19e+12  -1.0 5.13e+04    -  7.48e-01 5.58e-05h 14
    3856  0.0000000e+00 2.44e+04 6.19e+12  -1.0 5.13e+04    -  7.48e-01 5.58e-05h 14
    3857  0.0000000e+00 2.44e+04 6.19e+12  -1.0 5.12e+04    -  7.48e-01 5.58e-05h 14
    3858  0.0000000e+00 3.71e+05 1.26e+13  -1.0 5.13e+04    -  7.49e-01 4.57e-01w  1
    3859  0.0000000e+00 3.31e+05 6.14e+12  -1.0 5.02e+04    -  7.60e-01 3.88e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3860  0.0000000e+00 1.70e+05 4.72e+12  -1.0 5.20e+03    -  1.00e+00 4.95e-01w  1
    3861  0.0000000e+00 2.44e+04 6.19e+12  -1.0 2.82e+03    -  7.49e-01 5.58e-05h 13
    3862  0.0000000e+00 2.44e+04 6.19e+12  -1.0 5.12e+04    -  7.49e-01 5.58e-05h 14
    3863  0.0000000e+00 2.44e+04 6.19e+12  -1.0 5.13e+04    -  7.49e-01 5.58e-05h 14
    3864  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.13e+04    -  7.49e-01 5.58e-05h 14
    3865  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.11e+04    -  7.49e-01 5.58e-05h 14
    3866  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.13e+04    -  7.49e-01 5.58e-05h 14
    3867  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.13e+04    -  7.49e-01 5.58e-05h 14
    3868  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.13e+04    -  7.49e-01 5.58e-05h 14
    3869  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.13e+04    -  7.49e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3870  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.13e+04    -  7.49e-01 5.58e-05h 14
    3871  0.0000000e+00 3.70e+05 1.26e+13  -1.0 5.13e+04    -  7.49e-01 4.57e-01w  1
    3872  0.0000000e+00 3.30e+05 6.14e+12  -1.0 5.02e+04    -  7.60e-01 3.88e-01w  1
    3873  0.0000000e+00 1.70e+05 4.72e+12  -1.0 5.22e+03    -  1.00e+00 4.95e-01w  1
    3874  0.0000000e+00 2.43e+04 6.19e+12  -1.0 2.83e+03    -  7.49e-01 5.58e-05h 13
    3875  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.13e+04    -  7.49e-01 5.58e-05h 14
    3876  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.13e+04    -  7.49e-01 5.58e-05h 14
    3877  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.13e+04    -  7.49e-01 5.58e-05h 14
    3878  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.13e+04    -  7.49e-01 5.58e-05h 14
    3879  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.13e+04    -  7.50e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3880  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.13e+04    -  7.50e-01 5.58e-05h 14
    3881  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.13e+04    -  7.50e-01 5.58e-05h 14
    3882  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.50e-01 5.58e-05h 14
    3883  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.50e-01 5.58e-05h 14
    3884  0.0000000e+00 3.70e+05 1.26e+13  -1.0 5.12e+04    -  7.50e-01 4.57e-01w  1
    3885  0.0000000e+00 3.30e+05 6.14e+12  -1.0 5.01e+04    -  7.61e-01 3.88e-01w  1
    3886  0.0000000e+00 1.69e+05 4.73e+12  -1.0 5.22e+03    -  1.00e+00 4.95e-01w  1
    3887  0.0000000e+00 2.43e+04 6.19e+12  -1.0 2.83e+03    -  7.50e-01 5.58e-05h 13
    3888  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.13e+04    -  7.50e-01 5.58e-05h 14
    3889  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.50e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3890  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.50e-01 5.58e-05h 14
    3891  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.13e+04    -  7.50e-01 5.58e-05h 14
    3892  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.13e+04    -  7.50e-01 5.58e-05h 14
    3893  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.13e+04    -  7.50e-01 5.58e-05h 14
    3894  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.11e+04    -  7.50e-01 5.58e-05h 14
    3895  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.51e-01 5.58e-05h 14
    3896  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.51e-01 5.58e-05h 14
    3897  0.0000000e+00 3.69e+05 1.26e+13  -1.0 5.13e+04    -  7.51e-01 4.57e-01w  1
    3898  0.0000000e+00 3.29e+05 6.14e+12  -1.0 5.01e+04    -  7.62e-01 3.88e-01w  1
    3899  0.0000000e+00 1.69e+05 4.73e+12  -1.0 5.27e+03    -  1.00e+00 4.95e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3900  0.0000000e+00 2.43e+04 6.19e+12  -1.0 2.85e+03    -  7.51e-01 5.58e-05h 13
    3901  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.11e+04    -  7.51e-01 5.58e-05h 14
    3902  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.11e+04    -  7.51e-01 5.58e-05h 14
    3903  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.51e-01 5.58e-05h 14
    3904  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.51e-01 5.58e-05h 14
    3905  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.51e-01 5.58e-05h 14
    3906  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.11e+04    -  7.51e-01 5.58e-05h 14
    3907  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.51e-01 5.58e-05h 14
    3908  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.51e-01 5.58e-05h 14
    3909  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.51e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3910  0.0000000e+00 3.69e+05 1.26e+13  -1.0 5.12e+04    -  7.51e-01 4.57e-01w  1
    3911  0.0000000e+00 3.29e+05 6.14e+12  -1.0 5.00e+04    -  7.62e-01 3.89e-01w  1
    3912  0.0000000e+00 1.69e+05 4.73e+12  -1.0 5.28e+03    -  1.00e+00 4.95e-01w  1
    3913  0.0000000e+00 2.43e+04 6.19e+12  -1.0 2.86e+03    -  7.51e-01 5.58e-05h 13
    3914  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.51e-01 5.58e-05h 14
    3915  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.51e-01 5.58e-05h 14
    3916  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.52e-01 5.58e-05h 14
    3917  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.10e+04    -  7.52e-01 5.58e-05h 14
    3918  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.52e-01 5.58e-05h 14
    3919  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.52e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3920  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.52e-01 5.58e-05h 14
    3921  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.52e-01 5.58e-05h 14
    3922  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.11e+04    -  7.52e-01 5.58e-05h 14
    3923  0.0000000e+00 3.68e+05 1.26e+13  -1.0 5.12e+04    -  7.52e-01 4.57e-01w  1
    3924  0.0000000e+00 3.28e+05 6.14e+12  -1.0 5.00e+04    -  7.63e-01 3.89e-01w  1
    3925  0.0000000e+00 1.69e+05 4.74e+12  -1.0 5.31e+03    -  1.00e+00 4.95e-01w  1
    3926  0.0000000e+00 2.43e+04 6.19e+12  -1.0 2.88e+03    -  7.52e-01 5.58e-05h 13
    3927  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.52e-01 5.58e-05h 14
    3928  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.11e+04    -  7.52e-01 5.58e-05h 14
    3929  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.52e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3930  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.52e-01 5.58e-05h 14
    3931  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.52e-01 5.58e-05h 14
    3932  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.11e+04    -  7.52e-01 5.58e-05h 14
    3933  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.53e-01 5.58e-05h 14
    3934  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.53e-01 5.58e-05h 14
    3935  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.11e+04    -  7.53e-01 5.58e-05h 14
    3936  0.0000000e+00 3.68e+05 1.26e+13  -1.0 5.11e+04    -  7.53e-01 4.57e-01w  1
    3937  0.0000000e+00 3.28e+05 6.14e+12  -1.0 4.99e+04    -  7.63e-01 3.89e-01w  1
    3938  0.0000000e+00 1.68e+05 4.74e+12  -1.0 5.30e+03    -  1.00e+00 4.95e-01w  1
    3939  0.0000000e+00 2.43e+04 6.19e+12  -1.0 2.88e+03    -  7.53e-01 5.58e-05h 13
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3940  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.11e+04    -  7.53e-01 5.58e-05h 14
    3941  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.53e-01 5.58e-05h 14
    3942  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.10e+04    -  7.53e-01 5.58e-05h 14
    3943  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.11e+04    -  7.53e-01 5.58e-05h 14
    3944  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.11e+04    -  7.53e-01 5.58e-05h 14
    3945  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.11e+04    -  7.53e-01 5.58e-05h 14
    3946  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.53e-01 5.58e-05h 14
    3947  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.53e-01 5.58e-05h 14
    3948  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.11e+04    -  7.53e-01 5.58e-05h 14
    3949  0.0000000e+00 3.68e+05 1.26e+13  -1.0 5.11e+04    -  7.53e-01 4.57e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3950  0.0000000e+00 3.27e+05 6.14e+12  -1.0 4.98e+04    -  7.64e-01 3.89e-01w  1
    3951  0.0000000e+00 1.68e+05 4.75e+12  -1.0 5.32e+03    -  1.00e+00 4.95e-01w  1
    3952  0.0000000e+00 2.43e+04 6.19e+12  -1.0 2.88e+03    -  7.53e-01 5.58e-05h 13
    3953  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.11e+04    -  7.53e-01 5.58e-05h 14
    3954  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.12e+04    -  7.54e-01 5.58e-05h 14
    3955  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.11e+04    -  7.54e-01 5.58e-05h 14
    3956  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.11e+04    -  7.54e-01 5.58e-05h 14
    3957  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.10e+04    -  7.54e-01 5.58e-05h 14
    3958  0.0000000e+00 2.43e+04 6.19e+12  -1.0 5.11e+04    -  7.54e-01 5.58e-05h 14
    3959  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.54e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3960  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.11e+04    -  7.54e-01 5.58e-05h 14
    3961  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.54e-01 5.58e-05h 14
    3962  0.0000000e+00 3.67e+05 1.26e+13  -1.0 5.11e+04    -  7.54e-01 4.57e-01w  1
    3963  0.0000000e+00 3.27e+05 6.14e+12  -1.0 4.98e+04    -  7.64e-01 3.89e-01w  1
    3964  0.0000000e+00 1.68e+05 4.75e+12  -1.0 5.36e+03    -  1.00e+00 4.95e-01w  1
    3965  0.0000000e+00 2.42e+04 6.19e+12  -1.0 2.90e+03    -  7.54e-01 5.58e-05h 13
    3966  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.11e+04    -  7.54e-01 5.58e-05h 14
    3967  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.11e+04    -  7.54e-01 5.58e-05h 14
    3968  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.11e+04    -  7.54e-01 5.58e-05h 14
    3969  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.07e+04    -  7.54e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3970  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.11e+04    -  7.55e-01 5.58e-05h 14
    3971  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.11e+04    -  7.54e-01 5.58e-05h 14
    3972  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.11e+04    -  7.55e-01 5.58e-05h 14
    3973  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.11e+04    -  7.55e-01 5.58e-05h 14
    3974  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.11e+04    -  7.55e-01 5.58e-05h 14
    3975  0.0000000e+00 3.67e+05 1.26e+13  -1.0 5.11e+04    -  7.55e-01 4.57e-01w  1
    3976  0.0000000e+00 3.26e+05 6.14e+12  -1.0 4.97e+04    -  7.65e-01 3.90e-01w  1
    3977  0.0000000e+00 1.67e+05 4.75e+12  -1.0 5.39e+03    -  1.00e+00 4.95e-01w  1
    3978  0.0000000e+00 2.42e+04 6.19e+12  -1.0 2.92e+03    -  7.55e-01 5.58e-05h 13
    3979  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.11e+04    -  7.55e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3980  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.11e+04    -  7.55e-01 5.58e-05h 14
    3981  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.55e-01 5.58e-05h 14
    3982  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.09e+04    -  7.55e-01 5.58e-05h 14
    3983  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.11e+04    -  7.55e-01 5.58e-05h 14
    3984  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.11e+04    -  7.55e-01 5.58e-05h 14
    3985  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.55e-01 5.58e-05h 14
    3986  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.11e+04    -  7.55e-01 5.58e-05h 14
    3987  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.11e+04    -  7.55e-01 5.58e-05h 14
    3988  0.0000000e+00 3.66e+05 1.26e+13  -1.0 5.11e+04    -  7.55e-01 4.57e-01w  1
    3989  0.0000000e+00 3.26e+05 6.14e+12  -1.0 4.97e+04    -  7.65e-01 3.90e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    3990  0.0000000e+00 1.67e+05 4.76e+12  -1.0 5.41e+03    -  1.00e+00 4.95e-01w  1
    3991  0.0000000e+00 2.42e+04 6.19e+12  -1.0 2.93e+03    -  7.55e-01 5.58e-05h 13
    3992  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.11e+04    -  7.56e-01 5.58e-05h 14
    3993  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.11e+04    -  7.56e-01 5.58e-05h 14
    3994  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.56e-01 5.58e-05h 14
    3995  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.56e-01 5.58e-05h 14
    3996  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.11e+04    -  7.56e-01 5.58e-05h 14
    3997  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.56e-01 5.58e-05h 14
    3998  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.11e+04    -  7.56e-01 5.58e-05h 14
    3999  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.56e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4000  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.11e+04    -  7.56e-01 5.58e-05h 14
    4001  0.0000000e+00 3.66e+05 1.26e+13  -1.0 5.11e+04    -  7.56e-01 4.57e-01w  1
    4002  0.0000000e+00 3.25e+05 6.14e+12  -1.0 4.96e+04    -  7.66e-01 3.90e-01w  1
    4003  0.0000000e+00 1.67e+05 4.76e+12  -1.0 5.43e+03    -  1.00e+00 4.95e-01w  1
    4004  0.0000000e+00 2.42e+04 6.19e+12  -1.0 2.94e+03    -  7.56e-01 5.58e-05h 13
    4005  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.56e-01 5.58e-05h 14
    4006  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.09e+04    -  7.56e-01 5.58e-05h 14
    4007  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.56e-01 5.58e-05h 14
    4008  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.56e-01 5.58e-05h 14
    4009  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.57e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4010  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.57e-01 5.58e-05h 14
    4011  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.57e-01 5.58e-05h 14
    4012  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.57e-01 5.58e-05h 14
    4013  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.57e-01 5.58e-05h 14
    4014  0.0000000e+00 3.65e+05 1.26e+13  -1.0 5.10e+04    -  7.57e-01 4.57e-01w  1
    4015  0.0000000e+00 3.24e+05 6.14e+12  -1.0 4.95e+04    -  7.66e-01 3.90e-01w  1
    4016  0.0000000e+00 1.66e+05 4.77e+12  -1.0 5.44e+03    -  1.00e+00 4.95e-01w  1
    4017  0.0000000e+00 2.42e+04 6.19e+12  -1.0 2.95e+03    -  7.57e-01 5.58e-05h 13
    4018  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.57e-01 5.58e-05h 14
    4019  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.57e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4020  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.08e+04    -  7.57e-01 5.58e-05h 14
    4021  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.57e-01 5.58e-05h 14
    4022  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.57e-01 5.58e-05h 14
    4023  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.57e-01 5.58e-05h 14
    4024  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.57e-01 5.58e-05h 14
    4025  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.09e+04    -  7.57e-01 5.58e-05h 14
    4026  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.58e-01 5.58e-05h 14
    4027  0.0000000e+00 3.65e+05 1.26e+13  -1.0 5.09e+04    -  7.58e-01 4.57e-01w  1
    4028  0.0000000e+00 3.24e+05 6.14e+12  -1.0 4.94e+04    -  7.67e-01 3.90e-01w  1
    4029  0.0000000e+00 1.66e+05 4.77e+12  -1.0 5.45e+03    -  1.00e+00 4.95e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4030  0.0000000e+00 2.42e+04 6.19e+12  -1.0 2.94e+03    -  7.58e-01 5.58e-05h 13
    4031  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.58e-01 5.58e-05h 14
    4032  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.58e-01 5.58e-05h 14
    4033  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.09e+04    -  7.58e-01 5.58e-05h 14
    4034  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.58e-01 5.58e-05h 14
    4035  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.09e+04    -  7.58e-01 5.58e-05h 14
    4036  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.09e+04    -  7.58e-01 5.58e-05h 14
    4037  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.09e+04    -  7.58e-01 5.58e-05h 14
    4038  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.09e+04    -  7.58e-01 5.58e-05h 14
    4039  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.58e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4040  0.0000000e+00 3.65e+05 1.26e+13  -1.0 5.07e+04    -  7.58e-01 4.57e-01w  1
    4041  0.0000000e+00 3.23e+05 6.14e+12  -1.0 4.93e+04    -  7.67e-01 3.91e-01w  1
    4042  0.0000000e+00 1.66e+05 4.77e+12  -1.0 5.42e+03    -  1.00e+00 4.95e-01w  1
    4043  0.0000000e+00 2.42e+04 6.19e+12  -1.0 2.93e+03    -  7.58e-01 5.58e-05h 13
    4044  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.09e+04    -  7.58e-01 5.58e-05h 14
    4045  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.09e+04    -  7.58e-01 5.58e-05h 14
    4046  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.09e+04    -  7.58e-01 5.58e-05h 14
    4047  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.09e+04    -  7.59e-01 5.58e-05h 14
    4048  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.59e-01 5.58e-05h 14
    4049  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.10e+04    -  7.59e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4050  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.09e+04    -  7.59e-01 5.58e-05h 14
    4051  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.09e+04    -  7.59e-01 5.58e-05h 14
    4052  0.0000000e+00 2.42e+04 6.19e+12  -1.0 5.09e+04    -  7.59e-01 5.58e-05h 14
    4053  0.0000000e+00 3.64e+05 1.26e+13  -1.0 5.09e+04    -  7.59e-01 4.57e-01w  1
    4054  0.0000000e+00 3.23e+05 6.14e+12  -1.0 4.94e+04    -  7.68e-01 3.91e-01w  1
    4055  0.0000000e+00 1.66e+05 4.78e+12  -1.0 5.52e+03    -  1.00e+00 4.95e-01w  1
    4056  0.0000000e+00 2.42e+04 6.19e+12  -1.0 2.98e+03    -  7.59e-01 5.58e-05h 13
    4057  0.0000000e+00 2.41e+04 6.19e+12  -1.0 5.09e+04    -  7.59e-01 5.58e-05h 14
    4058  0.0000000e+00 2.41e+04 6.19e+12  -1.0 5.08e+04    -  7.59e-01 5.58e-05h 14
    4059  0.0000000e+00 2.41e+04 6.19e+12  -1.0 5.09e+04    -  7.59e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4060  0.0000000e+00 2.41e+04 6.19e+12  -1.0 5.09e+04    -  7.59e-01 5.58e-05h 14
    4061  0.0000000e+00 2.41e+04 6.19e+12  -1.0 5.09e+04    -  7.59e-01 5.58e-05h 14
    4062  0.0000000e+00 2.41e+04 6.19e+12  -1.0 5.09e+04    -  7.59e-01 5.58e-05h 14
    4063  0.0000000e+00 2.41e+04 6.19e+12  -1.0 5.09e+04    -  7.59e-01 5.58e-05h 14
    4064  0.0000000e+00 2.41e+04 6.19e+12  -1.0 5.09e+04    -  7.59e-01 5.58e-05h 14
    4065  0.0000000e+00 2.41e+04 6.19e+12  -1.0 5.09e+04    -  7.60e-01 5.58e-05h 14
    4066  0.0000000e+00 3.64e+05 1.26e+13  -1.0 5.09e+04    -  7.60e-01 4.57e-01w  1
    4067  0.0000000e+00 3.22e+05 6.14e+12  -1.0 4.93e+04    -  7.68e-01 3.91e-01w  1
    4068  0.0000000e+00 1.65e+05 4.78e+12  -1.0 5.54e+03    -  1.00e+00 4.95e-01w  1
    4069  0.0000000e+00 2.41e+04 6.19e+12  -1.0 2.99e+03    -  7.60e-01 5.58e-05h 13
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4070  0.0000000e+00 2.41e+04 6.19e+12  -1.0 5.09e+04    -  7.60e-01 5.58e-05h 14
    4071  0.0000000e+00 2.41e+04 6.19e+12  -1.0 5.09e+04    -  7.60e-01 5.58e-05h 14
    4072  0.0000000e+00 2.41e+04 6.19e+12  -1.0 5.09e+04    -  7.60e-01 5.58e-05h 14
    4073  0.0000000e+00 2.41e+04 6.19e+12  -1.0 5.09e+04    -  7.60e-01 5.58e-05h 14
    4074  0.0000000e+00 2.41e+04 6.19e+12  -1.0 5.08e+04    -  7.60e-01 5.58e-05h 14
    4075  0.0000000e+00 2.41e+04 6.19e+12  -1.0 5.09e+04    -  7.60e-01 5.58e-05h 14
    4076  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.60e-01 5.58e-05h 14
    4077  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.09e+04    -  7.60e-01 5.58e-05h 14
    4078  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.09e+04    -  7.60e-01 5.58e-05h 14
    4079  0.0000000e+00 3.63e+05 1.26e+13  -1.0 5.08e+04    -  7.60e-01 4.57e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4080  0.0000000e+00 3.22e+05 6.14e+12  -1.0 4.92e+04    -  7.69e-01 3.91e-01w  1
    4081  0.0000000e+00 1.65e+05 4.79e+12  -1.0 5.53e+03    -  1.00e+00 4.95e-01w  1
    4082  0.0000000e+00 2.41e+04 6.20e+12  -1.0 2.99e+03    -  7.60e-01 5.58e-05h 13
    4083  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.09e+04    -  7.60e-01 5.58e-05h 14
    4084  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.60e-01 5.58e-05h 14
    4085  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.09e+04    -  7.61e-01 5.58e-05h 14
    4086  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.09e+04    -  7.61e-01 5.58e-05h 14
    4087  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.09e+04    -  7.61e-01 5.58e-05h 14
    4088  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.09e+04    -  7.61e-01 5.58e-05h 14
    4089  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.09e+04    -  7.61e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4090  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.61e-01 5.58e-05h 14
    4091  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.09e+04    -  7.61e-01 5.58e-05h 14
    4092  0.0000000e+00 3.63e+05 1.27e+13  -1.0 5.09e+04    -  7.61e-01 4.57e-01w  1
    4093  0.0000000e+00 3.21e+05 6.14e+12  -1.0 4.92e+04    -  7.69e-01 3.92e-01w  1
    4094  0.0000000e+00 1.65e+05 4.79e+12  -1.0 5.58e+03    -  1.00e+00 4.95e-01w  1
    4095  0.0000000e+00 2.41e+04 6.20e+12  -1.0 3.01e+03    -  7.61e-01 5.58e-05h 13
    4096  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.61e-01 5.58e-05h 14
    4097  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.09e+04    -  7.61e-01 5.58e-05h 14
    4098  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.61e-01 5.58e-05h 14
    4099  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.61e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4100  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.61e-01 5.58e-05h 14
    4101  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.61e-01 5.58e-05h 14
    4102  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.05e+04    -  7.62e-01 5.58e-05h 14
    4103  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.62e-01 5.58e-05h 14
    4104  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.62e-01 5.58e-05h 14
    4105  0.0000000e+00 3.62e+05 1.27e+13  -1.0 5.08e+04    -  7.62e-01 4.57e-01w  1
    4106  0.0000000e+00 3.21e+05 6.14e+12  -1.0 4.91e+04    -  7.70e-01 3.92e-01w  1
    4107  0.0000000e+00 1.64e+05 4.79e+12  -1.0 5.60e+03    -  1.00e+00 4.95e-01w  1
    4108  0.0000000e+00 2.41e+04 6.20e+12  -1.0 3.02e+03    -  7.62e-01 5.58e-05h 13
    4109  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.62e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4110  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.62e-01 5.58e-05h 14
    4111  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.62e-01 5.58e-05h 14
    4112  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.62e-01 5.58e-05h 14
    4113  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.62e-01 5.58e-05h 14
    4114  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.62e-01 5.58e-05h 14
    4115  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.62e-01 5.58e-05h 14
    4116  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.62e-01 5.58e-05h 14
    4117  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.62e-01 5.58e-05h 14
    4118  0.0000000e+00 3.62e+05 1.27e+13  -1.0 5.08e+04    -  7.62e-01 4.57e-01w  1
    4119  0.0000000e+00 3.20e+05 6.14e+12  -1.0 4.91e+04    -  7.70e-01 3.92e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4120  0.0000000e+00 1.64e+05 4.80e+12  -1.0 5.61e+03    -  1.00e+00 4.95e-01w  1
    4121  0.0000000e+00 2.41e+04 6.20e+12  -1.0 3.03e+03    -  7.62e-01 5.58e-05h 13
    4122  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.62e-01 5.58e-05h 14
    4123  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.62e-01 5.58e-05h 14
    4124  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.07e+04    -  7.63e-01 5.58e-05h 14
    4125  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.63e-01 5.58e-05h 14
    4126  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.63e-01 5.58e-05h 14
    4127  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.07e+04    -  7.63e-01 5.58e-05h 14
    4128  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.63e-01 5.58e-05h 14
    4129  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.63e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4130  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.63e-01 5.58e-05h 14
    4131  0.0000000e+00 3.62e+05 1.27e+13  -1.0 5.07e+04    -  7.63e-01 4.58e-01w  1
    4132  0.0000000e+00 3.20e+05 6.14e+12  -1.0 4.90e+04    -  7.71e-01 3.92e-01w  1
    4133  0.0000000e+00 1.64e+05 4.80e+12  -1.0 5.61e+03    -  1.00e+00 4.95e-01w  1
    4134  0.0000000e+00 2.41e+04 6.20e+12  -1.0 3.03e+03    -  7.63e-01 5.58e-05h 13
    4135  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.07e+04    -  7.63e-01 5.58e-05h 14
    4136  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.63e-01 5.58e-05h 14
    4137  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.63e-01 5.58e-05h 14
    4138  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.06e+04    -  7.63e-01 5.59e-05h 14
    4139  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.07e+04    -  7.63e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4140  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.05e+04    -  7.64e-01 5.59e-05h 14
    4141  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.07e+04    -  7.64e-01 5.58e-05h 14
    4142  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.08e+04    -  7.64e-01 5.58e-05h 14
    4143  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.07e+04    -  7.64e-01 5.58e-05h 14
    4144  0.0000000e+00 3.61e+05 1.27e+13  -1.0 5.08e+04    -  7.64e-01 4.58e-01w  1
    4145  0.0000000e+00 3.19e+05 6.14e+12  -1.0 4.90e+04    -  7.71e-01 3.92e-01w  1
    4146  0.0000000e+00 1.64e+05 4.81e+12  -1.0 5.66e+03    -  1.00e+00 4.95e-01w  1
    4147  0.0000000e+00 2.41e+04 6.20e+12  -1.0 3.05e+03    -  7.64e-01 5.58e-05h 13
    4148  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.07e+04    -  7.64e-01 5.58e-05h 14
    4149  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.07e+04    -  7.64e-01 5.58e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4150  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.05e+04    -  7.64e-01 5.59e-05h 14
    4151  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.07e+04    -  7.64e-01 5.58e-05h 14
    4152  0.0000000e+00 2.41e+04 6.20e+12  -1.0 5.07e+04    -  7.64e-01 5.58e-05h 14
    4153  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.05e+04    -  7.64e-01 5.59e-05h 14
    4154  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.07e+04    -  7.64e-01 5.59e-05h 14
    4155  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.05e+04    -  7.64e-01 5.59e-05h 14
    4156  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.64e-01 5.59e-05h 14
    4157  0.0000000e+00 3.61e+05 1.27e+13  -1.0 5.07e+04    -  7.64e-01 4.58e-01w  1
    4158  0.0000000e+00 3.19e+05 6.14e+12  -1.0 4.89e+04    -  7.72e-01 3.92e-01w  1
    4159  0.0000000e+00 1.63e+05 4.81e+12  -1.0 5.66e+03    -  1.00e+00 4.95e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4160  0.0000000e+00 2.40e+04 6.20e+12  -1.0 3.06e+03    -  7.64e-01 5.59e-05h 13
    4161  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.07e+04    -  7.65e-01 5.59e-05h 14
    4162  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.07e+04    -  7.65e-01 5.59e-05h 14
    4163  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.07e+04    -  7.65e-01 5.59e-05h 14
    4164  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.07e+04    -  7.65e-01 5.59e-05h 14
    4165  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.65e-01 5.59e-05h 14
    4166  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.07e+04    -  7.65e-01 5.59e-05h 14
    4167  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.07e+04    -  7.65e-01 5.59e-05h 14
    4168  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.07e+04    -  7.65e-01 5.59e-05h 14
    4169  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.65e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4170  0.0000000e+00 3.60e+05 1.27e+13  -1.0 5.07e+04    -  7.65e-01 4.58e-01w  1
    4171  0.0000000e+00 3.18e+05 6.14e+12  -1.0 4.89e+04    -  7.72e-01 3.93e-01w  1
    4172  0.0000000e+00 1.63e+05 4.81e+12  -1.0 5.69e+03    -  1.00e+00 4.95e-01w  1
    4173  0.0000000e+00 2.40e+04 6.20e+12  -1.0 3.07e+03    -  7.65e-01 5.59e-05h 13
    4174  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.65e-01 5.59e-05h 14
    4175  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.65e-01 5.59e-05h 14
    4176  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.07e+04    -  7.65e-01 5.59e-05h 14
    4177  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.05e+04    -  7.65e-01 5.59e-05h 14
    4178  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.05e+04    -  7.66e-01 5.59e-05h 14
    4179  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.66e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4180  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.05e+04    -  7.66e-01 5.59e-05h 14
    4181  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.05e+04    -  7.66e-01 5.59e-05h 14
    4182  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.66e-01 5.59e-05h 14
    4183  0.0000000e+00 3.60e+05 1.27e+13  -1.0 5.06e+04    -  7.66e-01 4.58e-01w  1
    4184  0.0000000e+00 3.18e+05 6.14e+12  -1.0 4.88e+04    -  7.73e-01 3.93e-01w  1
    4185  0.0000000e+00 1.63e+05 4.82e+12  -1.0 5.69e+03    -  1.00e+00 4.95e-01w  1
    4186  0.0000000e+00 2.40e+04 6.20e+12  -1.0 3.07e+03    -  7.66e-01 5.59e-05h 13
    4187  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.66e-01 5.59e-05h 14
    4188  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.66e-01 5.59e-05h 14
    4189  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.07e+04    -  7.66e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4190  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.05e+04    -  7.66e-01 5.59e-05h 14
    4191  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.66e-01 5.59e-05h 14
    4192  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.07e+04    -  7.66e-01 5.59e-05h 14
    4193  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.05e+04    -  7.66e-01 5.59e-05h 14
    4194  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.66e-01 5.59e-05h 14
    4195  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.66e-01 5.59e-05h 14
    4196  0.0000000e+00 3.59e+05 1.27e+13  -1.0 5.06e+04    -  7.66e-01 4.58e-01w  1
    4197  0.0000000e+00 3.17e+05 6.14e+12  -1.0 4.88e+04    -  7.73e-01 3.93e-01w  1
    4198  0.0000000e+00 1.62e+05 4.82e+12  -1.0 5.74e+03    -  1.00e+00 4.95e-01w  1
    4199  0.0000000e+00 2.40e+04 6.20e+12  -1.0 3.09e+03    -  7.66e-01 5.59e-05h 13
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4200  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.67e-01 5.59e-05h 14
    4201  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.05e+04    -  7.67e-01 5.59e-05h 14
    4202  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.67e-01 5.59e-05h 14
    4203  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.67e-01 5.59e-05h 14
    4204  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.05e+04    -  7.67e-01 5.59e-05h 14
    4205  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.67e-01 5.59e-05h 14
    4206  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.05e+04    -  7.67e-01 5.59e-05h 14
    4207  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.67e-01 5.59e-05h 14
    4208  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.67e-01 5.59e-05h 14
    4209  0.0000000e+00 3.59e+05 1.27e+13  -1.0 5.06e+04    -  7.67e-01 4.58e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4210  0.0000000e+00 3.17e+05 6.15e+12  -1.0 4.87e+04    -  7.74e-01 3.93e-01w  1
    4211  0.0000000e+00 1.62e+05 4.83e+12  -1.0 5.77e+03    -  1.00e+00 4.95e-01w  1
    4212  0.0000000e+00 2.40e+04 6.20e+12  -1.0 3.10e+03    -  7.67e-01 5.59e-05h 13
    4213  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.67e-01 5.59e-05h 14
    4214  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.67e-01 5.59e-05h 14
    4215  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.67e-01 5.59e-05h 14
    4216  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.67e-01 5.59e-05h 14
    4217  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.68e-01 5.59e-05h 14
    4218  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.68e-01 5.59e-05h 14
    4219  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.04e+04    -  7.68e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4220  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.68e-01 5.59e-05h 14
    4221  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.68e-01 5.59e-05h 14
    4222  0.0000000e+00 3.59e+05 1.27e+13  -1.0 5.06e+04    -  7.68e-01 4.58e-01w  1
    4223  0.0000000e+00 3.16e+05 6.14e+12  -1.0 4.86e+04    -  7.74e-01 3.93e-01w  1
    4224  0.0000000e+00 1.62e+05 4.83e+12  -1.0 5.76e+03    -  1.00e+00 4.95e-01w  1
    4225  0.0000000e+00 2.40e+04 6.20e+12  -1.0 3.11e+03    -  7.68e-01 5.59e-05h 13
    4226  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.04e+04    -  7.68e-01 5.59e-05h 14
    4227  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.05e+04    -  7.68e-01 5.59e-05h 14
    4228  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.68e-01 5.59e-05h 14
    4229  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.05e+04    -  7.68e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4230  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.04e+04    -  7.68e-01 5.59e-05h 14
    4231  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.06e+04    -  7.68e-01 5.59e-05h 14
    4232  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.05e+04    -  7.68e-01 5.59e-05h 14
    4233  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.05e+04    -  7.68e-01 5.59e-05h 14
    4234  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.04e+04    -  7.69e-01 5.59e-05h 14
    4235  0.0000000e+00 3.58e+05 1.27e+13  -1.0 5.05e+04    -  7.69e-01 4.58e-01w  1
    4236  0.0000000e+00 3.16e+05 6.14e+12  -1.0 4.86e+04    -  7.75e-01 3.94e-01w  1
    4237  0.0000000e+00 1.61e+05 4.83e+12  -1.0 5.79e+03    -  1.00e+00 4.95e-01w  1
    4238  0.0000000e+00 2.40e+04 6.20e+12  -1.0 3.11e+03    -  7.69e-01 5.59e-05h 13
    4239  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.05e+04    -  7.69e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4240  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.04e+04    -  7.69e-01 5.59e-05h 14
    4241  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.05e+04    -  7.69e-01 5.59e-05h 14
    4242  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.01e+04    -  7.69e-01 5.59e-05h 14
    4243  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.01e+04    -  7.69e-01 5.59e-05h 14
    4244  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.04e+04    -  7.69e-01 5.59e-05h 14
    4245  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.04e+04    -  7.69e-01 5.59e-05h 14
    4246  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.05e+04    -  7.69e-01 5.59e-05h 14
    4247  0.0000000e+00 2.40e+04 6.20e+12  -1.0 5.05e+04    -  7.69e-01 5.59e-05h 14
    4248  0.0000000e+00 3.58e+05 1.27e+13  -1.0 5.05e+04    -  7.69e-01 4.58e-01w  1
    4249  0.0000000e+00 3.15e+05 6.14e+12  -1.0 4.85e+04    -  7.76e-01 3.94e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4250  0.0000000e+00 1.61e+05 4.84e+12  -1.0 5.81e+03    -  1.00e+00 4.95e-01w  1
    4251  0.0000000e+00 2.39e+04 6.20e+12  -1.0 3.13e+03    -  7.69e-01 5.59e-05h 13
    4252  0.0000000e+00 2.39e+04 6.20e+12  -1.0 5.05e+04    -  7.69e-01 5.59e-05h 14
    4253  0.0000000e+00 2.39e+04 6.20e+12  -1.0 5.05e+04    -  7.69e-01 5.59e-05h 14
    4254  0.0000000e+00 2.39e+04 6.20e+12  -1.0 5.05e+04    -  7.69e-01 5.59e-05h 14
    4255  0.0000000e+00 2.39e+04 6.20e+12  -1.0 5.05e+04    -  7.70e-01 5.59e-05h 14
    4256  0.0000000e+00 2.39e+04 6.20e+12  -1.0 5.05e+04    -  7.70e-01 5.59e-05h 14
    4257  0.0000000e+00 2.39e+04 6.20e+12  -1.0 5.05e+04    -  7.70e-01 5.59e-05h 14
    4258  0.0000000e+00 2.39e+04 6.20e+12  -1.0 5.05e+04    -  7.70e-01 5.59e-05h 14
    4259  0.0000000e+00 2.39e+04 6.20e+12  -1.0 5.04e+04    -  7.70e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4260  0.0000000e+00 2.39e+04 6.20e+12  -1.0 5.05e+04    -  7.70e-01 5.59e-05h 14
    4261  0.0000000e+00 3.57e+05 1.27e+13  -1.0 5.05e+04    -  7.70e-01 4.58e-01w  1
    4262  0.0000000e+00 3.15e+05 6.15e+12  -1.0 4.85e+04    -  7.76e-01 3.94e-01w  1
    4263  0.0000000e+00 1.61e+05 4.84e+12  -1.0 5.83e+03    -  1.00e+00 4.95e-01w  1
    4264  0.0000000e+00 2.39e+04 6.20e+12  -1.0 3.14e+03    -  7.70e-01 5.59e-05h 13
    4265  0.0000000e+00 2.39e+04 6.20e+12  -1.0 5.04e+04    -  7.70e-01 5.59e-05h 14
    4266  0.0000000e+00 2.39e+04 6.20e+12  -1.0 5.04e+04    -  7.70e-01 5.59e-05h 14
    4267  0.0000000e+00 2.39e+04 6.20e+12  -1.0 5.05e+04    -  7.70e-01 5.59e-05h 14
    4268  0.0000000e+00 2.39e+04 6.20e+12  -1.0 5.05e+04    -  7.70e-01 5.59e-05h 14
    4269  0.0000000e+00 2.39e+04 6.20e+12  -1.0 5.05e+04    -  7.70e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4270  0.0000000e+00 2.39e+04 6.20e+12  -1.0 5.04e+04    -  7.70e-01 5.59e-05h 14
    4271  0.0000000e+00 2.39e+04 6.20e+12  -1.0 5.04e+04    -  7.70e-01 5.59e-05h 14
    4272  0.0000000e+00 2.39e+04 6.20e+12  -1.0 5.04e+04    -  7.71e-01 5.59e-05h 14
    4273  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.03e+04    -  7.71e-01 5.59e-05h 14
    4274  0.0000000e+00 3.57e+05 1.27e+13  -1.0 5.04e+04    -  7.71e-01 4.58e-01w  1
    4275  0.0000000e+00 3.14e+05 6.15e+12  -1.0 4.84e+04    -  7.76e-01 3.94e-01w  1
    4276  0.0000000e+00 1.61e+05 4.85e+12  -1.0 5.84e+03    -  1.00e+00 4.95e-01w  1
    4277  0.0000000e+00 2.39e+04 6.21e+12  -1.0 3.15e+03    -  7.71e-01 5.59e-05h 13
    4278  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.05e+04    -  7.71e-01 5.59e-05h 14
    4279  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.04e+04    -  7.71e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4280  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.05e+04    -  7.71e-01 5.59e-05h 14
    4281  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.05e+04    -  7.71e-01 5.59e-05h 14
    4282  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.05e+04    -  7.71e-01 5.59e-05h 14
    4283  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.03e+04    -  7.71e-01 5.59e-05h 14
    4284  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.03e+04    -  7.71e-01 5.59e-05h 14
    4285  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.05e+04    -  7.71e-01 5.59e-05h 14
    4286  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.04e+04    -  7.71e-01 5.59e-05h 14
    4287  0.0000000e+00 3.57e+05 1.27e+13  -1.0 5.03e+04    -  7.71e-01 4.58e-01w  1
    4288  0.0000000e+00 3.14e+05 6.15e+12  -1.0 4.83e+04    -  7.77e-01 3.94e-01w  1
    4289  0.0000000e+00 1.60e+05 4.85e+12  -1.0 5.83e+03    -  1.00e+00 4.95e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4290  0.0000000e+00 2.39e+04 6.21e+12  -1.0 3.14e+03    -  7.71e-01 5.59e-05h 13
    4291  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.04e+04    -  7.71e-01 5.59e-05h 14
    4292  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.04e+04    -  7.71e-01 5.59e-05h 14
    4293  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.03e+04    -  7.71e-01 5.59e-05h 14
    4294  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.04e+04    -  7.72e-01 5.59e-05h 14
    4295  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.03e+04    -  7.72e-01 5.59e-05h 14
    4296  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.01e+04    -  7.72e-01 5.59e-05h 14
    4297  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.04e+04    -  7.72e-01 5.59e-05h 14
    4298  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.04e+04    -  7.72e-01 5.59e-05h 14
    4299  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.04e+04    -  7.72e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4300  0.0000000e+00 3.56e+05 1.27e+13  -1.0 5.04e+04    -  7.72e-01 4.58e-01w  1
    4301  0.0000000e+00 3.13e+05 6.15e+12  -1.0 4.83e+04    -  7.77e-01 3.95e-01w  1
    4302  0.0000000e+00 1.60e+05 4.86e+12  -1.0 5.89e+03    -  1.00e+00 4.95e-01w  1
    4303  0.0000000e+00 2.39e+04 6.21e+12  -1.0 3.17e+03    -  7.72e-01 5.59e-05h 13
    4304  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.04e+04    -  7.72e-01 5.59e-05h 14
    4305  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.04e+04    -  7.72e-01 5.59e-05h 14
    4306  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.03e+04    -  7.72e-01 5.59e-05h 14
    4307  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.03e+04    -  7.72e-01 5.59e-05h 14
    4308  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.03e+04    -  7.72e-01 5.59e-05h 14
    4309  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.03e+04    -  7.72e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4310  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.04e+04    -  7.72e-01 5.59e-05h 14
    4311  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.04e+04    -  7.73e-01 5.59e-05h 14
    4312  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.04e+04    -  7.73e-01 5.59e-05h 14
    4313  0.0000000e+00 3.56e+05 1.27e+13  -1.0 5.04e+04    -  7.73e-01 4.58e-01w  1
    4314  0.0000000e+00 3.13e+05 6.15e+12  -1.0 4.83e+04    -  7.78e-01 3.95e-01w  1
    4315  0.0000000e+00 1.60e+05 4.86e+12  -1.0 5.91e+03    -  1.00e+00 4.95e-01w  1
    4316  0.0000000e+00 2.39e+04 6.21e+12  -1.0 3.18e+03    -  7.73e-01 5.59e-05h 13
    4317  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.04e+04    -  7.73e-01 5.59e-05h 14
    4318  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.03e+04    -  7.73e-01 5.59e-05h 14
    4319  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.03e+04    -  7.73e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4320  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.04e+04    -  7.73e-01 5.59e-05h 14
    4321  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.04e+04    -  7.73e-01 5.59e-05h 14
    4322  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.04e+04    -  7.73e-01 5.59e-05h 14
    4323  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.04e+04    -  7.73e-01 5.59e-05h 14
    4324  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.03e+04    -  7.73e-01 5.59e-05h 14
    4325  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.04e+04    -  7.73e-01 5.59e-05h 14
    4326  0.0000000e+00 3.55e+05 1.27e+13  -1.0 5.02e+04    -  7.73e-01 4.58e-01w  1
    4327  0.0000000e+00 3.12e+05 6.15e+12  -1.0 4.82e+04    -  7.78e-01 3.95e-01w  1
    4328  0.0000000e+00 1.60e+05 4.86e+12  -1.0 5.88e+03    -  1.00e+00 4.95e-01w  1
    4329  0.0000000e+00 2.39e+04 6.21e+12  -1.0 3.17e+03    -  7.73e-01 5.59e-05h 13
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4330  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.03e+04    -  7.73e-01 5.59e-05h 14
    4331  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.04e+04    -  7.73e-01 5.59e-05h 14
    4332  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.02e+04    -  7.74e-01 5.59e-05h 14
    4333  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.04e+04    -  7.74e-01 5.59e-05h 14
    4334  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.04e+04    -  7.74e-01 5.59e-05h 14
    4335  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.04e+04    -  7.74e-01 5.59e-05h 14
    4336  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.03e+04    -  7.74e-01 5.59e-05h 14
    4337  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.02e+04    -  7.74e-01 5.59e-05h 14
    4338  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.03e+04    -  7.74e-01 5.59e-05h 14
    4339  0.0000000e+00 3.55e+05 1.27e+13  -1.0 5.03e+04    -  7.74e-01 4.58e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4340  0.0000000e+00 3.12e+05 6.15e+12  -1.0 4.82e+04    -  7.79e-01 3.95e-01w  1
    4341  0.0000000e+00 1.59e+05 4.87e+12  -1.0 5.94e+03    -  1.00e+00 4.95e-01w  1
    4342  0.0000000e+00 2.39e+04 6.21e+12  -1.0 3.20e+03    -  7.74e-01 5.59e-05h 13
    4343  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.03e+04    -  7.74e-01 5.59e-05h 14
    4344  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.02e+04    -  7.74e-01 5.59e-05h 14
    4345  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.03e+04    -  7.74e-01 5.59e-05h 14
    4346  0.0000000e+00 2.39e+04 6.21e+12  -1.0 5.02e+04    -  7.74e-01 5.59e-05h 14
    4347  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.03e+04    -  7.74e-01 5.59e-05h 14
    4348  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.03e+04    -  7.74e-01 5.59e-05h 14
    4349  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.03e+04    -  7.74e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4350  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.03e+04    -  7.75e-01 5.59e-05h 14
    4351  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.03e+04    -  7.75e-01 5.59e-05h 14
    4352  0.0000000e+00 3.54e+05 1.27e+13  -1.0 5.02e+04    -  7.75e-01 4.58e-01w  1
    4353  0.0000000e+00 3.11e+05 6.15e+12  -1.0 4.81e+04    -  7.79e-01 3.95e-01w  1
    4354  0.0000000e+00 1.59e+05 4.87e+12  -1.0 5.91e+03    -  1.00e+00 4.95e-01w  1
    4355  0.0000000e+00 2.38e+04 6.21e+12  -1.0 3.19e+03    -  7.75e-01 5.59e-05h 13
    4356  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.01e+04    -  7.75e-01 5.59e-05h 14
    4357  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.75e-01 5.59e-05h 14
    4358  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.03e+04    -  7.75e-01 5.59e-05h 14
    4359  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.03e+04    -  7.75e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4360  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.03e+04    -  7.75e-01 5.59e-05h 14
    4361  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.03e+04    -  7.75e-01 5.59e-05h 14
    4362  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.03e+04    -  7.75e-01 5.59e-05h 14
    4363  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.75e-01 5.59e-05h 14
    4364  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.03e+04    -  7.75e-01 5.59e-05h 14
    4365  0.0000000e+00 3.54e+05 1.27e+13  -1.0 5.03e+04    -  7.75e-01 4.58e-01w  1
    4366  0.0000000e+00 3.11e+05 6.15e+12  -1.0 4.80e+04    -  7.80e-01 3.96e-01w  1
    4367  0.0000000e+00 1.59e+05 4.88e+12  -1.0 5.98e+03    -  1.00e+00 4.95e-01w  1
    4368  0.0000000e+00 2.38e+04 6.21e+12  -1.0 3.22e+03    -  7.75e-01 5.59e-05h 13
    4369  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.75e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4370  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.03e+04    -  7.76e-01 5.59e-05h 14
    4371  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.76e-01 5.59e-05h 14
    4372  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.03e+04    -  7.76e-01 5.59e-05h 14
    4373  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.01e+04    -  7.76e-01 5.59e-05h 14
    4374  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.03e+04    -  7.76e-01 5.59e-05h 14
    4375  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.76e-01 5.59e-05h 14
    4376  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.00e+04    -  7.76e-01 5.59e-05h 14
    4377  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.03e+04    -  7.76e-01 5.59e-05h 14
    4378  0.0000000e+00 3.54e+05 1.27e+13  -1.0 5.01e+04    -  7.76e-01 4.58e-01w  1
    4379  0.0000000e+00 3.10e+05 6.15e+12  -1.0 4.80e+04    -  7.80e-01 3.96e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4380  0.0000000e+00 1.58e+05 4.88e+12  -1.0 5.96e+03    -  1.00e+00 4.95e-01w  1
    4381  0.0000000e+00 2.38e+04 6.21e+12  -1.0 3.21e+03    -  7.76e-01 5.59e-05h 13
    4382  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.76e-01 5.59e-05h 14
    4383  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.76e-01 5.59e-05h 14
    4384  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.76e-01 5.59e-05h 14
    4385  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.76e-01 5.59e-05h 14
    4386  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.76e-01 5.59e-05h 14
    4387  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.76e-01 5.59e-05h 14
    4388  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.77e-01 5.59e-05h 14
    4389  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.77e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4390  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.77e-01 5.59e-05h 14
    4391  0.0000000e+00 3.53e+05 1.27e+13  -1.0 5.02e+04    -  7.77e-01 4.58e-01w  1
    4392  0.0000000e+00 3.10e+05 6.15e+12  -1.0 4.80e+04    -  7.81e-01 3.96e-01w  1
    4393  0.0000000e+00 1.58e+05 4.88e+12  -1.0 6.02e+03    -  1.00e+00 4.95e-01w  1
    4394  0.0000000e+00 2.38e+04 6.21e+12  -1.0 3.24e+03    -  7.77e-01 5.59e-05h 13
    4395  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.77e-01 5.59e-05h 14
    4396  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.77e-01 5.59e-05h 14
    4397  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.77e-01 5.59e-05h 14
    4398  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.77e-01 5.59e-05h 14
    4399  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.77e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4400  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.77e-01 5.59e-05h 14
    4401  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.77e-01 5.59e-05h 14
    4402  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.77e-01 5.59e-05h 14
    4403  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.01e+04    -  7.77e-01 5.59e-05h 14
    4404  0.0000000e+00 3.53e+05 1.27e+13  -1.0 5.02e+04    -  7.77e-01 4.58e-01w  1
    4405  0.0000000e+00 3.09e+05 6.15e+12  -1.0 4.79e+04    -  7.81e-01 3.96e-01w  1
    4406  0.0000000e+00 1.58e+05 4.89e+12  -1.0 6.04e+03    -  1.00e+00 4.95e-01w  1
    4407  0.0000000e+00 2.38e+04 6.21e+12  -1.0 3.25e+03    -  7.77e-01 5.59e-05h 13
    4408  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.01e+04    -  7.77e-01 5.59e-05h 14
    4409  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.78e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4410  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.78e-01 5.59e-05h 14
    4411  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.00e+04    -  7.78e-01 5.59e-05h 14
    4412  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.78e-01 5.59e-05h 14
    4413  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.78e-01 5.59e-05h 14
    4414  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.78e-01 5.59e-05h 14
    4415  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.01e+04    -  7.78e-01 5.59e-05h 14
    4416  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.78e-01 5.59e-05h 14
    4417  0.0000000e+00 3.52e+05 1.27e+13  -1.0 5.02e+04    -  7.78e-01 4.58e-01w  1
    4418  0.0000000e+00 3.09e+05 6.15e+12  -1.0 4.78e+04    -  7.82e-01 3.96e-01w  1
    4419  0.0000000e+00 1.58e+05 4.89e+12  -1.0 6.05e+03    -  1.00e+00 4.95e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4420  0.0000000e+00 2.38e+04 6.21e+12  -1.0 3.25e+03    -  7.78e-01 5.59e-05h 13
    4421  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.78e-01 5.59e-05h 14
    4422  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.78e-01 5.59e-05h 14
    4423  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.78e-01 5.59e-05h 14
    4424  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.78e-01 5.59e-05h 14
    4425  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.78e-01 5.59e-05h 14
    4426  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.01e+04    -  7.79e-01 5.59e-05h 14
    4427  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.01e+04    -  7.79e-01 5.59e-05h 14
    4428  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.01e+04    -  7.79e-01 5.59e-05h 14
    4429  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.79e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4430  0.0000000e+00 3.52e+05 1.27e+13  -1.0 5.02e+04    -  7.79e-01 4.58e-01w  1
    4431  0.0000000e+00 3.08e+05 6.15e+12  -1.0 4.78e+04    -  7.82e-01 3.96e-01w  1
    4432  0.0000000e+00 1.57e+05 4.90e+12  -1.0 6.08e+03    -  1.00e+00 4.95e-01w  1
    4433  0.0000000e+00 2.38e+04 6.21e+12  -1.0 3.27e+03    -  7.79e-01 5.59e-05h 13
    4434  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.01e+04    -  7.79e-01 5.59e-05h 14
    4435  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.79e-01 5.59e-05h 14
    4436  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.02e+04    -  7.79e-01 5.59e-05h 14
    4437  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.01e+04    -  7.79e-01 5.59e-05h 14
    4438  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.01e+04    -  7.79e-01 5.59e-05h 14
    4439  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.01e+04    -  7.79e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4440  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.01e+04    -  7.79e-01 5.59e-05h 14
    4441  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.01e+04    -  7.79e-01 5.59e-05h 14
    4442  0.0000000e+00 2.38e+04 6.21e+12  -1.0 5.01e+04    -  7.79e-01 5.59e-05h 14
    4443  0.0000000e+00 3.52e+05 1.27e+13  -1.0 5.01e+04    -  7.79e-01 4.58e-01w  1
    4444  0.0000000e+00 3.08e+05 6.15e+12  -1.0 4.77e+04    -  7.83e-01 3.97e-01w  1
    4445  0.0000000e+00 1.57e+05 4.90e+12  -1.0 6.10e+03    -  1.00e+00 4.95e-01w  1
    4446  0.0000000e+00 2.37e+04 6.21e+12  -1.0 3.28e+03    -  7.79e-01 5.59e-05h 13
    4447  0.0000000e+00 2.37e+04 6.21e+12  -1.0 5.01e+04    -  7.79e-01 5.59e-05h 14
    4448  0.0000000e+00 2.37e+04 6.21e+12  -1.0 5.01e+04    -  7.80e-01 5.59e-05h 14
    4449  0.0000000e+00 2.37e+04 6.21e+12  -1.0 5.01e+04    -  7.80e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4450  0.0000000e+00 2.37e+04 6.21e+12  -1.0 5.01e+04    -  7.80e-01 5.59e-05h 14
    4451  0.0000000e+00 2.37e+04 6.21e+12  -1.0 5.01e+04    -  7.80e-01 5.59e-05h 14
    4452  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.01e+04    -  7.80e-01 5.59e-05h 14
    4453  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.01e+04    -  7.80e-01 5.59e-05h 14
    4454  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.80e-01 5.59e-05h 14
    4455  0.0000000e+00 2.37e+04 6.22e+12  -1.0 4.99e+04    -  7.80e-01 5.59e-05h 14
    4456  0.0000000e+00 3.51e+05 1.27e+13  -1.0 5.01e+04    -  7.80e-01 4.58e-01w  1
    4457  0.0000000e+00 3.07e+05 6.15e+12  -1.0 4.77e+04    -  7.83e-01 3.97e-01w  1
    4458  0.0000000e+00 1.57e+05 4.91e+12  -1.0 6.12e+03    -  1.00e+00 4.95e-01w  1
    4459  0.0000000e+00 2.37e+04 6.22e+12  -1.0 3.29e+03    -  7.80e-01 5.59e-05h 13
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4460  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.80e-01 5.59e-05h 14
    4461  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.01e+04    -  7.80e-01 5.59e-05h 14
    4462  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.01e+04    -  7.80e-01 5.59e-05h 14
    4463  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.80e-01 5.59e-05h 14
    4464  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.81e-01 5.59e-05h 14
    4465  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.01e+04    -  7.81e-01 5.59e-05h 14
    4466  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.01e+04    -  7.81e-01 5.59e-05h 14
    4467  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.81e-01 5.59e-05h 14
    4468  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.01e+04    -  7.81e-01 5.59e-05h 14
    4469  0.0000000e+00 3.51e+05 1.27e+13  -1.0 4.97e+04    -  7.81e-01 4.58e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4470  0.0000000e+00 3.07e+05 6.15e+12  -1.0 4.75e+04    -  7.84e-01 3.97e-01w  1
    4471  0.0000000e+00 1.57e+05 4.91e+12  -1.0 6.03e+03    -  1.00e+00 4.95e-01w  1
    4472  0.0000000e+00 2.37e+04 6.22e+12  -1.0 3.25e+03    -  7.81e-01 5.59e-05h 13
    4473  0.0000000e+00 2.37e+04 6.22e+12  -1.0 4.97e+04    -  7.81e-01 5.59e-05h 14
    4474  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.81e-01 5.59e-05h 14
    4475  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.01e+04    -  7.81e-01 5.59e-05h 14
    4476  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.01e+04    -  7.81e-01 5.59e-05h 14
    4477  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.01e+04    -  7.81e-01 5.59e-05h 14
    4478  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.81e-01 5.59e-05h 14
    4479  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.81e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4480  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.01e+04    -  7.81e-01 5.59e-05h 14
    4481  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.81e-01 5.59e-05h 14
    4482  0.0000000e+00 3.50e+05 1.27e+13  -1.0 5.00e+04    -  7.81e-01 4.58e-01w  1
    4483  0.0000000e+00 3.06e+05 6.16e+12  -1.0 4.76e+04    -  7.84e-01 3.97e-01w  1
    4484  0.0000000e+00 1.56e+05 4.91e+12  -1.0 6.15e+03    -  1.00e+00 4.95e-01w  1
    4485  0.0000000e+00 2.37e+04 6.22e+12  -1.0 3.30e+03    -  7.81e-01 5.59e-05h 13
    4486  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.82e-01 5.59e-05h 14
    4487  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.82e-01 5.59e-05h 14
    4488  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.82e-01 5.59e-05h 14
    4489  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.82e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4490  0.0000000e+00 2.37e+04 6.22e+12  -1.0 4.99e+04    -  7.82e-01 5.59e-05h 14
    4491  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.82e-01 5.59e-05h 14
    4492  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.82e-01 5.59e-05h 14
    4493  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.82e-01 5.59e-05h 14
    4494  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.82e-01 5.59e-05h 14
    4495  0.0000000e+00 3.50e+05 1.27e+13  -1.0 5.00e+04    -  7.82e-01 4.58e-01w  1
    4496  0.0000000e+00 3.06e+05 6.16e+12  -1.0 4.75e+04    -  7.85e-01 3.97e-01w  1
    4497  0.0000000e+00 1.56e+05 4.92e+12  -1.0 6.17e+03    -  1.00e+00 4.95e-01w  1
    4498  0.0000000e+00 2.37e+04 6.22e+12  -1.0 3.32e+03    -  7.82e-01 5.59e-05h 13
    4499  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.82e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4500  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.82e-01 5.59e-05h 14
    4501  0.0000000e+00 2.37e+04 6.22e+12  -1.0 4.98e+04    -  7.82e-01 5.59e-05h 14
    4502  0.0000000e+00 2.37e+04 6.22e+12  -1.0 4.98e+04    -  7.83e-01 5.59e-05h 14
    4503  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.83e-01 5.59e-05h 14
    4504  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.83e-01 5.59e-05h 14
    4505  0.0000000e+00 2.37e+04 6.22e+12  -1.0 4.97e+04    -  7.83e-01 5.59e-05h 14
    4506  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.83e-01 5.59e-05h 14
    4507  0.0000000e+00 2.37e+04 6.22e+12  -1.0 4.99e+04    -  7.83e-01 5.59e-05h 14
    4508  0.0000000e+00 3.49e+05 1.27e+13  -1.0 5.00e+04    -  7.83e-01 4.58e-01w  1
    4509  0.0000000e+00 3.05e+05 6.16e+12  -1.0 4.75e+04    -  7.85e-01 3.98e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4510  0.0000000e+00 1.56e+05 4.92e+12  -1.0 6.18e+03    -  1.00e+00 4.95e-01w  1
    4511  0.0000000e+00 2.37e+04 6.22e+12  -1.0 3.32e+03    -  7.83e-01 5.59e-05h 13
    4512  0.0000000e+00 2.37e+04 6.22e+12  -1.0 4.99e+04    -  7.83e-01 5.59e-05h 14
    4513  0.0000000e+00 2.37e+04 6.22e+12  -1.0 4.98e+04    -  7.83e-01 5.59e-05h 14
    4514  0.0000000e+00 2.37e+04 6.22e+12  -1.0 4.98e+04    -  7.83e-01 5.59e-05h 14
    4515  0.0000000e+00 2.37e+04 6.22e+12  -1.0 4.98e+04    -  7.83e-01 5.59e-05h 14
    4516  0.0000000e+00 2.37e+04 6.22e+12  -1.0 4.99e+04    -  7.83e-01 5.59e-05h 14
    4517  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.83e-01 5.59e-05h 14
    4518  0.0000000e+00 2.37e+04 6.22e+12  -1.0 4.98e+04    -  7.83e-01 5.59e-05h 14
    4519  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.83e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4520  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.83e-01 5.59e-05h 14
    4521  0.0000000e+00 3.49e+05 1.27e+13  -1.0 5.00e+04    -  7.84e-01 4.58e-01w  1
    4522  0.0000000e+00 3.05e+05 6.16e+12  -1.0 4.74e+04    -  7.86e-01 3.98e-01w  1
    4523  0.0000000e+00 1.56e+05 4.93e+12  -1.0 6.20e+03    -  1.00e+00 4.95e-01w  1
    4524  0.0000000e+00 2.37e+04 6.22e+12  -1.0 3.33e+03    -  7.84e-01 5.59e-05h 13
    4525  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.84e-01 5.59e-05h 14
    4526  0.0000000e+00 2.37e+04 6.22e+12  -1.0 4.99e+04    -  7.84e-01 5.59e-05h 14
    4527  0.0000000e+00 2.37e+04 6.22e+12  -1.0 4.99e+04    -  7.84e-01 5.59e-05h 14
    4528  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.84e-01 5.59e-05h 14
    4529  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.84e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4530  0.0000000e+00 2.37e+04 6.22e+12  -1.0 4.99e+04    -  7.84e-01 5.59e-05h 14
    4531  0.0000000e+00 2.37e+04 6.22e+12  -1.0 4.99e+04    -  7.84e-01 5.59e-05h 14
    4532  0.0000000e+00 2.37e+04 6.22e+12  -1.0 5.00e+04    -  7.84e-01 5.59e-05h 14
    4533  0.0000000e+00 2.37e+04 6.22e+12  -1.0 4.99e+04    -  7.84e-01 5.59e-05h 14
    4534  0.0000000e+00 3.49e+05 1.27e+13  -1.0 4.99e+04    -  7.84e-01 4.58e-01w  1
    4535  0.0000000e+00 3.04e+05 6.16e+12  -1.0 4.74e+04    -  7.86e-01 3.98e-01w  1
    4536  0.0000000e+00 1.55e+05 4.93e+12  -1.0 6.21e+03    -  1.00e+00 4.95e-01w  1
    4537  0.0000000e+00 2.37e+04 6.22e+12  -1.0 3.34e+03    -  7.84e-01 5.59e-05h 13
    4538  0.0000000e+00 2.37e+04 6.22e+12  -1.0 4.99e+04    -  7.84e-01 5.59e-05h 14
    4539  0.0000000e+00 2.37e+04 6.22e+12  -1.0 4.99e+04    -  7.84e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4540  0.0000000e+00 2.37e+04 6.22e+12  -1.0 4.99e+04    -  7.84e-01 5.59e-05h 14
    4541  0.0000000e+00 2.37e+04 6.22e+12  -1.0 4.99e+04    -  7.84e-01 5.59e-05h 14
    4542  0.0000000e+00 2.37e+04 6.22e+12  -1.0 4.99e+04    -  7.85e-01 5.59e-05h 14
    4543  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.99e+04    -  7.85e-01 5.59e-05h 14
    4544  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.99e+04    -  7.85e-01 5.59e-05h 14
    4545  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.99e+04    -  7.85e-01 5.59e-05h 14
    4546  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.99e+04    -  7.85e-01 5.59e-05h 14
    4547  0.0000000e+00 3.48e+05 1.27e+13  -1.0 4.99e+04    -  7.85e-01 4.58e-01w  1
    4548  0.0000000e+00 3.04e+05 6.16e+12  -1.0 4.73e+04    -  7.87e-01 3.98e-01w  1
    4549  0.0000000e+00 1.55e+05 4.93e+12  -1.0 6.24e+03    -  1.00e+00 4.95e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4550  0.0000000e+00 2.36e+04 6.22e+12  -1.0 3.35e+03    -  7.85e-01 5.59e-05h 13
    4551  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.99e+04    -  7.85e-01 5.59e-05h 14
    4552  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.99e+04    -  7.85e-01 5.59e-05h 14
    4553  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.97e+04    -  7.85e-01 5.59e-05h 14
    4554  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.99e+04    -  7.85e-01 5.59e-05h 14
    4555  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.99e+04    -  7.85e-01 5.59e-05h 14
    4556  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.85e-01 5.59e-05h 14
    4557  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.99e+04    -  7.85e-01 5.59e-05h 14
    4558  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.99e+04    -  7.85e-01 5.59e-05h 14
    4559  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.99e+04    -  7.85e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4560  0.0000000e+00 3.48e+05 1.27e+13  -1.0 4.98e+04    -  7.86e-01 4.58e-01w  1
    4561  0.0000000e+00 3.03e+05 6.17e+12  -1.0 4.73e+04    -  7.87e-01 3.98e-01w  1
    4562  0.0000000e+00 1.55e+05 4.94e+12  -1.0 6.24e+03    -  1.00e+00 4.95e-01w  1
    4563  0.0000000e+00 2.36e+04 6.22e+12  -1.0 3.35e+03    -  7.86e-01 5.59e-05h 13
    4564  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.99e+04    -  7.86e-01 5.59e-05h 14
    4565  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.99e+04    -  7.86e-01 5.59e-05h 14
    4566  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.86e-01 5.59e-05h 14
    4567  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.86e-01 5.59e-05h 14
    4568  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.86e-01 5.59e-05h 14
    4569  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.99e+04    -  7.86e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4570  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.97e+04    -  7.86e-01 5.59e-05h 14
    4571  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.99e+04    -  7.86e-01 5.59e-05h 14
    4572  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.86e-01 5.59e-05h 14
    4573  0.0000000e+00 3.47e+05 1.27e+13  -1.0 4.99e+04    -  7.86e-01 4.58e-01w  1
    4574  0.0000000e+00 3.03e+05 6.18e+12  -1.0 4.72e+04    -  7.88e-01 3.98e-01w  1
    4575  0.0000000e+00 1.54e+05 4.94e+12  -1.0 6.28e+03    -  1.00e+00 4.95e-01w  1
    4576  0.0000000e+00 2.36e+04 6.22e+12  -1.0 3.37e+03    -  7.86e-01 5.59e-05h 13
    4577  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.86e-01 5.59e-05h 14
    4578  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.86e-01 5.59e-05h 14
    4579  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.99e+04    -  7.86e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4580  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.99e+04    -  7.86e-01 5.59e-05h 14
    4581  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.99e+04    -  7.87e-01 5.59e-05h 14
    4582  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.87e-01 5.59e-05h 14
    4583  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.97e+04    -  7.87e-01 5.59e-05h 14
    4584  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.87e-01 5.59e-05h 14
    4585  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.87e-01 5.59e-05h 14
    4586  0.0000000e+00 3.47e+05 1.27e+13  -1.0 4.97e+04    -  7.87e-01 4.58e-01w  1
    4587  0.0000000e+00 3.02e+05 6.19e+12  -1.0 4.71e+04    -  7.88e-01 3.99e-01w  1
    4588  0.0000000e+00 1.54e+05 4.95e+12  -1.0 6.24e+03    -  1.00e+00 4.95e-01w  1
    4589  0.0000000e+00 2.36e+04 6.22e+12  -1.0 3.36e+03    -  7.87e-01 5.59e-05h 13
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4590  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.87e-01 5.59e-05h 14
    4591  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.87e-01 5.59e-05h 14
    4592  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.87e-01 5.59e-05h 14
    4593  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.97e+04    -  7.87e-01 5.59e-05h 14
    4594  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.87e-01 5.59e-05h 14
    4595  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.87e-01 5.59e-05h 14
    4596  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.87e-01 5.59e-05h 14
    4597  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.87e-01 5.59e-05h 14
    4598  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.87e-01 5.59e-05h 14
    4599  0.0000000e+00 3.47e+05 1.27e+13  -1.0 4.98e+04    -  7.88e-01 4.58e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4600  0.0000000e+00 3.02e+05 6.19e+12  -1.0 4.71e+04    -  7.89e-01 3.99e-01w  1
    4601  0.0000000e+00 1.54e+05 4.95e+12  -1.0 6.31e+03    -  1.00e+00 4.95e-01w  1
    4602  0.0000000e+00 2.36e+04 6.22e+12  -1.0 3.39e+03    -  7.88e-01 5.59e-05h 13
    4603  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.88e-01 5.59e-05h 14
    4604  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.88e-01 5.59e-05h 14
    4605  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.88e-01 5.59e-05h 14
    4606  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.88e-01 5.59e-05h 14
    4607  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.88e-01 5.59e-05h 14
    4608  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.88e-01 5.59e-05h 14
    4609  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.97e+04    -  7.88e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4610  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.88e-01 5.59e-05h 14
    4611  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.88e-01 5.59e-05h 14
    4612  0.0000000e+00 3.46e+05 1.27e+13  -1.0 4.96e+04    -  7.88e-01 4.58e-01w  1
    4613  0.0000000e+00 3.01e+05 6.20e+12  -1.0 4.70e+04    -  7.89e-01 3.99e-01w  1
    4614  0.0000000e+00 1.54e+05 4.96e+12  -1.0 6.27e+03    -  1.00e+00 4.95e-01w  1
    4615  0.0000000e+00 2.36e+04 6.22e+12  -1.0 3.37e+03    -  7.88e-01 5.59e-05h 13
    4616  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.95e+04    -  7.88e-01 5.59e-05h 14
    4617  0.0000000e+00 2.36e+04 6.22e+12  -1.0 4.98e+04    -  7.89e-01 5.59e-05h 14
    4618  0.0000000e+00 2.36e+04 6.23e+12  -1.0 4.97e+04    -  7.88e-01 5.59e-05h 14
    4619  0.0000000e+00 2.36e+04 6.23e+12  -1.0 4.98e+04    -  7.89e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4620  0.0000000e+00 2.36e+04 6.23e+12  -1.0 4.98e+04    -  7.89e-01 5.59e-05h 14
    4621  0.0000000e+00 2.36e+04 6.23e+12  -1.0 4.98e+04    -  7.89e-01 5.59e-05h 14
    4622  0.0000000e+00 2.36e+04 6.23e+12  -1.0 4.95e+04    -  7.89e-01 5.59e-05h 14
    4623  0.0000000e+00 2.36e+04 6.23e+12  -1.0 4.97e+04    -  7.89e-01 5.59e-05h 14
    4624  0.0000000e+00 2.36e+04 6.23e+12  -1.0 4.97e+04    -  7.89e-01 5.59e-05h 14
    4625  0.0000000e+00 3.46e+05 1.27e+13  -1.0 4.98e+04    -  7.89e-01 4.58e-01w  1
    4626  0.0000000e+00 3.01e+05 6.21e+12  -1.0 4.70e+04    -  7.90e-01 3.99e-01w  1
    4627  0.0000000e+00 1.53e+05 4.96e+12  -1.0 6.34e+03    -  1.00e+00 4.95e-01w  1
    4628  0.0000000e+00 2.36e+04 6.23e+12  -1.0 3.41e+03    -  7.89e-01 5.59e-05h 13
    4629  0.0000000e+00 2.36e+04 6.23e+12  -1.0 4.97e+04    -  7.89e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4630  0.0000000e+00 2.36e+04 6.23e+12  -1.0 4.97e+04    -  7.89e-01 5.59e-05h 14
    4631  0.0000000e+00 2.36e+04 6.23e+12  -1.0 4.97e+04    -  7.89e-01 5.59e-05h 14
    4632  0.0000000e+00 2.36e+04 6.23e+12  -1.0 4.97e+04    -  7.89e-01 5.59e-05h 14
    4633  0.0000000e+00 2.36e+04 6.23e+12  -1.0 4.97e+04    -  7.89e-01 5.59e-05h 14
    4634  0.0000000e+00 2.36e+04 6.23e+12  -1.0 4.97e+04    -  7.89e-01 5.59e-05h 14
    4635  0.0000000e+00 2.36e+04 6.23e+12  -1.0 4.97e+04    -  7.89e-01 5.59e-05h 14
    4636  0.0000000e+00 2.36e+04 6.23e+12  -1.0 4.97e+04    -  7.89e-01 5.59e-05h 14
    4637  0.0000000e+00 2.36e+04 6.23e+12  -1.0 4.97e+04    -  7.90e-01 5.59e-05h 14
    4638  0.0000000e+00 3.45e+05 1.27e+13  -1.0 4.97e+04    -  7.90e-01 4.58e-01w  1
    4639  0.0000000e+00 3.00e+05 6.22e+12  -1.0 4.70e+04    -  7.90e-01 3.99e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4640  0.0000000e+00 1.53e+05 4.96e+12  -1.0 6.36e+03    -  1.00e+00 4.95e-01w  1
    4641  0.0000000e+00 2.36e+04 6.23e+12  -1.0 3.42e+03    -  7.90e-01 5.59e-05h 13
    4642  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.97e+04    -  7.90e-01 5.59e-05h 14
    4643  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.97e+04    -  7.90e-01 5.59e-05h 14
    4644  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.97e+04    -  7.90e-01 5.59e-05h 14
    4645  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.97e+04    -  7.90e-01 5.59e-05h 14
    4646  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.97e+04    -  7.90e-01 5.59e-05h 14
    4647  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.97e+04    -  7.90e-01 5.59e-05h 14
    4648  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.90e-01 5.60e-05h 14
    4649  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.97e+04    -  7.90e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4650  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.97e+04    -  7.90e-01 5.59e-05h 14
    4651  0.0000000e+00 3.45e+05 1.27e+13  -1.0 4.97e+04    -  7.90e-01 4.58e-01w  1
    4652  0.0000000e+00 3.00e+05 6.23e+12  -1.0 4.69e+04    -  7.91e-01 3.99e-01w  1
    4653  0.0000000e+00 1.53e+05 4.97e+12  -1.0 6.38e+03    -  1.00e+00 4.95e-01w  1
    4654  0.0000000e+00 2.35e+04 6.23e+12  -1.0 3.42e+03    -  7.90e-01 5.59e-05h 13
    4655  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.97e+04    -  7.90e-01 5.59e-05h 14
    4656  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.97e+04    -  7.90e-01 5.59e-05h 14
    4657  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.97e+04    -  7.90e-01 5.59e-05h 14
    4658  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.97e+04    -  7.91e-01 5.59e-05h 14
    4659  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.97e+04    -  7.91e-01 5.59e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4660  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.97e+04    -  7.91e-01 5.59e-05h 14
    4661  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.91e-01 5.60e-05h 14
    4662  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.97e+04    -  7.91e-01 5.59e-05h 14
    4663  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.97e+04    -  7.91e-01 5.59e-05h 14
    4664  0.0000000e+00 3.45e+05 1.27e+13  -1.0 4.97e+04    -  7.91e-01 4.58e-01w  1
    4665  0.0000000e+00 2.99e+05 6.24e+12  -1.0 4.69e+04    -  7.91e-01 4.00e-01w  1
    4666  0.0000000e+00 1.53e+05 4.97e+12  -1.0 6.40e+03    -  1.00e+00 4.95e-01w  1
    4667  0.0000000e+00 2.35e+04 6.23e+12  -1.0 3.44e+03    -  7.91e-01 5.59e-05h 13
    4668  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.97e+04    -  7.91e-01 5.59e-05h 14
    4669  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.92e+04    -  7.91e-01 5.60e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4670  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.95e+04    -  7.91e-01 5.60e-05h 14
    4671  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.91e-01 5.60e-05h 14
    4672  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.97e+04    -  7.91e-01 5.60e-05h 14
    4673  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.97e+04    -  7.91e-01 5.60e-05h 14
    4674  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.93e+04    -  7.91e-01 5.60e-05h 14
    4675  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.92e-01 5.60e-05h 14
    4676  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.92e-01 5.60e-05h 14
    4677  0.0000000e+00 3.44e+05 1.27e+13  -1.0 4.96e+04    -  7.92e-01 4.58e-01w  1
    4678  0.0000000e+00 2.99e+05 6.25e+12  -1.0 4.68e+04    -  7.92e-01 4.00e-01w  1
    4679  0.0000000e+00 1.52e+05 4.98e+12  -1.0 6.41e+03    -  1.00e+00 4.95e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4680  0.0000000e+00 2.35e+04 6.23e+12  -1.0 3.44e+03    -  7.92e-01 5.60e-05h 13
    4681  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.97e+04    -  7.92e-01 5.60e-05h 14
    4682  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.92e-01 5.60e-05h 14
    4683  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.92e-01 5.60e-05h 14
    4684  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.92e-01 5.60e-05h 14
    4685  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.94e+04    -  7.92e-01 5.60e-05h 14
    4686  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.92e-01 5.60e-05h 14
    4687  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.92e-01 5.60e-05h 14
    4688  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.92e-01 5.60e-05h 14
    4689  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.95e+04    -  7.92e-01 5.60e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4690  0.0000000e+00 3.44e+05 1.27e+13  -1.0 4.96e+04    -  7.92e-01 4.58e-01w  1
    4691  0.0000000e+00 2.98e+05 6.26e+12  -1.0 4.68e+04    -  7.92e-01 4.00e-01w  1
    4692  0.0000000e+00 1.52e+05 4.98e+12  -1.0 6.41e+03    -  1.00e+00 4.95e-01w  1
    4693  0.0000000e+00 2.35e+04 6.23e+12  -1.0 3.44e+03    -  7.92e-01 5.60e-05h 13
    4694  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.95e+04    -  7.92e-01 5.60e-05h 14
    4695  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.95e+04    -  7.93e-01 5.60e-05h 14
    4696  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.92e-01 5.60e-05h 14
    4697  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.93e-01 5.60e-05h 14
    4698  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.93e-01 5.60e-05h 14
    4699  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.93e-01 5.60e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4700  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.95e+04    -  7.93e-01 5.60e-05h 14
    4701  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.93e-01 5.60e-05h 14
    4702  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.93e-01 5.60e-05h 14
    4703  0.0000000e+00 3.43e+05 1.27e+13  -1.0 4.96e+04    -  7.93e-01 4.58e-01w  1
    4704  0.0000000e+00 2.98e+05 6.27e+12  -1.0 4.67e+04    -  7.93e-01 4.00e-01w  1
    4705  0.0000000e+00 1.52e+05 4.99e+12  -1.0 6.43e+03    -  1.00e+00 4.95e-01w  1
    4706  0.0000000e+00 2.35e+04 6.23e+12  -1.0 3.45e+03    -  7.93e-01 5.60e-05h 13
    4707  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.93e-01 5.60e-05h 14
    4708  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.93e-01 5.60e-05h 14
    4709  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.93e-01 5.60e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4710  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.93e-01 5.60e-05h 14
    4711  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.95e+04    -  7.93e-01 5.60e-05h 14
    4712  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.93e-01 5.60e-05h 14
    4713  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.95e+04    -  7.93e-01 5.60e-05h 14
    4714  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.95e+04    -  7.94e-01 5.60e-05h 14
    4715  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.94e-01 5.60e-05h 14
    4716  0.0000000e+00 3.43e+05 1.27e+13  -1.0 4.95e+04    -  7.94e-01 4.58e-01w  1
    4717  0.0000000e+00 2.98e+05 6.28e+12  -1.0 4.67e+04    -  7.93e-01 4.00e-01w  1
    4718  0.0000000e+00 1.52e+05 4.99e+12  -1.0 6.46e+03    -  1.00e+00 4.95e-01w  1
    4719  0.0000000e+00 2.35e+04 6.23e+12  -1.0 3.47e+03    -  7.94e-01 5.60e-05h 13
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4720  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.94e-01 5.60e-05h 14
    4721  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.94e-01 5.60e-05h 14
    4722  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.93e+04    -  7.94e-01 5.60e-05h 14
    4723  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.95e+04    -  7.94e-01 5.60e-05h 14
    4724  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.95e+04    -  7.94e-01 5.60e-05h 14
    4725  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.94e-01 5.60e-05h 14
    4726  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.95e+04    -  7.94e-01 5.60e-05h 14
    4727  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.94e-01 5.60e-05h 14
    4728  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.94e-01 5.60e-05h 14
    4729  0.0000000e+00 3.43e+05 1.27e+13  -1.0 4.96e+04    -  7.94e-01 4.58e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4730  0.0000000e+00 2.97e+05 6.29e+12  -1.0 4.66e+04    -  7.94e-01 4.01e-01w  1
    4731  0.0000000e+00 1.51e+05 4.99e+12  -1.0 6.48e+03    -  1.00e+00 4.95e-01w  1
    4732  0.0000000e+00 2.35e+04 6.23e+12  -1.0 3.48e+03    -  7.94e-01 5.60e-05h 13
    4733  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.96e+04    -  7.94e-01 5.60e-05h 14
    4734  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.95e+04    -  7.94e-01 5.60e-05h 14
    4735  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.95e+04    -  7.94e-01 5.60e-05h 14
    4736  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.95e+04    -  7.95e-01 5.60e-05h 14
    4737  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.93e+04    -  7.95e-01 5.60e-05h 14
    4738  0.0000000e+00 2.35e+04 6.23e+12  -1.0 4.95e+04    -  7.95e-01 5.60e-05h 14
    4739  0.0000000e+00 2.34e+04 6.23e+12  -1.0 4.95e+04    -  7.95e-01 5.60e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4740  0.0000000e+00 2.34e+04 6.23e+12  -1.0 4.95e+04    -  7.95e-01 5.60e-05h 14
    4741  0.0000000e+00 2.34e+04 6.23e+12  -1.0 4.95e+04    -  7.95e-01 5.60e-05h 14
    4742  0.0000000e+00 3.42e+05 1.27e+13  -1.0 4.95e+04    -  7.95e-01 4.58e-01w  1
    4743  0.0000000e+00 2.97e+05 6.30e+12  -1.0 4.66e+04    -  7.94e-01 4.01e-01w  1
    4744  0.0000000e+00 1.51e+05 5.00e+12  -1.0 6.50e+03    -  1.00e+00 4.95e-01w  1
    4745  0.0000000e+00 2.34e+04 6.23e+12  -1.0 3.49e+03    -  7.95e-01 5.60e-05h 13
    4746  0.0000000e+00 2.34e+04 6.23e+12  -1.0 4.91e+04    -  7.95e-01 5.60e-05h 14
    4747  0.0000000e+00 2.34e+04 6.23e+12  -1.0 4.95e+04    -  7.95e-01 5.60e-05h 14
    4748  0.0000000e+00 2.34e+04 6.23e+12  -1.0 4.95e+04    -  7.95e-01 5.60e-05h 14
    4749  0.0000000e+00 2.34e+04 6.23e+12  -1.0 4.95e+04    -  7.95e-01 5.60e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4750  0.0000000e+00 2.34e+04 6.23e+12  -1.0 4.95e+04    -  7.95e-01 5.60e-05h 14
    4751  0.0000000e+00 2.34e+04 6.23e+12  -1.0 4.94e+04    -  7.95e-01 5.60e-05h 14
    4752  0.0000000e+00 2.34e+04 6.23e+12  -1.0 4.95e+04    -  7.95e-01 5.60e-05h 14
    4753  0.0000000e+00 2.34e+04 6.23e+12  -1.0 4.95e+04    -  7.96e-01 5.60e-05h 14
    4754  0.0000000e+00 2.34e+04 6.23e+12  -1.0 4.94e+04    -  7.96e-01 5.60e-05h 14
    4755  0.0000000e+00 3.42e+05 1.27e+13  -1.0 4.93e+04    -  7.96e-01 4.59e-01w  1
    4756  0.0000000e+00 2.96e+05 6.31e+12  -1.0 4.65e+04    -  7.94e-01 4.01e-01w  1
    4757  0.0000000e+00 1.51e+05 5.00e+12  -1.0 6.47e+03    -  1.00e+00 4.95e-01w  1
    4758  0.0000000e+00 2.34e+04 6.23e+12  -1.0 3.47e+03    -  7.96e-01 5.60e-05h 13
    4759  0.0000000e+00 2.34e+04 6.23e+12  -1.0 4.95e+04    -  7.96e-01 5.60e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4760  0.0000000e+00 2.34e+04 6.23e+12  -1.0 4.94e+04    -  7.96e-01 5.60e-05h 14
    4761  0.0000000e+00 2.34e+04 6.23e+12  -1.0 4.94e+04    -  7.96e-01 5.60e-05h 14
    4762  0.0000000e+00 2.34e+04 6.23e+12  -1.0 4.94e+04    -  7.96e-01 5.60e-05h 14
    4763  0.0000000e+00 2.34e+04 6.23e+12  -1.0 4.95e+04    -  7.96e-01 5.60e-05h 14
    4764  0.0000000e+00 2.34e+04 6.23e+12  -1.0 4.95e+04    -  7.96e-01 5.60e-05h 14
    4765  0.0000000e+00 2.34e+04 6.23e+12  -1.0 4.92e+04    -  7.96e-01 5.60e-05h 14
    4766  0.0000000e+00 2.34e+04 6.23e+12  -1.0 4.94e+04    -  7.96e-01 5.60e-05h 14
    4767  0.0000000e+00 2.34e+04 6.23e+12  -1.0 4.95e+04    -  7.96e-01 5.60e-05h 14
    4768  0.0000000e+00 3.41e+05 1.27e+13  -1.0 4.95e+04    -  7.96e-01 4.59e-01w  1
    4769  0.0000000e+00 2.96e+05 6.32e+12  -1.0 4.65e+04    -  7.95e-01 4.01e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4770  0.0000000e+00 1.51e+05 5.01e+12  -1.0 6.54e+03    -  1.00e+00 4.95e-01w  1
    4771  0.0000000e+00 2.34e+04 6.23e+12  -1.0 3.50e+03    -  7.96e-01 5.60e-05h 13
    4772  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.94e+04    -  7.96e-01 5.60e-05h 14
    4773  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.95e+04    -  7.96e-01 5.60e-05h 14
    4774  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.95e+04    -  7.96e-01 5.60e-05h 14
    4775  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.94e+04    -  7.97e-01 5.60e-05h 14
    4776  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.95e+04    -  7.97e-01 5.60e-05h 14
    4777  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.93e+04    -  7.97e-01 5.60e-05h 14
    4778  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.94e+04    -  7.97e-01 5.60e-05h 14
    4779  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.94e+04    -  7.97e-01 5.60e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4780  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.94e+04    -  7.97e-01 5.60e-05h 14
    4781  0.0000000e+00 3.41e+05 1.27e+13  -1.0 4.94e+04    -  7.97e-01 4.59e-01w  1
    4782  0.0000000e+00 2.95e+05 6.32e+12  -1.0 4.65e+04    -  7.95e-01 4.01e-01w  1
    4783  0.0000000e+00 1.50e+05 5.01e+12  -1.0 6.55e+03    -  1.00e+00 4.95e-01w  1
    4784  0.0000000e+00 2.34e+04 6.24e+12  -1.0 3.51e+03    -  7.97e-01 5.60e-05h 13
    4785  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.94e+04    -  7.97e-01 5.60e-05h 14
    4786  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.94e+04    -  7.97e-01 5.60e-05h 14
    4787  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.94e+04    -  7.97e-01 5.60e-05h 14
    4788  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.93e+04    -  7.97e-01 5.60e-05h 14
    4789  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.94e+04    -  7.97e-01 5.60e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4790  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.93e+04    -  7.97e-01 5.60e-05h 14
    4791  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.93e+04    -  7.97e-01 5.60e-05h 14
    4792  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.93e+04    -  7.98e-01 5.60e-05h 14
    4793  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.94e+04    -  7.98e-01 5.60e-05h 14
    4794  0.0000000e+00 3.41e+05 1.27e+13  -1.0 4.94e+04    -  7.98e-01 4.59e-01w  1
    4795  0.0000000e+00 2.95e+05 6.33e+12  -1.0 4.64e+04    -  7.96e-01 4.01e-01w  1
    4796  0.0000000e+00 1.50e+05 5.02e+12  -1.0 6.56e+03    -  1.00e+00 4.95e-01w  1
    4797  0.0000000e+00 2.34e+04 6.24e+12  -1.0 3.52e+03    -  7.98e-01 5.60e-05h 13
    4798  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.94e+04    -  7.98e-01 5.60e-05h 14
    4799  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.94e+04    -  7.98e-01 5.60e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4800  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.94e+04    -  7.98e-01 5.60e-05h 14
    4801  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.94e+04    -  7.98e-01 5.60e-05h 14
    4802  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.93e+04    -  7.98e-01 5.60e-05h 14
    4803  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.94e+04    -  7.98e-01 5.60e-05h 14
    4804  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.94e+04    -  7.98e-01 5.60e-05h 14
    4805  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.93e+04    -  7.98e-01 5.60e-05h 14
    4806  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.94e+04    -  7.98e-01 5.60e-05h 14
    4807  0.0000000e+00 3.40e+05 1.27e+13  -1.0 4.93e+04    -  7.98e-01 4.59e-01w  1
    4808  0.0000000e+00 2.94e+05 6.34e+12  -1.0 4.63e+04    -  7.96e-01 4.02e-01w  1
    4809  0.0000000e+00 1.50e+05 5.02e+12  -1.0 6.56e+03    -  1.00e+00 4.95e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4810  0.0000000e+00 2.34e+04 6.24e+12  -1.0 3.52e+03    -  7.98e-01 5.60e-05h 13
    4811  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.94e+04    -  7.98e-01 5.60e-05h 14
    4812  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.93e+04    -  7.98e-01 5.60e-05h 14
    4813  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.94e+04    -  7.99e-01 5.60e-05h 14
    4814  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.94e+04    -  7.99e-01 5.60e-05h 14
    4815  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.94e+04    -  7.99e-01 5.60e-05h 14
    4816  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.94e+04    -  7.99e-01 5.60e-05h 14
    4817  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.93e+04    -  7.99e-01 5.60e-05h 14
    4818  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.93e+04    -  7.99e-01 5.60e-05h 14
    4819  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.94e+04    -  7.99e-01 5.60e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4820  0.0000000e+00 3.40e+05 1.27e+13  -1.0 4.93e+04    -  7.99e-01 4.59e-01w  1
    4821  0.0000000e+00 2.94e+05 6.35e+12  -1.0 4.63e+04    -  7.97e-01 4.02e-01w  1
    4822  0.0000000e+00 1.50e+05 5.03e+12  -1.0 6.59e+03    -  1.00e+00 4.95e-01w  1
    4823  0.0000000e+00 2.34e+04 6.24e+12  -1.0 3.53e+03    -  7.99e-01 5.60e-05h 13
    4824  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.94e+04    -  7.99e-01 5.60e-05h 14
    4825  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.93e+04    -  7.99e-01 5.60e-05h 14
    4826  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.94e+04    -  7.99e-01 5.60e-05h 14
    4827  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.93e+04    -  7.99e-01 5.60e-05h 14
    4828  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.93e+04    -  7.99e-01 5.60e-05h 14
    4829  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.93e+04    -  7.99e-01 5.60e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4830  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.93e+04    -  7.99e-01 5.60e-05h 14
    4831  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.93e+04    -  8.00e-01 5.60e-05h 14
    4832  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.93e+04    -  8.00e-01 5.60e-05h 14
    4833  0.0000000e+00 3.39e+05 1.27e+13  -1.0 4.93e+04    -  8.00e-01 4.59e-01w  1
    4834  0.0000000e+00 2.93e+05 6.36e+12  -1.0 4.63e+04    -  7.97e-01 4.02e-01w  1
    4835  0.0000000e+00 1.49e+05 5.03e+12  -1.0 6.62e+03    -  1.00e+00 4.95e-01w  1
    4836  0.0000000e+00 2.34e+04 6.24e+12  -1.0 3.55e+03    -  8.00e-01 5.60e-05h 13
    4837  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.92e+04    -  8.00e-01 5.60e-05h 14
    4838  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.93e+04    -  8.00e-01 5.60e-05h 14
    4839  0.0000000e+00 2.34e+04 6.24e+12  -1.0 4.93e+04    -  8.00e-01 5.60e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4840  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.93e+04    -  8.00e-01 5.60e-05h 14
    4841  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.00e-01 5.60e-05h 14
    4842  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.93e+04    -  8.00e-01 5.60e-05h 14
    4843  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.93e+04    -  8.00e-01 5.60e-05h 14
    4844  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.93e+04    -  8.00e-01 5.60e-05h 14
    4845  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.93e+04    -  8.00e-01 5.60e-05h 14
    4846  0.0000000e+00 3.39e+05 1.27e+13  -1.0 4.93e+04    -  8.00e-01 4.59e-01w  1
    4847  0.0000000e+00 2.93e+05 6.37e+12  -1.0 4.62e+04    -  7.98e-01 4.02e-01w  1
    4848  0.0000000e+00 1.49e+05 5.03e+12  -1.0 6.63e+03    -  1.00e+00 4.95e-01w  1
    4849  0.0000000e+00 2.33e+04 6.24e+12  -1.0 3.55e+03    -  8.00e-01 5.60e-05h 13
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4850  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.93e+04    -  8.00e-01 5.60e-05h 14
    4851  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.93e+04    -  8.00e-01 5.60e-05h 14
    4852  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.01e-01 5.60e-05h 14
    4853  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.01e-01 5.60e-05h 14
    4854  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.01e-01 5.60e-05h 14
    4855  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.01e-01 5.60e-05h 14
    4856  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.91e+04    -  8.01e-01 5.60e-05h 14
    4857  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.93e+04    -  8.01e-01 5.60e-05h 14
    4858  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.01e-01 5.60e-05h 14
    4859  0.0000000e+00 3.39e+05 1.27e+13  -1.0 4.92e+04    -  8.01e-01 4.59e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4860  0.0000000e+00 2.93e+05 6.38e+12  -1.0 4.61e+04    -  7.98e-01 4.02e-01w  1
    4861  0.0000000e+00 1.49e+05 5.04e+12  -1.0 6.63e+03    -  1.00e+00 4.95e-01w  1
    4862  0.0000000e+00 2.33e+04 6.24e+12  -1.0 3.56e+03    -  8.01e-01 5.60e-05h 13
    4863  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.93e+04    -  8.01e-01 5.60e-05h 14
    4864  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.93e+04    -  8.01e-01 5.60e-05h 14
    4865  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.93e+04    -  8.01e-01 5.60e-05h 14
    4866  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.91e+04    -  8.01e-01 5.60e-05h 14
    4867  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.01e-01 5.60e-05h 14
    4868  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.01e-01 5.60e-05h 14
    4869  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.93e+04    -  8.01e-01 5.60e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4870  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.93e+04    -  8.02e-01 5.60e-05h 14
    4871  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.02e-01 5.60e-05h 14
    4872  0.0000000e+00 3.38e+05 1.27e+13  -1.0 4.93e+04    -  8.02e-01 4.59e-01w  1
    4873  0.0000000e+00 2.92e+05 6.39e+12  -1.0 4.61e+04    -  7.99e-01 4.02e-01w  1
    4874  0.0000000e+00 1.49e+05 5.04e+12  -1.0 6.67e+03    -  1.00e+00 4.95e-01w  1
    4875  0.0000000e+00 2.33e+04 6.24e+12  -1.0 3.58e+03    -  8.02e-01 5.60e-05h 13
    4876  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.93e+04    -  8.02e-01 5.60e-05h 14
    4877  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.02e-01 5.60e-05h 14
    4878  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.93e+04    -  8.02e-01 5.60e-05h 14
    4879  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.91e+04    -  8.02e-01 5.60e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4880  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.89e+04    -  8.02e-01 5.60e-05h 14
    4881  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.02e-01 5.60e-05h 14
    4882  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.02e-01 5.60e-05h 14
    4883  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.02e-01 5.60e-05h 14
    4884  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.91e+04    -  8.02e-01 5.60e-05h 14
    4885  0.0000000e+00 3.38e+05 1.27e+13  -1.0 4.92e+04    -  8.02e-01 4.59e-01w  1
    4886  0.0000000e+00 2.92e+05 6.40e+12  -1.0 4.61e+04    -  7.99e-01 4.03e-01w  1
    4887  0.0000000e+00 1.48e+05 5.05e+12  -1.0 6.68e+03    -  1.00e+00 4.95e-01w  1
    4888  0.0000000e+00 2.33e+04 6.24e+12  -1.0 3.58e+03    -  8.02e-01 5.60e-05h 13
    4889  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.02e-01 5.60e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4890  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.91e+04    -  8.03e-01 5.60e-05h 14
    4891  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.91e+04    -  8.03e-01 5.60e-05h 14
    4892  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.03e-01 5.60e-05h 14
    4893  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.03e-01 5.60e-05h 14
    4894  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.03e-01 5.60e-05h 14
    4895  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.91e+04    -  8.03e-01 5.60e-05h 14
    4896  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.03e-01 5.60e-05h 14
    4897  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.91e+04    -  8.03e-01 5.60e-05h 14
    4898  0.0000000e+00 3.37e+05 1.27e+13  -1.0 4.92e+04    -  8.03e-01 4.59e-01w  1
    4899  0.0000000e+00 2.91e+05 6.41e+12  -1.0 4.60e+04    -  8.00e-01 4.03e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4900  0.0000000e+00 1.48e+05 5.05e+12  -1.0 6.70e+03    -  1.00e+00 4.95e-01w  1
    4901  0.0000000e+00 2.33e+04 6.24e+12  -1.0 3.59e+03    -  8.03e-01 5.60e-05h 13
    4902  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.03e-01 5.60e-05h 14
    4903  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.03e-01 5.60e-05h 14
    4904  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.03e-01 5.60e-05h 14
    4905  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.03e-01 5.60e-05h 14
    4906  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.03e-01 5.60e-05h 14
    4907  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.03e-01 5.60e-05h 14
    4908  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.91e+04    -  8.03e-01 5.60e-05h 14
    4909  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.04e-01 5.60e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4910  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.91e+04    -  8.04e-01 5.60e-05h 14
    4911  0.0000000e+00 3.37e+05 1.27e+13  -1.0 4.91e+04    -  8.04e-01 4.59e-01w  1
    4912  0.0000000e+00 2.91e+05 6.42e+12  -1.0 4.60e+04    -  8.00e-01 4.03e-01w  1
    4913  0.0000000e+00 1.48e+05 5.06e+12  -1.0 6.71e+03    -  1.00e+00 4.95e-01w  1
    4914  0.0000000e+00 2.33e+04 6.24e+12  -1.0 3.59e+03    -  8.04e-01 5.60e-05h 13
    4915  0.0000000e+00 2.33e+04 6.24e+12  -1.0 4.92e+04    -  8.04e-01 5.60e-05h 14
    4916  0.0000000e+00 2.33e+04 6.25e+12  -1.0 4.92e+04    -  8.04e-01 5.60e-05h 14
    4917  0.0000000e+00 2.33e+04 6.25e+12  -1.0 4.91e+04    -  8.04e-01 5.60e-05h 14
    4918  0.0000000e+00 2.33e+04 6.25e+12  -1.0 4.92e+04    -  8.04e-01 5.60e-05h 14
    4919  0.0000000e+00 2.33e+04 6.25e+12  -1.0 4.91e+04    -  8.04e-01 5.60e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4920  0.0000000e+00 2.33e+04 6.25e+12  -1.0 4.91e+04    -  8.04e-01 5.60e-05h 14
    4921  0.0000000e+00 2.33e+04 6.25e+12  -1.0 4.91e+04    -  8.04e-01 5.60e-05h 14
    4922  0.0000000e+00 2.33e+04 6.25e+12  -1.0 4.91e+04    -  8.04e-01 5.60e-05h 14
    4923  0.0000000e+00 2.33e+04 6.25e+12  -1.0 4.91e+04    -  8.04e-01 5.60e-05h 14
    4924  0.0000000e+00 3.37e+05 1.27e+13  -1.0 4.92e+04    -  8.04e-01 4.59e-01w  1
    4925  0.0000000e+00 2.90e+05 6.43e+12  -1.0 4.59e+04    -  8.01e-01 4.03e-01w  1
    4926  0.0000000e+00 1.48e+05 5.06e+12  -1.0 6.73e+03    -  1.00e+00 4.95e-01w  1
    4927  0.0000000e+00 2.33e+04 6.25e+12  -1.0 3.61e+03    -  8.04e-01 5.60e-05h 13
    4928  0.0000000e+00 2.33e+04 6.25e+12  -1.0 4.91e+04    -  8.04e-01 5.60e-05h 14
    4929  0.0000000e+00 2.33e+04 6.25e+12  -1.0 4.91e+04    -  8.05e-01 5.60e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4930  0.0000000e+00 2.33e+04 6.25e+12  -1.0 4.91e+04    -  8.05e-01 5.60e-05h 14
    4931  0.0000000e+00 2.33e+04 6.25e+12  -1.0 4.91e+04    -  8.05e-01 5.60e-05h 14
    4932  0.0000000e+00 2.33e+04 6.25e+12  -1.0 4.91e+04    -  8.05e-01 5.60e-05h 14
    4933  0.0000000e+00 2.33e+04 6.25e+12  -1.0 4.91e+04    -  8.05e-01 5.60e-05h 14
    4934  0.0000000e+00 2.33e+04 6.25e+12  -1.0 4.91e+04    -  8.05e-01 5.60e-05h 14
    4935  0.0000000e+00 2.33e+04 6.25e+12  -1.0 4.91e+04    -  8.05e-01 5.60e-05h 14
    4936  0.0000000e+00 2.33e+04 6.25e+12  -1.0 4.91e+04    -  8.05e-01 5.60e-05h 14
    4937  0.0000000e+00 3.36e+05 1.27e+13  -1.0 4.90e+04    -  8.05e-01 4.59e-01w  1
    4938  0.0000000e+00 2.90e+05 6.44e+12  -1.0 4.59e+04    -  8.01e-01 4.03e-01w  1
    4939  0.0000000e+00 1.47e+05 5.07e+12  -1.0 6.72e+03    -  1.00e+00 4.95e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4940  0.0000000e+00 2.32e+04 6.25e+12  -1.0 3.61e+03    -  8.05e-01 5.60e-05h 13
    4941  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.91e+04    -  8.05e-01 5.60e-05h 14
    4942  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.91e+04    -  8.05e-01 5.60e-05h 14
    4943  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.91e+04    -  8.05e-01 5.60e-05h 14
    4944  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.91e+04    -  8.05e-01 5.60e-05h 14
    4945  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.91e+04    -  8.05e-01 5.60e-05h 14
    4946  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.91e+04    -  8.05e-01 5.60e-05h 14
    4947  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.91e+04    -  8.05e-01 5.60e-05h 14
    4948  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.91e+04    -  8.06e-01 5.60e-05h 14
    4949  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.91e+04    -  8.06e-01 5.60e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4950  0.0000000e+00 3.36e+05 1.27e+13  -1.0 4.91e+04    -  8.06e-01 4.59e-01w  1
    4951  0.0000000e+00 2.89e+05 6.45e+12  -1.0 4.59e+04    -  8.01e-01 4.03e-01w  1
    4952  0.0000000e+00 1.47e+05 5.07e+12  -1.0 6.76e+03    -  1.00e+00 4.95e-01w  1
    4953  0.0000000e+00 2.32e+04 6.25e+12  -1.0 3.63e+03    -  8.06e-01 5.60e-05h 13
    4954  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.90e+04    -  8.06e-01 5.60e-05h 14
    4955  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.91e+04    -  8.06e-01 5.60e-05h 14
    4956  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.91e+04    -  8.06e-01 5.60e-05h 14
    4957  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.91e+04    -  8.06e-01 5.60e-05h 14
    4958  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.91e+04    -  8.06e-01 5.60e-05h 14
    4959  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.91e+04    -  8.06e-01 5.60e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4960  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.90e+04    -  8.06e-01 5.60e-05h 14
    4961  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.89e+04    -  8.06e-01 5.60e-05h 14
    4962  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.91e+04    -  8.06e-01 5.60e-05h 14
    4963  0.0000000e+00 3.35e+05 1.27e+13  -1.0 4.91e+04    -  8.06e-01 4.59e-01w  1
    4964  0.0000000e+00 2.89e+05 6.46e+12  -1.0 4.58e+04    -  8.02e-01 4.04e-01w  1
    4965  0.0000000e+00 1.47e+05 5.07e+12  -1.0 6.78e+03    -  1.00e+00 4.95e-01w  1
    4966  0.0000000e+00 2.32e+04 6.25e+12  -1.0 3.63e+03    -  8.06e-01 5.60e-05h 13
    4967  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.90e+04    -  8.06e-01 5.60e-05h 14
    4968  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.91e+04    -  8.07e-01 5.60e-05h 14
    4969  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.91e+04    -  8.06e-01 5.60e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4970  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.90e+04    -  8.07e-01 5.60e-05h 14
    4971  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.91e+04    -  8.07e-01 5.60e-05h 14
    4972  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.91e+04    -  8.07e-01 5.60e-05h 14
    4973  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.90e+04    -  8.07e-01 5.60e-05h 14
    4974  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.90e+04    -  8.07e-01 5.60e-05h 14
    4975  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.91e+04    -  8.07e-01 5.60e-05h 14
    4976  0.0000000e+00 3.35e+05 1.27e+13  -1.0 4.90e+04    -  8.07e-01 4.59e-01w  1
    4977  0.0000000e+00 2.88e+05 6.47e+12  -1.0 4.58e+04    -  8.02e-01 4.04e-01w  1
    4978  0.0000000e+00 1.47e+05 5.08e+12  -1.0 6.79e+03    -  1.00e+00 4.95e-01w  1
    4979  0.0000000e+00 2.32e+04 6.25e+12  -1.0 3.64e+03    -  8.07e-01 5.60e-05h 13
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4980  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.90e+04    -  8.07e-01 5.60e-05h 14
    4981  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.90e+04    -  8.07e-01 5.60e-05h 14
    4982  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.90e+04    -  8.07e-01 5.60e-05h 14
    4983  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.90e+04    -  8.07e-01 5.60e-05h 14
    4984  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.90e+04    -  8.07e-01 5.60e-05h 14
    4985  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.90e+04    -  8.07e-01 5.60e-05h 14
    4986  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.89e+04    -  8.07e-01 5.60e-05h 14
    4987  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.89e+04    -  8.08e-01 5.60e-05h 14
    4988  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.90e+04    -  8.08e-01 5.60e-05h 14
    4989  0.0000000e+00 3.35e+05 1.27e+13  -1.0 4.90e+04    -  8.08e-01 4.59e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    4990  0.0000000e+00 2.88e+05 6.48e+12  -1.0 4.57e+04    -  8.03e-01 4.04e-01w  1
    4991  0.0000000e+00 1.46e+05 5.08e+12  -1.0 6.81e+03    -  1.00e+00 4.95e-01w  1
    4992  0.0000000e+00 2.32e+04 6.25e+12  -1.0 3.65e+03    -  8.08e-01 5.60e-05h 13
    4993  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.90e+04    -  8.08e-01 5.60e-05h 14
    4994  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.90e+04    -  8.08e-01 5.60e-05h 14
    4995  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.90e+04    -  8.08e-01 5.60e-05h 14
    4996  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.90e+04    -  8.08e-01 5.60e-05h 14
    4997  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.90e+04    -  8.08e-01 5.60e-05h 14
    4998  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.90e+04    -  8.08e-01 5.60e-05h 14
    4999  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.90e+04    -  8.08e-01 5.60e-05h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    5000  0.0000000e+00 2.32e+04 6.25e+12  -1.0 4.90e+04    -  8.08e-01 5.60e-05h 14
    
    Number of Iterations....: 5000
    
                                       (scaled)                 (unscaled)
    Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    Dual infeasibility......:   6.2513164942996592e+12    6.2513164942996592e+12
    Constraint violation....:   2.0953258594124179e-02    2.3187475840772902e+04
    Complementarity.........:   1.0000000000000000e+09    1.0000000000000000e+09
    Overall NLP error.......:   7.9330703659079532e+01    6.2513164942996592e+12
    
    
    Number of objective function evaluations             = 50631
    Number of objective gradient evaluations             = 4986
    Number of equality constraint evaluations            = 50653
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 5048
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 5000
    Total CPU secs in IPOPT (w/o function evaluations)   =      7.623
    Total CPU secs in NLP function evaluations           =      1.370
    
    EXIT: Maximum Number of Iterations Exceeded.
    WARNING: Loading a SolverResults object with a warning status into
        model=unknown;
            message from solver=Ipopt 3.12.8\x3a Maximum Number of Iterations
            Exceeded.


.. code:: ipython3

    # For testing purposes
    from pyomo.environ import TerminationCondition
    assert results.solver.termination_condition == TerminationCondition.optimal


::


    ---------------------------------------------------------------------------

    AssertionError                            Traceback (most recent call last)

    <ipython-input-38-27ff7383fadd> in <module>
          1 # For testing purposes
          2 from pyomo.environ import TerminationCondition
    ----> 3 assert results.solver.termination_condition == TerminationCondition.optimal
    

    AssertionError: 


Analyze the results of the square problem
-----------------------------------------

What is the total operating cost?

.. code:: ipython3

    print('operating cost = $', value(m.fs.operating_cost))


.. parsed-literal::

    operating cost = $ 404944.6501353814


For this operating cost, what is the amount of benzene we are able to
produce and what purity we are able to achieve?

.. code:: ipython3

    m.fs.F102.report()
    
    print()
    print('benzene purity = ', value(m.fs.purity))


.. parsed-literal::

    
    ====================================================================================
    Unit : fs.F102                                                             Time: 0.0
    ------------------------------------------------------------------------------------
        Unit Performance
    
        Variables: 
    
        Key             : Value       : Fixed : Bounds
              Heat Duty :      3878.5 : False : (None, None)
        Pressure Change : -2.0000e+05 :  True : (None, None)
    
    ------------------------------------------------------------------------------------
        Stream Table
                                                   Inlet    Vapor Outlet  Liquid Outlet
        flow_mol_phase_comp ('Liq', 'benzene')     0.18713   1.0000e-08       0.13938  
        flow_mol_phase_comp ('Liq', 'toluene')    0.081976   1.0000e-08      0.071799  
        flow_mol_phase_comp ('Liq', 'hydrogen') 2.6910e-07   1.0000e-08    2.1118e-07  
        flow_mol_phase_comp ('Liq', 'methane')  2.6910e-07   1.0000e-08    2.1118e-07  
        flow_mol_phase_comp ('Vap', 'benzene')  1.0000e-08     0.047743    1.0000e-08  
        flow_mol_phase_comp ('Vap', 'toluene')  1.0000e-08     0.010177    1.0000e-08  
        flow_mol_phase_comp ('Vap', 'hydrogen') 1.0000e-08   6.7920e-08    1.0000e-08  
        flow_mol_phase_comp ('Vap', 'methane')  1.0000e-08   6.7920e-08    1.0000e-08  
        temperature                                 325.00       375.00        375.00  
        pressure                                3.5000e+05   1.5000e+05    1.5000e+05  
    ====================================================================================
    
    benzene purity =  0.8242960919828991


Next, let’s look at how much benzene we are loosing with the light gases
out of F101. IDAES has tools for creating stream tables based on the
``Arcs`` and/or ``Ports`` in a flowsheet. Let us create and print a
simple stream table showing the stream leaving the reactor and the vapor
stream from F101.

.. raw:: html

   <div class="alert alert-block alert-info">

Inline Exercise: How much benzene are we losing in the F101 vapor outlet
stream?

.. raw:: html

   </div>

.. code:: ipython3

    from idaes.core.util.tables import create_stream_table_dataframe, stream_table_dataframe_to_string
    
    st = create_stream_table_dataframe({"Reactor": m.fs.s05, "Light Gases": m.fs.s06})
    st




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>Reactor</th>
          <th>Light Gases</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td>flow_mol_phase_comp ('Liq', 'benzene')</td>
          <td>0.167639</td>
          <td>1.000000e-08</td>
        </tr>
        <tr>
          <td>flow_mol_phase_comp ('Liq', 'toluene')</td>
          <td>0.023130</td>
          <td>1.000000e-08</td>
        </tr>
        <tr>
          <td>flow_mol_phase_comp ('Liq', 'hydrogen')</td>
          <td>0.207166</td>
          <td>1.000000e-08</td>
        </tr>
        <tr>
          <td>flow_mol_phase_comp ('Liq', 'methane')</td>
          <td>0.487846</td>
          <td>1.000000e-08</td>
        </tr>
        <tr>
          <td>flow_mol_phase_comp ('Vap', 'benzene')</td>
          <td>0.154129</td>
          <td>1.346430e-01</td>
        </tr>
        <tr>
          <td>flow_mol_phase_comp ('Vap', 'toluene')</td>
          <td>0.079048</td>
          <td>2.020218e-02</td>
        </tr>
        <tr>
          <td>flow_mol_phase_comp ('Vap', 'hydrogen')</td>
          <td>0.222913</td>
          <td>4.300786e-01</td>
        </tr>
        <tr>
          <td>flow_mol_phase_comp ('Vap', 'methane')</td>
          <td>0.682373</td>
          <td>1.170219e+00</td>
        </tr>
        <tr>
          <td>temperature</td>
          <td>830.063048</td>
          <td>3.250000e+02</td>
        </tr>
        <tr>
          <td>pressure</td>
          <td>349999.999191</td>
          <td>3.500000e+05</td>
        </tr>
      </tbody>
    </table>
    </div>



.. raw:: html

   <div class="alert alert-block alert-info">

Inline Exercise: You can querry additional variables here if you like.

Use Shift+Enter to run the cell once you have typed in your code.

.. raw:: html

   </div>
