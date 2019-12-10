
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

    2019-12-10 09:05:32 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 1 Complete.
    2019-12-10 09:05:32 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 2 Complete.
    2019-12-10 09:05:32 - INFO - idaes.core.unit_model - fs.H101 Initialisation Complete.
    2019-12-10 09:05:32 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 1 Complete.
    2019-12-10 09:05:32 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 2 Complete.
    2019-12-10 09:05:32 - INFO - idaes.core.unit_model - fs.R101 Initialisation Complete.
    2019-12-10 09:05:32 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 1 Complete.
    2019-12-10 09:05:32 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 2 Complete.
    2019-12-10 09:05:32 - INFO - idaes.core.unit_model - fs.F101 Initialisation Complete.
    2019-12-10 09:05:32 - INFO - idaes.unit_models.separator - fs.S101 Initialisation Complete.
    2019-12-10 09:05:32 - INFO - idaes.core.unit_model - fs.F102 Initialisation Step 1 Complete.
    2019-12-10 09:05:32 - INFO - idaes.core.unit_model - fs.F102 Initialisation Step 2 Complete.
    2019-12-10 09:05:32 - INFO - idaes.core.unit_model - fs.F102 Initialisation Complete.
    2019-12-10 09:05:32 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 1 Complete.
    2019-12-10 09:05:32 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 2 Complete.
    2019-12-10 09:05:32 - INFO - idaes.core.unit_model - fs.C101 Initialisation Complete.
    2019-12-10 09:05:33 - INFO - idaes.unit_models.mixer - fs.M101 Initialisation Complete.
    2019-12-10 09:05:33 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 1 Complete.
    2019-12-10 09:05:33 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 2 Complete.
    2019-12-10 09:05:33 - INFO - idaes.core.unit_model - fs.H101 Initialisation Complete.
    2019-12-10 09:05:33 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 1 Complete.
    2019-12-10 09:05:33 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 2 Complete.
    2019-12-10 09:05:33 - INFO - idaes.core.unit_model - fs.R101 Initialisation Complete.
    2019-12-10 09:05:33 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 1 Complete.
    2019-12-10 09:05:33 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 2 Complete.
    2019-12-10 09:05:33 - INFO - idaes.core.unit_model - fs.F101 Initialisation Complete.
    2019-12-10 09:05:33 - INFO - idaes.unit_models.separator - fs.S101 Initialisation Complete.
    2019-12-10 09:05:33 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 1 Complete.
    2019-12-10 09:05:33 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 2 Complete.
    2019-12-10 09:05:33 - INFO - idaes.core.unit_model - fs.C101 Initialisation Complete.
    2019-12-10 09:05:33 - INFO - idaes.unit_models.mixer - fs.M101 Initialisation Complete.
    2019-12-10 09:05:33 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 1 Complete.
    2019-12-10 09:05:33 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 2 Complete.
    2019-12-10 09:05:33 - INFO - idaes.core.unit_model - fs.H101 Initialisation Complete.
    2019-12-10 09:05:33 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 1 Complete.
    2019-12-10 09:05:33 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 2 Complete.
    2019-12-10 09:05:33 - INFO - idaes.core.unit_model - fs.R101 Initialisation Complete.
    2019-12-10 09:05:33 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 1 Complete.
    2019-12-10 09:05:34 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 2 Complete.
    2019-12-10 09:05:34 - INFO - idaes.core.unit_model - fs.F101 Initialisation Complete.
    2019-12-10 09:05:34 - INFO - idaes.unit_models.separator - fs.S101 Initialisation Complete.
    2019-12-10 09:05:34 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 1 Complete.
    2019-12-10 09:05:34 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 2 Complete.
    2019-12-10 09:05:34 - INFO - idaes.core.unit_model - fs.C101 Initialisation Complete.
    2019-12-10 09:05:34 - INFO - idaes.unit_models.mixer - fs.M101 Initialisation Complete.
    2019-12-10 09:05:34 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 1 Complete.
    2019-12-10 09:05:34 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 2 Complete.
    2019-12-10 09:05:34 - INFO - idaes.core.unit_model - fs.H101 Initialisation Complete.
    2019-12-10 09:05:34 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 1 Complete.
    2019-12-10 09:05:34 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 2 Complete.
    2019-12-10 09:05:34 - INFO - idaes.core.unit_model - fs.R101 Initialisation Complete.
    2019-12-10 09:05:34 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 1 Complete.
    2019-12-10 09:05:34 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 2 Complete.
    2019-12-10 09:05:34 - INFO - idaes.core.unit_model - fs.F101 Initialisation Complete.
    2019-12-10 09:05:34 - INFO - idaes.unit_models.separator - fs.S101 Initialisation Complete.
    2019-12-10 09:05:34 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 1 Complete.
    2019-12-10 09:05:34 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 2 Complete.
    2019-12-10 09:05:34 - INFO - idaes.core.unit_model - fs.C101 Initialisation Complete.
    2019-12-10 09:05:34 - INFO - idaes.unit_models.mixer - fs.M101 Initialisation Complete.
    2019-12-10 09:05:34 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 1 Complete.
    2019-12-10 09:05:34 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 2 Complete.
    2019-12-10 09:05:34 - INFO - idaes.core.unit_model - fs.H101 Initialisation Complete.
    2019-12-10 09:05:34 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 1 Complete.
    2019-12-10 09:05:34 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 2 Complete.
    2019-12-10 09:05:34 - INFO - idaes.core.unit_model - fs.R101 Initialisation Complete.
    2019-12-10 09:05:35 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 1 Complete.
    2019-12-10 09:05:35 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 2 Complete.
    2019-12-10 09:05:35 - INFO - idaes.core.unit_model - fs.F101 Initialisation Complete.
    2019-12-10 09:05:35 - INFO - idaes.unit_models.separator - fs.S101 Initialisation Complete.
    2019-12-10 09:05:35 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 1 Complete.
    2019-12-10 09:05:35 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 2 Complete.
    2019-12-10 09:05:35 - INFO - idaes.core.unit_model - fs.C101 Initialisation Complete.
    2019-12-10 09:05:35 - INFO - idaes.unit_models.mixer - fs.M101 Initialisation Complete.
    2019-12-10 09:05:35 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 1 Complete.
    2019-12-10 09:05:35 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 2 Complete.
    2019-12-10 09:05:35 - INFO - idaes.core.unit_model - fs.H101 Initialisation Complete.
    2019-12-10 09:05:35 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 1 Complete.
    2019-12-10 09:05:35 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 2 Complete.
    2019-12-10 09:05:35 - INFO - idaes.core.unit_model - fs.R101 Initialisation Complete.
    2019-12-10 09:05:35 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 1 Complete.
    2019-12-10 09:05:35 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 2 Complete.
    2019-12-10 09:05:35 - INFO - idaes.core.unit_model - fs.F101 Initialisation Complete.
    2019-12-10 09:05:35 - INFO - idaes.unit_models.separator - fs.S101 Initialisation Complete.
    2019-12-10 09:05:35 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 1 Complete.
    2019-12-10 09:05:35 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 2 Complete.
    2019-12-10 09:05:35 - INFO - idaes.core.unit_model - fs.C101 Initialisation Complete.
    2019-12-10 09:05:35 - INFO - idaes.unit_models.mixer - fs.M101 Initialisation Complete.
    2019-12-10 09:05:35 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 1 Complete.
    2019-12-10 09:05:35 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 2 Complete.
    2019-12-10 09:05:35 - INFO - idaes.core.unit_model - fs.H101 Initialisation Complete.
    2019-12-10 09:05:35 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 1 Complete.
    2019-12-10 09:05:35 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 2 Complete.
    2019-12-10 09:05:35 - INFO - idaes.core.unit_model - fs.R101 Initialisation Complete.
    2019-12-10 09:05:35 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 1 Complete.
    2019-12-10 09:05:36 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 2 Complete.
    2019-12-10 09:05:36 - INFO - idaes.core.unit_model - fs.F101 Initialisation Complete.
    2019-12-10 09:05:36 - INFO - idaes.unit_models.separator - fs.S101 Initialisation Complete.
    2019-12-10 09:05:36 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 1 Complete.
    2019-12-10 09:05:36 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 2 Complete.
    2019-12-10 09:05:36 - INFO - idaes.core.unit_model - fs.C101 Initialisation Complete.
    2019-12-10 09:05:36 - INFO - idaes.unit_models.mixer - fs.M101 Initialisation Complete.
    WARNING: Wegstein failed to converge in 5 iterations
    2019-12-10 09:05:36 - INFO - idaes.core.unit_model - fs.F102 Initialisation Step 1 Complete.
    2019-12-10 09:05:36 - INFO - idaes.core.unit_model - fs.F102 Initialisation Step 2 Complete.
    2019-12-10 09:05:36 - INFO - idaes.core.unit_model - fs.F102 Initialisation Complete.


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

    Ipopt 3.12.13: tol=1e-06
    max_iter=5000
    
    
    ******************************************************************************
    This program contains Ipopt, a library for large-scale nonlinear optimization.
     Ipopt is released as open source code under the Eclipse Public License (EPL).
             For more information visit http://projects.coin-or.org/Ipopt
    
    This version of Ipopt was compiled from source code available at
        https://github.com/IDAES/Ipopt as part of the Institute for the Design of
        Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
        Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    
    This version of Ipopt was compiled using HSL, a collection of Fortran codes
        for large-scale scientific computation.  All technical papers, sales and
        publicity material resulting from use of the HSL codes within IPOPT must
        contain the following acknowledgement:
            HSL, a collection of Fortran codes for large-scale scientific
            computation. See http://www.hsl.rl.ac.uk.
    ******************************************************************************
    
    This is Ipopt version 3.12.13, running with linear solver ma27.
    
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
       1  0.0000000e+00 3.05e+04 1.21e+02  -1.0 2.26e+04    -  9.15e-01 1.86e-01h  1
       2  0.0000000e+00 1.41e+04 5.72e+02  -1.0 1.66e+04    -  9.89e-01 4.77e-01h  1
       3  0.0000000e+00 1.38e+04 9.90e+04  -1.0 1.60e+04    -  9.90e-01 1.55e-02h  6
       4  0.0000000e+00 1.37e+04 3.00e+05  -1.0 1.35e+04    -  1.00e+00 7.75e-03h  7
       5  0.0000000e+00 1.36e+04 7.03e+05  -1.0 1.16e+04    -  1.00e+00 7.76e-03h  7
       6  0.0000000e+00 1.36e+04 1.51e+06  -1.0 1.27e+04    -  1.00e+00 3.88e-03h  8
       7  0.0000000e+00 1.35e+04 3.11e+06  -1.0 1.46e+04    -  1.00e+00 3.89e-03h  8
       8  0.0000000e+00 1.35e+04 6.28e+06  -1.0 1.77e+04    -  1.00e+00 1.95e-03h  9
       9  0.0000000e+00 1.35e+04 1.25e+07  -1.0 2.85e+04    -  1.00e+00 1.97e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  0.0000000e+00 1.35e+04 2.44e+07  -1.0 5.98e+04    -  1.00e+00 9.94e-04h 10
      11  0.0000000e+00 1.34e+04 4.69e+07  -1.0 9.91e+04    -  1.00e+00 4.83e-04h 11
      12  0.0000000e+00 1.34e+04 8.93e+07  -1.0 1.24e+05    -  1.00e+00 4.07e-04h 11
      13  0.0000000e+00 7.48e+05 1.01e+08  -1.0 1.28e+05    -  1.00e+00 4.07e-01w  1
      14  0.0000000e+00 3.93e+05 1.91e+08  -1.0 2.70e+04    -  6.78e-01 4.91e-01w  1
      15  0.0000000e+00 2.54e+05 1.96e+08  -1.0 6.63e+04    -  4.89e-01 4.37e-01w  1
      16  0.0000000e+00 1.34e+04 1.70e+08  -1.0 1.77e+04    -  1.00e+00 3.97e-04h 10
      17  0.0000000e+00 1.34e+04 3.23e+08  -1.0 1.23e+05    -  1.00e+00 4.08e-04h 11
      18  0.0000000e+00 1.34e+04 6.15e+08  -1.0 1.15e+05    -  1.00e+00 4.29e-04h 11
      19  0.0000000e+00 1.34e+04 1.12e+09  -1.0 1.08e+05    -  9.03e-01 4.51e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      20  0.0000000e+00 1.34e+04 2.15e+09  -1.0 1.03e+05    -  1.00e+00 4.68e-04h 11
      21  0.0000000e+00 1.34e+04 2.92e+09  -1.0 9.87e+04    -  3.93e-01 4.81e-04h 11
      22  0.0000000e+00 1.34e+04 5.60e+09  -1.0 9.73e+04    -  1.00e+00 4.86e-04h 11
      23  0.0000000e+00 1.34e+04 7.69e+09  -1.0 9.53e+04    -  4.07e-01 4.93e-04h 11
      24  0.0000000e+00 1.34e+04 1.48e+10  -1.0 9.45e+04    -  1.00e+00 4.96e-04h 11
      25  0.0000000e+00 1.34e+04 2.04e+10  -1.0 9.36e+04    -  4.13e-01 4.99e-04h 11
      26  0.0000000e+00 8.25e+05 1.91e+10  -1.0 9.31e+04    -  1.00e+00 5.12e-01w  1
      27  0.0000000e+00 4.02e+05 4.76e+10  -1.0 2.28e+04    -  6.55e-01 5.20e-01w  1
      28  0.0000000e+00 1.09e+05 1.49e+12  -1.0 1.14e+04    -  1.97e-01 7.36e-01w  1
      29  0.0000000e+00 1.34e+04 3.92e+10  -1.0 4.00e+03    -  1.00e+00 5.00e-04h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      30  0.0000000e+00 1.34e+04 5.42e+10  -1.0 9.27e+04    -  4.16e-01 5.01e-04h 11
      31  0.0000000e+00 1.34e+04 1.04e+11  -1.0 9.24e+04    -  1.00e+00 5.02e-04h 11
      32  0.0000000e+00 1.34e+04 1.44e+11  -1.0 9.22e+04    -  4.18e-01 5.02e-04h 11
      33  0.0000000e+00 1.33e+04 2.77e+11  -1.0 9.19e+04    -  1.00e+00 5.03e-04h 11
      34  0.0000000e+00 1.33e+04 3.84e+11  -1.0 9.19e+04    -  4.19e-01 5.03e-04h 11
      35  0.0000000e+00 1.33e+04 7.39e+11  -1.0 9.16e+04    -  1.00e+00 5.03e-04h 11
      36  0.0000000e+00 1.33e+04 1.03e+12  -1.0 9.16e+04    -  4.20e-01 5.03e-04h 11
      37  0.0000000e+00 1.33e+04 1.97e+12  -1.0 9.14e+04    -  1.00e+00 5.03e-04h 11
      38  0.0000000e+00 1.33e+04 2.74e+12  -1.0 9.13e+04    -  4.21e-01 5.03e-04h 11
      39  0.0000000e+00 8.13e+05 2.55e+12  -1.0 9.11e+04    -  1.00e+00 5.15e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      40  0.0000000e+00 3.78e+05 6.76e+12  -1.0 2.26e+04    -  6.54e-01 5.43e-01w  1
      41  0.0000000e+00 2.66e+05 2.46e+13  -1.0 7.80e+03    -  1.20e-01 2.98e-01w  1
      42  0.0000000e+00 1.33e+04 5.26e+12  -1.0 7.84e+03    -  1.00e+00 5.03e-04h 10
      43  0.0000000e+00 1.33e+04 7.30e+12  -1.0 9.11e+04    -  4.22e-01 5.03e-04h 11
      44  0.0000000e+00 1.33e+04 1.40e+13  -1.0 9.09e+04    -  1.00e+00 5.03e-04h 11
      45  0.0000000e+00 1.33e+04 1.95e+13  -1.0 9.09e+04    -  4.22e-01 5.03e-04h 11
      46  0.0000000e+00 1.33e+04 3.75e+13  -1.0 9.07e+04    -  1.00e+00 5.03e-04h 11
      47  0.0000000e+00 1.33e+04 5.22e+13  -1.0 9.07e+04    -  4.23e-01 5.03e-04h 11
      48  0.0000000e+00 1.33e+04 1.00e+14  -1.0 9.05e+04    -  1.00e+00 5.03e-04h 11
      49  0.0000000e+00 1.33e+04 1.40e+14  -1.0 9.05e+04    -  4.24e-01 5.03e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      50  0.0000000e+00 1.33e+04 2.68e+14  -1.0 9.03e+04    -  1.00e+00 5.03e-04h 11
      51  0.0000000e+00 1.32e+04 3.73e+14  -1.0 9.03e+04    -  4.24e-01 5.03e-04h 11
      52  0.0000000e+00 7.98e+05 3.49e+14  -1.0 9.01e+04    -  1.00e+00 5.15e-01w  1
      53  0.0000000e+00 3.71e+05 9.23e+14  -1.0 2.25e+04    -  6.55e-01 5.43e-01w  1
      54  0.0000000e+00 2.62e+05 3.29e+15  -1.0 7.71e+03    -  1.19e-01 2.95e-01w  1
      55  0.0000000e+00 1.32e+04 7.18e+14  -1.0 7.82e+03    -  1.00e+00 5.03e-04h 10
      56  0.0000000e+00 1.32e+04 1.00e+15  -1.0 9.01e+04    -  4.25e-01 5.03e-04h 11
      57  0.0000000e+00 1.32e+04 1.13e+15  -1.0 8.99e+04    -  1.00e+00 5.03e-04h 11
      58  0.0000000e+00 1.32e+04 1.13e+15  -1.0 3.94e+04    -  8.00e-01 1.97e-03h  9
      59  0.0000000e+00 1.32e+04 1.13e+15  -1.0 3.15e+04    -  1.00e+00 1.97e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      60  0.0000000e+00 1.31e+04 1.13e+15  -1.0 2.96e+04    -  9.52e-01 1.97e-03h  9
      61  0.0000000e+00 1.31e+04 1.13e+15  -1.0 2.94e+04    -  1.00e+00 1.97e-03h  9
      62  0.0000000e+00 1.31e+04 1.13e+15  -1.0 2.94e+04    -  9.56e-01 1.97e-03h  9
      63  0.0000000e+00 1.31e+04 1.13e+15  -1.0 2.93e+04    -  1.00e+00 1.97e-03h  9
      64  0.0000000e+00 1.30e+04 1.13e+15  -1.0 2.92e+04    -  9.60e-01 1.97e-03h  9
      65  0.0000000e+00 2.85e+05 1.11e+15  -1.0 2.91e+04    -  1.00e+00 5.04e-01w  1
      66  0.0000000e+00 1.32e+05 2.90e+15  -1.0 3.76e+04    -  6.63e-01 5.44e-01w  1
      67  0.0000000e+00 4.04e+05 4.78e+16  -1.0 6.00e+04    -  1.19e-01 5.55e-01w  1
      68  0.0000000e+00 1.30e+04 1.13e+15  -1.0 4.31e+03    -  1.00e+00 1.97e-03h  8
      69  0.0000000e+00 1.30e+04 1.13e+15  -1.0 2.90e+04    -  9.64e-01 1.97e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      70  0.0000000e+00 1.30e+04 1.13e+15  -1.0 2.90e+04    -  1.00e+00 1.97e-03h  9
      71  0.0000000e+00 1.29e+04 1.13e+15  -1.0 2.89e+04    -  9.68e-01 1.97e-03h  9
      72  0.0000000e+00 1.29e+04 1.13e+15  -1.0 2.88e+04    -  1.00e+00 1.97e-03h  9
      73  0.0000000e+00 1.29e+04 1.13e+15  -1.0 2.87e+04    -  9.72e-01 1.97e-03h  9
      74  0.0000000e+00 1.29e+04 1.13e+15  -1.0 2.86e+04    -  1.00e+00 1.97e-03h  9
      75  0.0000000e+00 1.28e+04 1.13e+15  -1.0 2.86e+04    -  9.76e-01 1.97e-03h  9
      76  0.0000000e+00 1.28e+04 1.13e+15  -1.0 2.85e+04    -  1.00e+00 1.97e-03h  9
      77  0.0000000e+00 1.28e+04 1.13e+15  -1.0 2.84e+04    -  9.80e-01 1.97e-03h  9
      78  0.0000000e+00 2.81e+05 1.14e+15  -1.0 2.83e+04    -  1.00e+00 5.04e-01w  1
      79  0.0000000e+00 1.30e+05 3.02e+15  -1.0 3.72e+04    -  6.59e-01 5.45e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      80  0.0000000e+00 3.85e+05 4.73e+16  -1.0 6.09e+04    -  1.17e-01 5.42e-01w  1
      81  0.0000000e+00 1.28e+04 1.13e+15  -1.0 4.31e+03    -  1.00e+00 1.97e-03h  8
      82  0.0000000e+00 1.27e+04 1.13e+15  -1.0 2.83e+04    -  9.84e-01 1.97e-03h  9
      83  0.0000000e+00 1.27e+04 1.13e+15  -1.0 2.82e+04    -  1.00e+00 1.97e-03h  9
      84  0.0000000e+00 1.27e+04 1.13e+15  -1.0 2.81e+04    -  9.88e-01 1.97e-03h  9
      85  0.0000000e+00 1.27e+04 1.13e+15  -1.0 2.80e+04    -  1.00e+00 1.97e-03h  9
      86  0.0000000e+00 1.26e+04 1.13e+15  -1.0 2.79e+04    -  9.92e-01 1.97e-03h  9
      87  0.0000000e+00 1.26e+04 1.13e+15  -1.0 2.79e+04    -  1.00e+00 1.97e-03h  9
      88  0.0000000e+00 1.26e+04 1.13e+15  -1.0 2.78e+04    -  9.96e-01 1.97e-03h  9
      89  0.0000000e+00 1.26e+04 1.13e+15  -1.0 2.77e+04    -  1.00e+00 1.97e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      90  0.0000000e+00 1.25e+04 1.13e+15  -1.0 2.76e+04    -  9.99e-01 1.97e-03h  9
      91  0.0000000e+00 2.69e+05 1.16e+15  -1.0 2.76e+04    -  1.00e+00 5.04e-01w  1
      92  0.0000000e+00 1.24e+05 3.14e+15  -1.0 3.67e+04    -  6.55e-01 5.47e-01w  1
      93  0.0000000e+00 3.72e+05 4.70e+16  -1.0 6.17e+04    -  1.14e-01 5.30e-01w  1
      94  0.0000000e+00 1.25e+04 1.13e+15  -1.0 4.31e+03    -  1.00e+00 1.97e-03h  8
      95  0.0000000e+00 1.25e+04 1.13e+15  -1.0 2.75e+04    -  1.00e+00 1.97e-03h  9
      96  0.0000000e+00 1.25e+04 1.13e+15  -1.0 2.74e+04    -  1.00e+00 1.97e-03h  9
      97  0.0000000e+00 1.24e+04 1.13e+15  -1.0 2.73e+04    -  1.00e+00 1.97e-03h  9
      98  0.0000000e+00 1.24e+04 1.13e+15  -1.0 2.73e+04    -  1.00e+00 1.97e-03h  9
      99  0.0000000e+00 1.24e+04 1.13e+15  -1.0 2.72e+04    -  1.00e+00 1.97e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     100  0.0000000e+00 1.24e+04 1.13e+15  -1.0 2.71e+04    -  1.00e+00 1.97e-03h  9
     101  0.0000000e+00 1.23e+04 1.13e+15  -1.0 2.70e+04    -  1.00e+00 1.97e-03h  9
     102  0.0000000e+00 1.23e+04 1.13e+15  -1.0 2.70e+04    -  1.00e+00 1.97e-03h  9
     103  0.0000000e+00 1.23e+04 1.13e+15  -1.0 2.69e+04    -  1.00e+00 1.97e-03h  9
     104  0.0000000e+00 2.58e+05 1.19e+15  -1.0 2.68e+04    -  1.00e+00 5.03e-01w  1
     105  0.0000000e+00 1.19e+05 3.27e+15  -1.0 3.63e+04    -  6.51e-01 5.48e-01w  1
     106  0.0000000e+00 3.59e+05 4.73e+16  -1.0 6.22e+04    -  1.12e-01 5.20e-01w  1
     107  0.0000000e+00 1.23e+04 1.13e+15  -1.0 4.30e+03    -  1.00e+00 1.97e-03h  8
     108  0.0000000e+00 1.23e+04 1.13e+15  -1.0 2.67e+04    -  1.00e+00 1.97e-03h  9
     109  0.0000000e+00 1.22e+04 1.13e+15  -1.0 2.67e+04    -  1.00e+00 1.97e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     110  0.0000000e+00 1.22e+04 1.13e+15  -1.0 2.66e+04    -  1.00e+00 1.97e-03h  9
     111  0.0000000e+00 1.22e+04 1.13e+15  -1.0 2.65e+04    -  1.00e+00 1.97e-03h  9
     112  0.0000000e+00 1.22e+04 1.13e+15  -1.0 2.65e+04    -  1.00e+00 1.97e-03h  9
     113  0.0000000e+00 1.21e+04 1.13e+15  -1.0 2.64e+04    -  1.00e+00 1.97e-03h  9
     114  0.0000000e+00 1.21e+04 1.14e+15  -1.0 2.63e+04    -  1.00e+00 1.97e-03h  9
     115  0.0000000e+00 1.21e+04 1.14e+15  -1.0 2.62e+04    -  1.00e+00 1.97e-03h  9
     116  0.0000000e+00 1.21e+04 1.14e+15  -1.0 2.62e+04    -  1.00e+00 1.97e-03h  9
     117  0.0000000e+00 2.48e+05 1.21e+15  -1.0 2.61e+04    -  1.00e+00 5.03e-01w  1
     118  0.0000000e+00 1.14e+05 3.41e+15  -1.0 3.59e+04    -  6.47e-01 5.50e-01w  1
     119  0.0000000e+00 3.46e+05 5.02e+16  -1.0 6.12e+04    -  1.10e-01 5.21e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     120  0.0000000e+00 1.20e+04 1.14e+15  -1.0 4.25e+03    -  1.00e+00 1.97e-03h  8
     121  0.0000000e+00 1.20e+04 1.14e+15  -1.0 2.60e+04    -  1.00e+00 1.97e-03h  9
     122  0.0000000e+00 1.20e+04 1.14e+15  -1.0 2.60e+04    -  1.00e+00 1.97e-03h  9
     123  0.0000000e+00 1.20e+04 1.14e+15  -1.0 2.59e+04    -  1.00e+00 1.97e-03h  9
     124  0.0000000e+00 1.19e+04 1.14e+15  -1.0 2.58e+04    -  1.00e+00 1.97e-03h  9
     125  0.0000000e+00 1.19e+04 1.14e+15  -1.0 2.57e+04    -  1.00e+00 1.97e-03h  9
     126  0.0000000e+00 1.19e+04 1.14e+15  -1.0 2.57e+04    -  1.00e+00 1.97e-03h  9
     127  0.0000000e+00 1.19e+04 1.14e+15  -1.0 2.56e+04    -  1.00e+00 1.97e-03h  9
     128  0.0000000e+00 1.19e+04 1.14e+15  -1.0 2.55e+04    -  1.00e+00 1.97e-03h  9
     129  0.0000000e+00 1.18e+04 1.14e+15  -1.0 2.55e+04    -  1.00e+00 1.97e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     130  0.0000000e+00 2.38e+05 1.24e+15  -1.0 2.54e+04    -  1.00e+00 5.03e-01w  1
     131  0.0000000e+00 1.09e+05 3.56e+15  -1.0 3.56e+04    -  6.43e-01 5.51e-01w  1
     132  0.0000000e+00 3.31e+05 6.55e+16  -1.0 5.53e+04    -  1.08e-01 5.63e-01w  1
     133  0.0000000e+00 1.18e+04 1.14e+15  -1.0 4.08e+03    -  1.00e+00 1.97e-03h  8
     134  0.0000000e+00 1.18e+04 1.14e+15  -1.0 2.53e+04    -  1.00e+00 1.97e-03h  9
     135  0.0000000e+00 1.18e+04 1.14e+15  -1.0 2.52e+04    -  1.00e+00 1.97e-03h  9
     136  0.0000000e+00 1.17e+04 1.14e+15  -1.0 2.52e+04    -  1.00e+00 1.97e-03h  9
     137  0.0000000e+00 1.17e+04 1.14e+15  -1.0 2.51e+04    -  1.00e+00 1.97e-03h  9
     138  0.0000000e+00 1.17e+04 1.14e+15  -1.0 2.50e+04    -  1.00e+00 1.97e-03h  9
     139  0.0000000e+00 1.17e+04 1.15e+15  -1.0 2.49e+04    -  1.00e+00 1.97e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     140  0.0000000e+00 1.16e+04 1.15e+15  -1.0 2.49e+04    -  1.00e+00 1.97e-03h  9
     141  0.0000000e+00 1.16e+04 1.15e+15  -1.0 2.48e+04    -  1.00e+00 1.97e-03h  9
     142  0.0000000e+00 1.16e+04 1.15e+15  -1.0 2.47e+04    -  1.00e+00 1.97e-03h  9
     143  0.0000000e+00 2.28e+05 1.27e+15  -1.0 2.46e+04    -  1.00e+00 5.03e-01w  1
     144  0.0000000e+00 1.04e+05 3.72e+15  -1.0 3.53e+04    -  6.37e-01 5.53e-01w  1
     145  0.0000000e+00 3.08e+05 1.85e+17  -1.0 3.85e+04    -  1.07e-01 7.55e-01w  1
     146  0.0000000e+00 1.16e+04 1.15e+15  -1.0 4.52e+03    -  1.00e+00 1.97e-03h  8
     147  0.0000000e+00 1.16e+04 1.15e+15  -1.0 2.46e+04    -  1.00e+00 1.97e-03h  9
     148  0.0000000e+00 1.15e+04 1.15e+15  -1.0 2.45e+04    -  1.00e+00 1.97e-03h  9
     149  0.0000000e+00 1.15e+04 1.15e+15  -1.0 2.44e+04    -  1.00e+00 1.97e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     150  0.0000000e+00 1.15e+04 1.15e+15  -1.0 2.43e+04    -  1.00e+00 1.97e-03h  9
     151  0.0000000e+00 1.15e+04 1.15e+15  -1.0 2.42e+04    -  1.00e+00 1.97e-03h  9
     152  0.0000000e+00 1.14e+04 1.15e+15  -1.0 2.41e+04    -  1.00e+00 1.97e-03h  9
     153  0.0000000e+00 1.14e+04 1.15e+15  -1.0 2.40e+04    -  1.00e+00 1.97e-03h  9
     154  0.0000000e+00 1.14e+04 1.15e+15  -1.0 2.40e+04    -  1.00e+00 1.97e-03h  9
     155  0.0000000e+00 1.14e+04 1.15e+15  -1.0 2.39e+04    -  1.00e+00 1.96e-03h  9
     156  0.0000000e+00 2.17e+05 1.30e+15  -1.0 2.37e+04    -  1.00e+00 5.03e-01w  1
     157  0.0000000e+00 9.99e+04 3.90e+15  -1.0 3.57e+04    -  6.24e-01 5.54e-01w  1
     158  0.0000000e+00 1.31e+05 7.87e+17  -1.0 1.57e+04    -  1.05e-01 9.31e-01w  1
     159  0.0000000e+00 1.14e+04 1.15e+15  -1.0 6.99e+05    -  1.00e+00 1.96e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     160  0.0000000e+00 1.13e+04 1.16e+15  -1.0 2.36e+04    -  1.00e+00 1.96e-03h  9
     161  0.0000000e+00 1.13e+04 1.16e+15  -1.0 2.35e+04    -  1.00e+00 1.96e-03h  9
     162  0.0000000e+00 1.13e+04 1.16e+15  -1.0 2.34e+04    -  1.00e+00 1.96e-03h  9
     163  0.0000000e+00 1.13e+04 1.16e+15  -1.0 2.33e+04    -  1.00e+00 1.96e-03h  9
     164  0.0000000e+00 1.12e+04 1.16e+15  -1.0 2.31e+04    -  1.00e+00 1.96e-03h  9
     165  0.0000000e+00 1.12e+04 1.16e+15  -1.0 2.29e+04    -  1.00e+00 1.96e-03h  9
     166  0.0000000e+00 1.12e+04 1.16e+15  -1.0 2.28e+04    -  1.00e+00 1.96e-03h  9
     167  0.0000000e+00 1.12e+04 1.16e+15  -1.0 2.26e+04    -  1.00e+00 1.96e-03h  9
     168  0.0000000e+00 1.12e+04 1.16e+15  -1.0 2.24e+04    -  1.00e+00 1.96e-03h  9
     169  0.0000000e+00 2.02e+05 1.34e+15  -1.0 2.21e+04    -  1.00e+00 5.03e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     170  0.0000000e+00 9.58e+04 4.16e+15  -1.0 3.80e+04    -  5.87e-01 5.55e-01w  1
     171  0.0000000e+00 8.97e+04 2.18e+18  -1.0 4.09e+03    -  1.03e-01 9.81e-01w  1
     172  0.0000000e+00 1.11e+04 1.16e+15  -1.0 4.36e+06    -  1.00e+00 1.96e-03h  8
     173  0.0000000e+00 1.11e+04 1.16e+15  -1.0 2.18e+04    -  1.00e+00 1.96e-03h  9
     174  0.0000000e+00 1.11e+04 1.16e+15  -1.0 2.15e+04    -  1.00e+00 1.96e-03h  9
     175  0.0000000e+00 1.11e+04 1.17e+15  -1.0 2.11e+04    -  1.00e+00 1.96e-03h  9
     176  0.0000000e+00 1.10e+04 1.17e+15  -1.0 2.07e+04    -  1.00e+00 1.96e-03h  9
     177  0.0000000e+00 1.10e+04 1.17e+15  -1.0 2.02e+04    -  1.00e+00 1.96e-03h  9
     178  0.0000000e+00 1.10e+04 1.17e+15  -1.0 1.96e+04    -  1.00e+00 1.96e-03h  9
     179  0.0000000e+00 1.10e+04 1.17e+15  -1.0 1.89e+04    -  1.00e+00 1.96e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     180  0.0000000e+00 1.10e+04 1.17e+15  -1.0 1.80e+04    -  1.00e+00 1.96e-03h  9
     181  0.0000000e+00 1.09e+04 1.17e+15  -1.0 1.72e+04    -  1.00e+00 3.92e-03h  8
     182  0.0000000e+00 1.59e+05 1.38e+15  -1.0 1.67e+04    -  1.00e+00 5.02e-01w  1
     183  0.0000000e+00 9.33e+04 4.59e+15  -1.0 4.75e+04    -  5.08e-01 5.55e-01w  1
     184  0.0000000e+00 9.75e+04 2.25e+18  -1.0 2.67e+03    -  1.02e-01 9.81e-01w  1
     185  0.0000000e+00 1.09e+04 1.17e+15  -1.0 3.80e+06    -  1.00e+00 3.92e-03h  7
     186  0.0000000e+00 1.08e+04 1.18e+15  -1.0 1.61e+04    -  1.00e+00 3.92e-03h  8
     187  0.0000000e+00 1.08e+04 1.18e+15  -1.0 1.54e+04    -  1.00e+00 3.91e-03h  8
     188  0.0000000e+00 2.03e+04 8.92e+15  -1.0 1.43e+04    -  7.56e-01 7.98e-01H  1
     189  0.0000000e+00 1.74e+04 6.11e+15  -1.0 1.09e+04    -  9.32e-01 1.24e-01h  3
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     190  0.0000000e+00 2.18e+03 3.24e+15  -1.0 7.69e+03    -  9.67e-01 5.03e-01h  1
     191  0.0000000e+00 2.18e+03 1.06e+16  -1.0 2.51e+04    -  2.14e-01 9.54e-04h 11
     192  0.0000000e+00 2.18e+03 1.06e+16  -1.0 2.56e+04    -  2.15e-01 9.54e-04h 11
     193  0.0000000e+00 2.17e+03 1.06e+16  -1.0 2.56e+04    -  2.16e-01 9.54e-04h 11
     194  0.0000000e+00 2.17e+03 1.06e+16  -1.0 2.55e+04    -  2.32e-01 9.54e-04h 11
     195  0.0000000e+00 2.17e+03 1.06e+16  -1.0 2.55e+04    -  1.00e+00 9.54e-04h 11
     196  0.0000000e+00 2.17e+03 1.06e+16  -1.0 2.53e+04    -  2.93e-01 9.54e-04h 11
     197  0.0000000e+00 2.16e+03 1.06e+16  -1.0 2.52e+04    -  1.00e+00 9.54e-04h 11
     198  0.0000000e+00 2.16e+03 1.06e+16  -1.0 2.50e+04    -  2.88e-01 1.91e-03h 10
     199  0.0000000e+00 2.16e+03 1.06e+16  -1.0 2.49e+04    -  1.00e+00 1.91e-03h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     200  0.0000000e+00 2.15e+03 1.05e+16  -1.0 2.47e+04    -  2.96e-01 1.91e-03h 10
     201  0.0000000e+00 2.00e+05 2.31e+17  -1.0 2.46e+04    -  1.00e+00 9.77e-01w  1
     202  0.0000000e+00 1.96e+05 2.28e+17  -1.0 8.75e+04    -  5.63e-02 1.85e-02w  1
     203  0.0000000e+00 5.12e+05 2.22e+17  -1.0 2.37e+05    -  3.03e-02 1.69e-01w  1
     204  0.0000000e+00 2.15e+03 1.05e+16  -1.0 2.44e+05    -  1.00e+00 1.91e-03h  9
     205  0.0000000e+00 2.14e+03 1.05e+16  -1.0 2.44e+04    -  3.04e-01 1.91e-03h 10
     206  0.0000000e+00 2.14e+03 1.05e+16  -1.0 2.43e+04    -  1.00e+00 1.91e-03h 10
     207  0.0000000e+00 2.14e+03 1.05e+16  -1.0 2.41e+04    -  3.12e-01 1.91e-03h 10
     208  0.0000000e+00 2.13e+03 1.04e+16  -1.0 2.40e+04    -  1.00e+00 1.91e-03h 10
     209  0.0000000e+00 2.13e+03 1.04e+16  -1.0 2.38e+04    -  3.20e-01 1.91e-03h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     210  0.0000000e+00 2.12e+03 1.04e+16  -1.0 2.37e+04    -  1.00e+00 1.91e-03h 10
     211  0.0000000e+00 2.12e+03 1.04e+16  -1.0 2.35e+04    -  3.28e-01 1.91e-03h 10
     212  0.0000000e+00 2.12e+03 1.04e+16  -1.0 2.34e+04    -  1.00e+00 1.91e-03h 10
     213  0.0000000e+00 2.11e+03 1.03e+16  -1.0 2.32e+04    -  3.37e-01 1.91e-03h 10
     214  0.0000000e+00 1.78e+05 2.36e+17  -1.0 2.31e+04    -  1.00e+00 9.78e-01w  1
     215  0.0000000e+00 1.76e+05 2.33e+17  -1.0 9.51e+04    -  4.85e-02 1.67e-02w  1
     216  0.0000000e+00 4.92e+05 2.27e+17  -1.0 2.41e+05    -  2.84e-02 1.65e-01w  1
     217  0.0000000e+00 2.11e+03 1.03e+16  -1.0 2.49e+05    -  1.00e+00 1.91e-03h  9
     218  0.0000000e+00 2.10e+03 1.03e+16  -1.0 2.29e+04    -  3.45e-01 1.91e-03h 10
     219  0.0000000e+00 2.10e+03 1.03e+16  -1.0 2.28e+04    -  1.00e+00 1.91e-03h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     220  0.0000000e+00 2.10e+03 1.03e+16  -1.0 2.26e+04    -  3.54e-01 1.91e-03h 10
     221  0.0000000e+00 2.09e+03 1.02e+16  -1.0 2.25e+04    -  1.00e+00 1.91e-03h 10
     222  0.0000000e+00 2.09e+03 1.02e+16  -1.0 2.23e+04    -  3.62e-01 1.91e-03h 10
     223  0.0000000e+00 2.08e+03 1.02e+16  -1.0 2.22e+04    -  1.00e+00 1.91e-03h 10
     224  0.0000000e+00 2.08e+03 1.02e+16  -1.0 2.20e+04    -  3.71e-01 1.91e-03h 10
     225  0.0000000e+00 2.08e+03 1.02e+16  -1.0 2.20e+04    -  1.00e+00 1.91e-03h 10
     226  0.0000000e+00 2.07e+03 1.01e+16  -1.0 2.18e+04    -  3.81e-01 1.91e-03h 10
     227  0.0000000e+00 1.59e+05 2.41e+17  -1.0 2.17e+04    -  1.00e+00 9.79e-01w  1
     228  0.0000000e+00 1.57e+05 2.38e+17  -1.0 1.03e+05    -  4.21e-02 1.50e-02w  1
     229  0.0000000e+00 4.74e+05 2.33e+17  -1.0 2.45e+05    -  2.68e-02 1.60e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     230  0.0000000e+00 2.07e+03 1.01e+16  -1.0 2.54e+05    -  1.00e+00 1.91e-03h  9
     231  0.0000000e+00 2.06e+03 1.01e+16  -1.0 2.15e+04    -  3.90e-01 1.91e-03h 10
     232  0.0000000e+00 2.06e+03 1.01e+16  -1.0 2.14e+04    -  1.00e+00 1.91e-03h 10
     233  0.0000000e+00 2.06e+03 1.01e+16  -1.0 2.12e+04    -  3.99e-01 1.91e-03h 10
     234  0.0000000e+00 2.05e+03 1.00e+16  -1.0 2.11e+04    -  1.00e+00 1.91e-03h 10
     235  0.0000000e+00 2.05e+03 1.00e+16  -1.0 2.09e+04    -  4.09e-01 1.91e-03h 10
     236  0.0000000e+00 2.04e+03 1.00e+16  -1.0 2.09e+04    -  1.00e+00 1.91e-03h 10
     237  0.0000000e+00 2.04e+03 9.98e+15  -1.0 2.07e+04    -  4.19e-01 1.91e-03h 10
     238  0.0000000e+00 2.04e+03 9.96e+15  -1.0 2.06e+04    -  1.00e+00 1.91e-03h 10
     239  0.0000000e+00 2.03e+03 9.94e+15  -1.0 2.04e+04    -  4.29e-01 1.91e-03h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     240  0.0000000e+00 1.42e+05 2.45e+17  -1.0 2.03e+04    -  1.00e+00 9.79e-01w  1
     241  0.0000000e+00 1.41e+05 2.43e+17  -1.0 1.13e+05    -  3.67e-02 1.35e-02w  1
     242  0.0000000e+00 4.57e+05 2.38e+17  -1.0 2.50e+05    -  2.53e-02 1.56e-01w  1
     243  0.0000000e+00 2.03e+03 9.93e+15  -1.0 2.59e+05    -  1.00e+00 1.91e-03h  9
     244  0.0000000e+00 2.02e+03 9.91e+15  -1.0 2.02e+04    -  4.39e-01 1.91e-03h 10
     245  0.0000000e+00 2.02e+03 9.89e+15  -1.0 2.01e+04    -  1.00e+00 1.91e-03h 10
     246  0.0000000e+00 2.02e+03 9.87e+15  -1.0 1.99e+04    -  4.49e-01 1.91e-03h 10
     247  0.0000000e+00 2.01e+03 9.85e+15  -1.0 1.98e+04    -  1.00e+00 1.91e-03h 10
     248  0.0000000e+00 2.01e+03 9.83e+15  -1.0 1.97e+04    -  4.60e-01 1.91e-03h 10
     249  0.0000000e+00 2.00e+03 9.81e+15  -1.0 1.96e+04    -  1.00e+00 1.91e-03h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     250  0.0000000e+00 2.00e+03 9.79e+15  -1.0 1.94e+04    -  4.70e-01 1.91e-03h 10
     251  0.0000000e+00 2.00e+03 9.77e+15  -1.0 1.93e+04    -  1.00e+00 1.91e-03h 10
     252  0.0000000e+00 1.99e+03 9.76e+15  -1.0 1.92e+04    -  4.81e-01 1.91e-03h 10
     253  0.0000000e+00 1.27e+05 2.50e+17  -1.0 1.91e+04    -  1.00e+00 9.80e-01w  1
     254  0.0000000e+00 1.26e+05 2.48e+17  -1.0 1.22e+05    -  3.21e-02 1.22e-02w  1
     255  0.0000000e+00 4.42e+05 2.43e+17  -1.0 2.54e+05    -  2.40e-02 1.52e-01w  1
     256  0.0000000e+00 1.99e+03 9.74e+15  -1.0 2.64e+05    -  1.00e+00 1.91e-03h  9
     257  0.0000000e+00 1.99e+03 9.72e+15  -1.0 1.89e+04    -  4.92e-01 1.91e-03h 10
     258  0.0000000e+00 1.98e+03 9.70e+15  -1.0 1.88e+04    -  1.00e+00 1.91e-03h 10
     259  0.0000000e+00 1.98e+03 9.68e+15  -1.0 1.87e+04    -  5.04e-01 1.91e-03h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     260  0.0000000e+00 1.97e+03 9.66e+15  -1.0 1.86e+04    -  1.00e+00 1.91e-03h 10
     261  0.0000000e+00 1.97e+03 9.64e+15  -1.0 1.84e+04    -  5.15e-01 1.91e-03h 10
     262  0.0000000e+00 1.97e+03 9.63e+15  -1.0 1.84e+04    -  1.00e+00 1.91e-03h 10
     263  0.0000000e+00 1.96e+03 9.61e+15  -1.0 1.82e+04    -  5.27e-01 1.91e-03h 10
     264  0.0000000e+00 1.96e+03 9.59e+15  -1.0 1.81e+04    -  1.00e+00 1.92e-03h 10
     265  0.0000000e+00 3.55e+03 9.39e+15  -1.0 1.80e+04    -  3.06e-01 1.97e-02H  1
     266  0.0000000e+00 3.31e+03 8.76e+15  -1.0 1.41e+04    -  2.72e-03 6.14e-02f  5
     267  0.0000000e+00 3.21e+03 8.48e+15  -1.0 1.38e+04    -  1.00e+00 3.07e-02f  6
     268  0.0000000e+00 3.10e+03 8.21e+15  -1.0 1.37e+04    -  1.25e-01 3.07e-02h  6
     269  0.0000000e+00 3.05e+03 8.08e+15  -1.0 1.35e+04    -  1.00e+00 1.54e-02h  7
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     270  0.0000000e+00 3.01e+03 7.95e+15  -1.0 1.34e+04    -  1.89e-01 1.54e-02h  7
     271  0.0000000e+00 3.89e+03 7.87e+15  -1.0 1.34e+04    -  3.43e-01 1.03e-02H  1
     272  0.0000000e+00 3.63e+03 7.34e+15  -1.0 1.16e+04    -  1.96e-03 6.15e-02f  5
     273  0.0000000e+00 3.52e+03 7.10e+15  -1.0 1.15e+04    -  1.00e+00 3.08e-02f  6
     274  0.0000000e+00 3.40e+03 6.88e+15  -1.0 1.14e+04    -  1.14e-01 3.08e-02h  6
     275  0.0000000e+00 4.12e+03 6.82e+15  -1.0 1.14e+04    -  4.08e-01 8.24e-03H  1
     276  0.0000000e+00 3.85e+03 6.36e+15  -1.0 1.01e+04    -  1.40e-03 6.16e-02f  5
     277  0.0000000e+00 4.23e+03 6.33e+15  -1.0 1.01e+04    -  1.00e+00 4.23e-03H  1
     278  0.0000000e+00 3.94e+03 5.90e+15  -1.0 8.86e+03    -  6.82e-04 6.16e-02f  5
     279  0.0000000e+00 4.31e+03 5.88e+15  -1.0 8.96e+03    -  4.59e-01 4.14e-03H  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     280  0.0000000e+00 4.18e+03 5.69e+15  -1.0 8.21e+03    -  6.73e-04 3.08e-02f  6
     281  0.0000000e+00 4.04e+03 5.51e+15  -1.0 8.27e+03    -  1.00e+00 3.08e-02f  6
     282  0.0000000e+00 3.91e+03 5.33e+15  -1.0 8.34e+03    -  7.09e-02 3.08e-02f  6
     283  0.0000000e+00 3.85e+03 5.24e+15  -1.0 8.39e+03    -  1.00e+00 1.54e-02h  7
     284  0.0000000e+00 3.79e+03 5.16e+15  -1.0 8.42e+03    -  1.34e-01 1.54e-02h  7
     285  0.0000000e+00 3.54e+03 4.81e+15  -1.0 8.44e+03    -  1.00e+00 6.16e-02h  5
     286  0.0000000e+00 4.57e+03 4.75e+15  -1.0 8.51e+03    -  4.16e-01 1.22e-02H  1
     287  0.0000000e+00 4.58e+03 2.97e+15  -1.0 7.25e+03    -  1.00e+00 1.16e-04H  1
    In iteration 287, 2 Slacks too small, adjusting variable bounds
     288  0.0000000e+00 4.58e+03 7.72e+20  -1.0 6.48e+03    -  4.95e-01 1.14e-06H  1
     289  0.0000000e+00 4.58e+03 6.27e+20  -1.0 5.31e+07    -  1.89e-07 9.97e-08f 13
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     290  0.0000000e+00 4.58e+03 4.98e+20  -1.0 4.16e+07    -  2.88e-07 1.27e-07f 13
     291  0.0000000e+00 4.58e+03 3.99e+20  -1.0 4.18e+07    -  4.14e-07 6.33e-08f 14
     292  0.0000000e+00 4.58e+03 3.99e+20  -1.0 3.89e+07    -  4.76e-07 6.81e-08f 14
     293  0.0000000e+00 4.58e+03 3.99e+20  -1.0 3.74e+07    -  5.44e-07 7.07e-08f 14
     294  0.0000000e+00 8.22e+04 9.13e+21  -1.0 3.68e+07    -  6.19e-07 5.89e-04f  1
     295r 0.0000000e+00 8.22e+04 1.00e+03   1.3 0.00e+00  -2.0 0.00e+00 3.59e-07R  6
     296r 0.0000000e+00 2.45e+04 1.10e+03   1.3 2.75e+04    -  4.44e-04 9.93e-04f  1
     297  0.0000000e+00 7.00e+04 3.83e+03  -1.0 9.26e+03    -  9.90e-01 2.78e-01h  1
     298  0.0000000e+00 4.58e+04 1.31e+04  -1.0 1.70e+04    -  9.90e-01 3.04e-01h  2
     299  0.0000000e+00 4.25e+04 4.87e+04  -1.0 1.39e+04    -  9.94e-01 6.95e-02h  4
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     300  0.0000000e+00 4.16e+04 1.02e+05  -1.0 2.69e+04    -  1.00e+00 2.10e-02h  6
     301  0.0000000e+00 4.11e+04 1.53e+05  -1.0 3.81e+04    -  1.00e+00 1.26e-02h  7
     302  0.0000000e+00 4.08e+04 1.98e+05  -1.0 4.37e+04    -  1.00e+00 7.00e-03h  8
     303  0.0000000e+00 5.32e+04 3.29e+05  -1.0 4.67e+04    -  1.00e+00 3.75e-01H  1
     304  0.0000000e+00 1.62e+04 7.70e+05  -1.0 3.30e+04    -  1.35e-01 6.14e-01H  1
     305  0.0000000e+00 1.61e+04 6.42e+06  -1.0 3.16e+04    -  1.00e+00 7.89e-03h  7
     306  0.0000000e+00 1.56e+04 1.02e+07  -1.0 7.07e+04    -  3.44e-01 3.17e-02h  5
     307  0.0000000e+00 1.51e+04 2.68e+07  -1.0 6.16e+04    -  1.00e+00 3.17e-02h  5
     308  0.0000000e+00 1.46e+04 4.65e+07  -1.0 4.82e+04    -  6.26e-01 3.17e-02h  5
     309  0.0000000e+00 1.41e+04 9.88e+07  -1.0 5.21e+04    -  1.00e+00 3.17e-02h  5
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     310  0.0000000e+00 1.36e+04 1.68e+08  -1.0 5.04e+04    -  7.00e-01 3.17e-02h  5
     311  0.0000000e+00 1.35e+04 3.48e+08  -1.0 5.24e+04    -  1.00e+00 7.93e-03h  7
     312  0.0000000e+00 1.34e+04 5.97e+08  -1.0 5.39e+04    -  7.20e-01 7.93e-03h  7
     313  0.0000000e+00 1.33e+04 1.18e+09  -1.0 6.02e+04    -  1.00e+00 7.94e-03h  7
     314  0.0000000e+00 1.32e+04 1.79e+09  -1.0 6.93e+04    -  5.41e-01 7.95e-03h  7
     315  0.0000000e+00 7.29e+05 1.86e+09  -1.0 8.38e+04    -  1.00e+00 5.10e-01w  1
     316  0.0000000e+00 3.49e+05 3.76e+09  -1.0 8.49e+03    -  1.00e+00 5.26e-01w  1
     317  0.0000000e+00 7.49e+05 4.59e+11  -1.0 5.26e+04    -  8.58e-01 9.93e-01w  1
     318  0.0000000e+00 1.31e+04 3.52e+09  -1.0 2.36e+04    -  1.00e+00 3.98e-03h  7
     319  0.0000000e+00 3.27e+05 2.78e+09  -1.0 1.15e+05    -  3.04e-01 2.56e-01h  2
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     320  0.0000000e+00 3.19e+05 4.46e+09  -1.0 1.17e+05    -  4.79e-01 3.40e-02h  4
     321  0.0000000e+00 3.19e+05 9.97e+09  -1.0 1.38e+05    -  1.00e+00 1.83e-03h  8
     322  0.0000000e+00 3.19e+05 1.59e+10  -1.0 1.42e+05    -  5.81e-01 4.47e-04h 10
     323  0.0000000e+00 3.19e+05 3.11e+10  -1.0 1.35e+05    -  1.00e+00 2.33e-04h 11
     324  0.0000000e+00 3.19e+05 5.12e+10  -1.0 1.16e+05    -  7.00e-01 3.30e-05h 14
     325  0.0000000e+00 3.15e+05 9.04e+10  -1.0 9.56e+04    -  1.00e+00 7.96e-02h  3
     326  0.0000000e+00 3.11e+05 1.39e+11  -1.0 5.52e+04    -  1.00e+00 2.38e-01h  2
     327  0.0000000e+00 2.95e+05 1.73e+11  -1.0 6.03e+04    -  4.21e-01 1.28e-01h  3
     328  0.0000000e+00 1.43e+05 2.04e+11  -1.0 7.88e+01  -2.5 1.00e+00 5.14e-01h  1
     329  0.0000000e+00 4.75e+04 4.58e+11  -1.0 3.83e+01  -3.0 1.00e+00 6.68e-01h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     330  0.0000000e+00 3.86e+04 5.36e+13  -1.0 1.27e+01  -3.4 6.67e-02 9.63e-01h  1
     331  0.0000000e+00 1.80e+04 5.23e+12  -1.0 7.89e+00  -3.9 1.00e+00 1.00e+00h  1
     332  0.0000000e+00 2.08e+04 9.78e+07  -1.0 1.90e+00  -4.4 1.00e+00 1.00e+00H  1
     333  0.0000000e+00 3.36e+04 5.43e+05  -1.0 5.38e+00  -4.9 1.00e+00 1.00e+00H  1
     334  0.0000000e+00 4.64e+04 8.02e+10  -2.5 2.96e+00  -5.3 2.57e-01 1.00e+00H  1
     335  0.0000000e+00 5.09e+04 7.24e+10  -2.5 2.12e+02  -5.8 9.72e-02 1.00e+00H  1
     336  0.0000000e+00 4.59e+04 7.09e+10  -2.5 2.40e+04  -5.4 2.56e-02 1.00e+00H  1
     337  0.0000000e+00 4.59e+04 6.97e+10  -2.5 1.37e+05  -4.1 2.02e-02 3.10e-04h  2
     338  0.0000000e+00 4.58e+04 6.93e+10  -2.5 4.29e+05  -3.6 7.56e-03 8.98e-05h  2
     339  0.0000000e+00 4.58e+04 4.33e+10  -2.5 2.16e+04  -3.2 3.84e-01 5.43e-04h  2
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     340  0.0000000e+00 4.55e+04 3.64e+10  -2.5 6.37e+04  -3.7 1.52e-01 5.68e-03h  4
     341  0.0000000e+00 4.53e+04 2.93e+10  -2.5 6.57e+04  -4.2 1.87e-01 4.82e-03h  4
     342  0.0000000e+00 4.51e+04 2.14e+10  -2.5 6.57e+04  -4.6 2.59e-01 4.33e-03h  4
     343  0.0000000e+00 4.48e+04 1.34e+10  -2.5 6.55e+04  -5.1 3.60e-01 6.24e-03h  4
     344  0.0000000e+00 4.46e+04 9.71e+09  -2.5 6.51e+04  -5.6 2.65e-01 6.10e-03h  7
     345  0.0000000e+00 4.41e+04 3.89e+09  -2.5 6.47e+04  -6.1 5.84e-01 1.00e-02h  7
     346  0.0000000e+00 4.28e+04 2.34e+09  -2.5 6.41e+04  -6.5 4.03e-01 2.91e-02h  6
     347  0.0000000e+00 4.67e+03 8.31e+08  -2.5 6.23e+04  -7.0 6.41e-01 1.00e+00w  1
     348  0.0000000e+00 3.94e+02 1.22e+08  -2.5 2.12e+03  -7.5 8.24e-01 1.00e+00w  1
     349  0.0000000e+00 4.04e+03 3.14e+07  -2.5 2.26e+03  -8.0 7.39e-01 8.01e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     350  0.0000000e+00 4.15e+04 9.19e+08  -2.5 2.84e+03  -8.5 6.41e-01 3.12e-02h  5
     351  0.0000000e+00 4.14e+04 9.13e+08  -2.5 2.66e+04  -2.6 1.21e-02 9.29e-04h  2
     352  0.0000000e+00 1.70e+03 4.06e+11  -2.5 1.18e+04  -3.1 6.48e-02 6.27e-03h  1
     353  0.0000000e+00 1.37e+02 2.17e+10  -2.5 6.97e+01  -3.6 8.87e-01 1.00e+00h  1
     354  0.0000000e+00 1.07e+00 3.99e+08  -2.5 3.26e+00  -4.0 1.00e+00 1.00e+00h  1
     355  0.0000000e+00 2.04e-04 4.94e+04  -2.5 7.34e-01  -4.5 1.00e+00 1.00e+00h  1
     356  0.0000000e+00 1.31e-03 5.18e-01  -2.5 2.17e+00  -5.0 1.00e+00 1.00e+00h  1
     357  0.0000000e+00 2.24e-08 1.74e+01  -7.0 2.22e-04  -5.5 1.00e+00 1.00e+00h  1
     358  0.0000000e+00 2.24e-08 6.99e-10  -7.0 6.28e-04  -6.0 1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 358
    
                                       (scaled)                 (unscaled)
    Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    Dual infeasibility......:   6.9888942673322408e-10    6.9888942673322408e-10
    Constraint violation....:   2.9103830456733704e-11    2.2351741790771484e-08
    Complementarity.........:   9.0909093616086122e-08    9.0909093616086122e-08
    Overall NLP error.......:   6.6368354273499748e-09    9.0909093616086122e-08
    
    
    Number of objective function evaluations             = 3103
    Number of objective gradient evaluations             = 359
    Number of equality constraint evaluations            = 3103
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 360
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 358
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.283
    Total CPU secs in NLP function evaluations           =      0.066
    
    EXIT: Optimal Solution Found.


.. code:: ipython3

    # For testing purposes
    # from pyomo.environ import TerminationCondition
    # assert results.solver.termination_condition == TerminationCondition.optimal

Analyze the results of the square problem
-----------------------------------------

What is the total operating cost?

.. code:: ipython3

    print('operating cost = $', value(m.fs.operating_cost))


.. parsed-literal::

    operating cost = $ 364461.48384096316


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
              Heat Duty :      2263.9 : False : (None, None)
        Pressure Change : -2.0000e+05 :  True : (None, None)
    
    ------------------------------------------------------------------------------------
        Stream Table
                                                   Inlet    Vapor Outlet  Liquid Outlet
        flow_mol_phase_comp ('Liq', 'benzene')     0.10734   1.0000e-08       0.10734  
        flow_mol_phase_comp ('Liq', 'toluene')     0.17017   1.0000e-08       0.17017  
        flow_mol_phase_comp ('Liq', 'hydrogen') 2.7751e-07   1.0000e-08    2.7751e-07  
        flow_mol_phase_comp ('Liq', 'methane')  2.7751e-07   1.0000e-08    2.7751e-07  
        flow_mol_phase_comp ('Vap', 'benzene')  1.0000e-08   4.8301e-08    1.0000e-08  
        flow_mol_phase_comp ('Vap', 'toluene')  1.0000e-08   3.1686e-08    1.0000e-08  
        flow_mol_phase_comp ('Vap', 'hydrogen') 1.0000e-08   1.0000e-08    1.0000e-08  
        flow_mol_phase_comp ('Vap', 'methane')  1.0000e-08   1.0000e-08    1.0000e-08  
        temperature                                 325.00       375.00        375.00  
        pressure                                3.5000e+05   1.5000e+05    1.5000e+05  
    ====================================================================================
    
    benzene purity =  0.6038615270937168


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
          <th>flow_mol_phase_comp ('Liq', 'benzene')</th>
          <td>0.019349</td>
          <td>1.000000e-08</td>
        </tr>
        <tr>
          <th>flow_mol_phase_comp ('Liq', 'toluene')</th>
          <td>0.126917</td>
          <td>1.000000e-08</td>
        </tr>
        <tr>
          <th>flow_mol_phase_comp ('Liq', 'hydrogen')</th>
          <td>0.096911</td>
          <td>1.000000e-08</td>
        </tr>
        <tr>
          <th>flow_mol_phase_comp ('Liq', 'methane')</th>
          <td>0.404752</td>
          <td>1.000000e-08</td>
        </tr>
        <tr>
          <th>flow_mol_phase_comp ('Vap', 'benzene')</th>
          <td>0.161092</td>
          <td>7.309994e-02</td>
        </tr>
        <tr>
          <th>flow_mol_phase_comp ('Vap', 'toluene')</th>
          <td>0.082946</td>
          <td>3.969182e-02</td>
        </tr>
        <tr>
          <th>flow_mol_phase_comp ('Vap', 'hydrogen')</th>
          <td>0.793635</td>
          <td>8.905453e-01</td>
        </tr>
        <tr>
          <th>flow_mol_phase_comp ('Vap', 'methane')</th>
          <td>0.305000</td>
          <td>7.097520e-01</td>
        </tr>
        <tr>
          <th>temperature</th>
          <td>767.666119</td>
          <td>3.250000e+02</td>
        </tr>
        <tr>
          <th>pressure</th>
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
