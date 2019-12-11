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


::


    ---------------------------------------------------------------------------

    KeyError                                  Traceback (most recent call last)

    ~/miniconda3/envs/sphinx_gallery/lib/python3.8/site-packages/idaes/dmf/workspace.py in __init__(self, path, create, add_defaults, html_paths)
        209         try:
    --> 210             self._id = self.meta[self.ID_FIELD]
        211         except KeyError:


    KeyError: '_id'

    
    During handling of the above exception, another exception occurred:


    WorkspaceConfMissingField                 Traceback (most recent call last)

    <ipython-input-6-5551fb573031> in <module>
    ----> 1 import hda_ideal_VLE as thermo_props
          2 import hda_reaction as reaction_props


    ~/src/idaes/dangunter/examples-dev/src/properties/Workshop_Module_2_DMF/hda_ideal_VLE.py in <module>
         60 
         61 # Set up DMF (use working directory, wherever that is)
    ---> 62 _dmf = DMF(".")
         63 
         64 


    ~/miniconda3/envs/sphinx_gallery/lib/python3.8/site-packages/idaes/dmf/dmfbase.py in __init__(self, path, name, desc, create, save_path, **ws_kwargs)
        201         ws_kwargs['create'] = create
        202         try:
    --> 203             super(DMF, self).__init__(path, **ws_kwargs)
        204         except (errors.ParseError, ValueError) as err:
        205             msg = 'Configuration "{}", parse error: {}'.format(path, err)


    ~/miniconda3/envs/sphinx_gallery/lib/python3.8/site-packages/idaes/dmf/workspace.py in __init__(self, path, create, add_defaults, html_paths)
        210             self._id = self.meta[self.ID_FIELD]
        211         except KeyError:
    --> 212             raise WorkspaceConfMissingField(path, self.ID_FIELD, 'ID field')
        213         if html_paths:
        214             self.set_doc_paths(html_paths)


    WorkspaceConfMissingField: Workspace config at path "/home/dang/src/idaes/dangunter/examples-dev/src/properties/Workshop_Module_2_DMF" missing ID field field "_id"


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
          <td>0</td>
          <td>benzene</td>
          <td>33870.0</td>
        </tr>
        <tr>
          <td>1</td>
          <td>toluene</td>
          <td>38262.0</td>
        </tr>
        <tr>
          <td>2</td>
          <td>hydrogen</td>
          <td>0.0</td>
        </tr>
        <tr>
          <td>3</td>
          <td>methane</td>
          <td>0.0</td>
        </tr>
        <tr>
          <td>4</td>
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

    2019-09-08 11:30:07 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.H101.control_volume.properties_in initialisation
    2019-09-08 11:30:07 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.H101.control_volume.properties_out initialisation
    2019-09-08 11:30:07 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 1 Complete.
    2019-09-08 11:30:07 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 2 Complete.
    2019-09-08 11:30:07 - INFO - idaes.core.unit_model - fs.H101 Initialisation Complete.
    2019-09-08 11:30:07 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.R101.control_volume.properties_in initialisation
    2019-09-08 11:30:07 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.R101.control_volume.properties_out initialisation
    2019-09-08 11:30:07 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 1 Complete.
    2019-09-08 11:30:07 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 2 Complete.
    2019-09-08 11:30:07 - INFO - idaes.core.unit_model - fs.R101 Initialisation Complete.
    2019-09-08 11:30:07 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.F101.control_volume.properties_in initialisation
    2019-09-08 11:30:07 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.F101.control_volume.properties_out initialisation
    2019-09-08 11:30:08 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 1 Complete.
    2019-09-08 11:30:08 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 2 Complete.
    2019-09-08 11:30:08 - INFO - idaes.core.unit_model - fs.F101 Initialisation Complete.
    2019-09-08 11:30:08 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.S101.mixed_state initialisation
    2019-09-08 11:30:08 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.S101.purge_state initialisation
    2019-09-08 11:30:08 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.S101.recycle_state initialisation
    2019-09-08 11:30:08 - INFO - idaes.unit_models.separator - fs.S101 Initialisation Complete.
    2019-09-08 11:30:08 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.F102.control_volume.properties_in initialisation
    2019-09-08 11:30:08 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.F102.control_volume.properties_out initialisation
    2019-09-08 11:30:08 - INFO - idaes.core.unit_model - fs.F102 Initialisation Step 1 Complete.
    2019-09-08 11:30:08 - INFO - idaes.core.unit_model - fs.F102 Initialisation Step 2 Complete.
    2019-09-08 11:30:08 - INFO - idaes.core.unit_model - fs.F102 Initialisation Complete.
    2019-09-08 11:30:08 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.C101.control_volume.properties_in initialisation
    2019-09-08 11:30:08 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.C101.control_volume.properties_out initialisation
    2019-09-08 11:30:08 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 1 Complete.
    2019-09-08 11:30:08 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 2 Complete.
    2019-09-08 11:30:08 - INFO - idaes.core.unit_model - fs.C101 Initialisation Complete.
    2019-09-08 11:30:08 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.toluene_feed_state initialisation
    2019-09-08 11:30:08 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.hydrogen_feed_state initialisation
    2019-09-08 11:30:08 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.vapor_recycle_state initialisation
    2019-09-08 11:30:08 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.mixed_state initialisation
    2019-09-08 11:30:08 - INFO - idaes.unit_models.mixer - fs.M101 Initialisation Complete.
    2019-09-08 11:30:08 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.H101.control_volume.properties_in initialisation
    2019-09-08 11:30:08 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.H101.control_volume.properties_out initialisation
    2019-09-08 11:30:08 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 1 Complete.
    2019-09-08 11:30:08 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 2 Complete.
    2019-09-08 11:30:08 - INFO - idaes.core.unit_model - fs.H101 Initialisation Complete.
    2019-09-08 11:30:08 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.R101.control_volume.properties_in initialisation
    2019-09-08 11:30:08 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.R101.control_volume.properties_out initialisation
    2019-09-08 11:30:08 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 1 Complete.
    2019-09-08 11:30:08 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 2 Complete.
    2019-09-08 11:30:08 - INFO - idaes.core.unit_model - fs.R101 Initialisation Complete.
    2019-09-08 11:30:08 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.F101.control_volume.properties_in initialisation
    2019-09-08 11:30:08 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.F101.control_volume.properties_out initialisation
    2019-09-08 11:30:08 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 1 Complete.
    2019-09-08 11:30:08 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 2 Complete.
    2019-09-08 11:30:08 - INFO - idaes.core.unit_model - fs.F101 Initialisation Complete.
    2019-09-08 11:30:08 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.S101.mixed_state initialisation
    2019-09-08 11:30:08 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.S101.purge_state initialisation
    2019-09-08 11:30:08 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.S101.recycle_state initialisation
    2019-09-08 11:30:08 - INFO - idaes.unit_models.separator - fs.S101 Initialisation Complete.
    2019-09-08 11:30:08 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.C101.control_volume.properties_in initialisation
    2019-09-08 11:30:09 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.C101.control_volume.properties_out initialisation
    2019-09-08 11:30:09 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 1 Complete.
    2019-09-08 11:30:09 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 2 Complete.
    2019-09-08 11:30:09 - INFO - idaes.core.unit_model - fs.C101 Initialisation Complete.
    2019-09-08 11:30:09 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.toluene_feed_state initialisation
    2019-09-08 11:30:09 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.hydrogen_feed_state initialisation
    2019-09-08 11:30:09 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.vapor_recycle_state initialisation
    2019-09-08 11:30:09 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.mixed_state initialisation
    2019-09-08 11:30:09 - INFO - idaes.unit_models.mixer - fs.M101 Initialisation Complete.
    2019-09-08 11:30:09 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.H101.control_volume.properties_in initialisation
    2019-09-08 11:30:09 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.H101.control_volume.properties_out initialisation
    2019-09-08 11:30:09 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 1 Complete.
    2019-09-08 11:30:09 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 2 Complete.
    2019-09-08 11:30:09 - INFO - idaes.core.unit_model - fs.H101 Initialisation Complete.
    2019-09-08 11:30:09 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.R101.control_volume.properties_in initialisation
    2019-09-08 11:30:09 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.R101.control_volume.properties_out initialisation
    2019-09-08 11:30:09 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 1 Complete.
    2019-09-08 11:30:09 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 2 Complete.
    2019-09-08 11:30:09 - INFO - idaes.core.unit_model - fs.R101 Initialisation Complete.
    2019-09-08 11:30:09 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.F101.control_volume.properties_in initialisation
    2019-09-08 11:30:09 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.F101.control_volume.properties_out initialisation
    2019-09-08 11:30:09 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 1 Complete.
    2019-09-08 11:30:09 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 2 Complete.
    2019-09-08 11:30:09 - INFO - idaes.core.unit_model - fs.F101 Initialisation Complete.
    2019-09-08 11:30:09 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.S101.mixed_state initialisation
    2019-09-08 11:30:09 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.S101.purge_state initialisation
    2019-09-08 11:30:09 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.S101.recycle_state initialisation
    2019-09-08 11:30:09 - INFO - idaes.unit_models.separator - fs.S101 Initialisation Complete.
    2019-09-08 11:30:09 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.C101.control_volume.properties_in initialisation
    2019-09-08 11:30:09 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.C101.control_volume.properties_out initialisation
    2019-09-08 11:30:09 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 1 Complete.
    2019-09-08 11:30:09 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 2 Complete.
    2019-09-08 11:30:09 - INFO - idaes.core.unit_model - fs.C101 Initialisation Complete.
    2019-09-08 11:30:09 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.toluene_feed_state initialisation
    2019-09-08 11:30:09 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.hydrogen_feed_state initialisation
    2019-09-08 11:30:09 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.vapor_recycle_state initialisation
    2019-09-08 11:30:09 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.mixed_state initialisation
    2019-09-08 11:30:09 - INFO - idaes.unit_models.mixer - fs.M101 Initialisation Complete.
    2019-09-08 11:30:09 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.H101.control_volume.properties_in initialisation
    2019-09-08 11:30:09 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.H101.control_volume.properties_out initialisation
    2019-09-08 11:30:09 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 1 Complete.
    2019-09-08 11:30:10 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 2 Complete.
    2019-09-08 11:30:10 - INFO - idaes.core.unit_model - fs.H101 Initialisation Complete.
    2019-09-08 11:30:10 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.R101.control_volume.properties_in initialisation
    2019-09-08 11:30:10 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.R101.control_volume.properties_out initialisation
    2019-09-08 11:30:10 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 1 Complete.
    2019-09-08 11:30:10 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 2 Complete.
    2019-09-08 11:30:10 - INFO - idaes.core.unit_model - fs.R101 Initialisation Complete.
    2019-09-08 11:30:10 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.F101.control_volume.properties_in initialisation
    2019-09-08 11:30:10 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.F101.control_volume.properties_out initialisation
    2019-09-08 11:30:10 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 1 Complete.
    2019-09-08 11:30:10 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 2 Complete.
    2019-09-08 11:30:10 - INFO - idaes.core.unit_model - fs.F101 Initialisation Complete.
    2019-09-08 11:30:10 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.S101.mixed_state initialisation
    2019-09-08 11:30:10 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.S101.purge_state initialisation
    2019-09-08 11:30:10 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.S101.recycle_state initialisation
    2019-09-08 11:30:10 - INFO - idaes.unit_models.separator - fs.S101 Initialisation Complete.
    2019-09-08 11:30:10 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.C101.control_volume.properties_in initialisation
    2019-09-08 11:30:10 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.C101.control_volume.properties_out initialisation
    2019-09-08 11:30:10 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 1 Complete.
    2019-09-08 11:30:10 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 2 Complete.
    2019-09-08 11:30:10 - INFO - idaes.core.unit_model - fs.C101 Initialisation Complete.
    2019-09-08 11:30:10 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.toluene_feed_state initialisation
    2019-09-08 11:30:10 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.hydrogen_feed_state initialisation
    2019-09-08 11:30:10 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.vapor_recycle_state initialisation
    2019-09-08 11:30:10 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.mixed_state initialisation
    2019-09-08 11:30:10 - INFO - idaes.unit_models.mixer - fs.M101 Initialisation Complete.
    2019-09-08 11:30:10 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.H101.control_volume.properties_in initialisation
    2019-09-08 11:30:10 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.H101.control_volume.properties_out initialisation
    2019-09-08 11:30:10 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 1 Complete.
    2019-09-08 11:30:10 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 2 Complete.
    2019-09-08 11:30:10 - INFO - idaes.core.unit_model - fs.H101 Initialisation Complete.
    2019-09-08 11:30:10 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.R101.control_volume.properties_in initialisation
    2019-09-08 11:30:10 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.R101.control_volume.properties_out initialisation
    2019-09-08 11:30:10 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 1 Complete.
    2019-09-08 11:30:10 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 2 Complete.
    2019-09-08 11:30:10 - INFO - idaes.core.unit_model - fs.R101 Initialisation Complete.
    2019-09-08 11:30:10 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.F101.control_volume.properties_in initialisation
    2019-09-08 11:30:10 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.F101.control_volume.properties_out initialisation
    2019-09-08 11:30:10 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 1 Complete.
    2019-09-08 11:30:10 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 2 Complete.
    2019-09-08 11:30:10 - INFO - idaes.core.unit_model - fs.F101 Initialisation Complete.
    2019-09-08 11:30:10 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.S101.mixed_state initialisation
    2019-09-08 11:30:10 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.S101.purge_state initialisation
    2019-09-08 11:30:10 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.S101.recycle_state initialisation
    2019-09-08 11:30:10 - INFO - idaes.unit_models.separator - fs.S101 Initialisation Complete.
    2019-09-08 11:30:11 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.C101.control_volume.properties_in initialisation
    2019-09-08 11:30:11 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.C101.control_volume.properties_out initialisation
    2019-09-08 11:30:11 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 1 Complete.
    2019-09-08 11:30:11 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 2 Complete.
    2019-09-08 11:30:11 - INFO - idaes.core.unit_model - fs.C101 Initialisation Complete.
    2019-09-08 11:30:11 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.toluene_feed_state initialisation
    2019-09-08 11:30:11 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.hydrogen_feed_state initialisation
    2019-09-08 11:30:11 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.vapor_recycle_state initialisation
    2019-09-08 11:30:11 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.mixed_state initialisation
    2019-09-08 11:30:11 - INFO - idaes.unit_models.mixer - fs.M101 Initialisation Complete.
    2019-09-08 11:30:11 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.H101.control_volume.properties_in initialisation
    2019-09-08 11:30:11 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.H101.control_volume.properties_out initialisation
    2019-09-08 11:30:11 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 1 Complete.
    2019-09-08 11:30:11 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 2 Complete.
    2019-09-08 11:30:11 - INFO - idaes.core.unit_model - fs.H101 Initialisation Complete.
    2019-09-08 11:30:11 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.R101.control_volume.properties_in initialisation
    2019-09-08 11:30:11 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.R101.control_volume.properties_out initialisation
    2019-09-08 11:30:11 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 1 Complete.
    2019-09-08 11:30:11 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 2 Complete.
    2019-09-08 11:30:11 - INFO - idaes.core.unit_model - fs.R101 Initialisation Complete.
    2019-09-08 11:30:11 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.F101.control_volume.properties_in initialisation
    2019-09-08 11:30:11 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.F101.control_volume.properties_out initialisation
    2019-09-08 11:30:11 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 1 Complete.
    2019-09-08 11:30:11 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 2 Complete.
    2019-09-08 11:30:11 - INFO - idaes.core.unit_model - fs.F101 Initialisation Complete.
    2019-09-08 11:30:11 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.S101.mixed_state initialisation
    2019-09-08 11:30:11 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.S101.purge_state initialisation
    2019-09-08 11:30:11 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.S101.recycle_state initialisation
    2019-09-08 11:30:11 - INFO - idaes.unit_models.separator - fs.S101 Initialisation Complete.
    2019-09-08 11:30:11 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.C101.control_volume.properties_in initialisation
    2019-09-08 11:30:11 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.C101.control_volume.properties_out initialisation
    2019-09-08 11:30:11 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 1 Complete.
    2019-09-08 11:30:11 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 2 Complete.
    2019-09-08 11:30:11 - INFO - idaes.core.unit_model - fs.C101 Initialisation Complete.
    2019-09-08 11:30:11 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.toluene_feed_state initialisation
    2019-09-08 11:30:11 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.hydrogen_feed_state initialisation
    2019-09-08 11:30:11 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.vapor_recycle_state initialisation
    2019-09-08 11:30:11 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.mixed_state initialisation
    2019-09-08 11:30:11 - INFO - idaes.unit_models.mixer - fs.M101 Initialisation Complete.
    2019-09-08 11:30:11 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.H101.control_volume.properties_in initialisation
    2019-09-08 11:30:11 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.H101.control_volume.properties_out initialisation
    2019-09-08 11:30:11 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 1 Complete.
    2019-09-08 11:30:11 - INFO - idaes.core.unit_model - fs.H101 Initialisation Step 2 Complete.
    2019-09-08 11:30:11 - INFO - idaes.core.unit_model - fs.H101 Initialisation Complete.
    2019-09-08 11:30:12 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.R101.control_volume.properties_in initialisation
    2019-09-08 11:30:12 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.R101.control_volume.properties_out initialisation
    2019-09-08 11:30:12 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 1 Complete.
    2019-09-08 11:30:12 - INFO - idaes.core.unit_model - fs.R101 Initialisation Step 2 Complete.
    2019-09-08 11:30:12 - INFO - idaes.core.unit_model - fs.R101 Initialisation Complete.
    2019-09-08 11:30:12 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.F101.control_volume.properties_in initialisation
    2019-09-08 11:30:12 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.F101.control_volume.properties_out initialisation
    2019-09-08 11:30:12 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 1 Complete.
    2019-09-08 11:30:12 - INFO - idaes.core.unit_model - fs.F101 Initialisation Step 2 Complete.
    2019-09-08 11:30:12 - INFO - idaes.core.unit_model - fs.F101 Initialisation Complete.
    2019-09-08 11:30:12 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.S101.mixed_state initialisation
    2019-09-08 11:30:12 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.S101.purge_state initialisation
    2019-09-08 11:30:12 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.S101.recycle_state initialisation
    2019-09-08 11:30:12 - INFO - idaes.unit_models.separator - fs.S101 Initialisation Complete.
    2019-09-08 11:30:12 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.C101.control_volume.properties_in initialisation
    2019-09-08 11:30:12 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.C101.control_volume.properties_out initialisation
    2019-09-08 11:30:12 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 1 Complete.
    2019-09-08 11:30:12 - INFO - idaes.core.unit_model - fs.C101 Initialisation Step 2 Complete.
    2019-09-08 11:30:12 - INFO - idaes.core.unit_model - fs.C101 Initialisation Complete.
    2019-09-08 11:30:12 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.toluene_feed_state initialisation
    2019-09-08 11:30:12 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.hydrogen_feed_state initialisation
    2019-09-08 11:30:12 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.vapor_recycle_state initialisation
    2019-09-08 11:30:12 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.M101.mixed_state initialisation
    2019-09-08 11:30:12 - INFO - idaes.unit_models.mixer - fs.M101 Initialisation Complete.
    WARNING: Wegstein failed to converge in 5 iterations
    2019-09-08 11:30:12 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.F102.control_volume.properties_in initialisation
    2019-09-08 11:30:12 - INFO - idaes.examples.properties.Workshop_Module_2.hda_ideal_VLE - Starting fs.F102.control_volume.properties_out initialisation
    2019-09-08 11:30:12 - INFO - idaes.core.unit_model - fs.F102 Initialisation Step 1 Complete.
    2019-09-08 11:30:12 - INFO - idaes.core.unit_model - fs.F102 Initialisation Step 2 Complete.
    2019-09-08 11:30:12 - INFO - idaes.core.unit_model - fs.F102 Initialisation Complete.


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
       0  0.0000000e+00 4.23e+04 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  0.0000000e+00 2.49e+04 1.18e+02  -1.0 1.73e+04    -  9.41e-01 1.84e-01h  1
       2  0.0000000e+00 1.11e+04 5.87e+02  -1.0 2.94e+04    -  9.86e-01 4.77e-01h  1
       3  0.0000000e+00 1.11e+04 1.00e+05  -1.0 2.07e+04    -  9.90e-01 3.88e-03h  8
       4  0.0000000e+00 1.11e+04 3.00e+05  -1.0 2.32e+04    -  1.00e+00 3.88e-03h  8
       5  0.0000000e+00 1.10e+04 7.00e+05  -1.0 2.60e+04    -  1.00e+00 1.94e-03h  9
       6  0.0000000e+00 1.10e+04 1.49e+06  -1.0 3.01e+04    -  1.00e+00 1.95e-03h  9
       7  0.0000000e+00 1.10e+04 3.07e+06  -1.0 3.66e+04    -  1.00e+00 1.95e-03h  9
       8  0.0000000e+00 1.10e+04 6.16e+06  -1.0 4.82e+04    -  1.00e+00 1.96e-03h  9
       9  0.0000000e+00 1.10e+04 1.22e+07  -1.0 6.87e+04    -  1.00e+00 9.86e-04h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  0.0000000e+00 1.10e+04 2.37e+07  -1.0 1.02e+05    -  1.00e+00 4.99e-04h 11
      11  0.0000000e+00 1.10e+04 4.54e+07  -1.0 1.43e+05    -  1.00e+00 4.71e-04h 11
      12  0.0000000e+00 1.10e+04 8.61e+07  -1.0 1.65e+05    -  1.00e+00 4.07e-04h 11
      13  0.0000000e+00 1.26e+06 9.62e+07  -1.0 1.68e+05    -  1.00e+00 4.10e-01w  1
      14  0.0000000e+00 6.69e+05 1.84e+08  -1.0 8.28e+03    -  6.72e-01 4.91e-01w  1
      15  0.0000000e+00 3.31e+05 2.72e+08  -1.0 4.82e+04    -  5.43e-01 5.47e-01w  1
      16  0.0000000e+00 1.09e+04 1.63e+08  -1.0 8.17e+03    -  1.00e+00 4.00e-04h 10
      17  0.0000000e+00 1.09e+04 3.09e+08  -1.0 1.63e+05    -  1.00e+00 4.12e-04h 11
      18  0.0000000e+00 1.09e+04 5.86e+08  -1.0 1.55e+05    -  1.00e+00 4.32e-04h 11
      19  0.0000000e+00 1.09e+04 7.93e+08  -1.0 1.47e+05    -  3.90e-01 4.54e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      20  0.0000000e+00 1.09e+04 1.51e+09  -1.0 1.44e+05    -  1.00e+00 4.63e-04h 11
      21  0.0000000e+00 1.09e+04 2.03e+09  -1.0 1.39e+05    -  3.78e-01 4.78e-04h 11
      22  0.0000000e+00 1.09e+04 3.89e+09  -1.0 1.38e+05    -  1.00e+00 4.84e-04h 11
      23  0.0000000e+00 1.09e+04 5.27e+09  -1.0 1.35e+05    -  3.92e-01 4.92e-04h 11
      24  0.0000000e+00 1.09e+04 1.01e+10  -1.0 1.34e+05    -  1.00e+00 4.95e-04h 11
      25  0.0000000e+00 1.09e+04 1.38e+10  -1.0 1.33e+05    -  3.99e-01 4.99e-04h 11
      26  0.0000000e+00 1.39e+06 1.29e+10  -1.0 1.33e+05    -  1.00e+00 5.13e-01w  1
      27  0.0000000e+00 6.91e+05 3.18e+10  -1.0 6.38e+03    -  6.49e-01 5.14e-01w  1
      28  0.0000000e+00 1.80e+05 9.30e+11  -1.0 3.27e+03    -  2.35e-01 7.51e-01w  1
      29  0.0000000e+00 1.09e+04 2.64e+10  -1.0 1.49e+03    -  1.00e+00 5.01e-04h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      30  0.0000000e+00 1.09e+04 3.61e+10  -1.0 1.32e+05    -  4.02e-01 5.02e-04h 11
      31  0.0000000e+00 1.09e+04 6.92e+10  -1.0 1.32e+05    -  1.00e+00 5.03e-04h 11
      32  0.0000000e+00 1.09e+04 9.48e+10  -1.0 1.31e+05    -  4.04e-01 5.04e-04h 11
      33  0.0000000e+00 1.09e+04 1.82e+11  -1.0 1.31e+05    -  1.00e+00 5.04e-04h 11
      34  0.0000000e+00 1.09e+04 2.49e+11  -1.0 1.31e+05    -  4.06e-01 5.04e-04h 11
      35  0.0000000e+00 1.09e+04 4.78e+11  -1.0 1.31e+05    -  1.00e+00 5.04e-04h 11
      36  0.0000000e+00 1.09e+04 6.55e+11  -1.0 1.31e+05    -  4.07e-01 5.04e-04h 11
      37  0.0000000e+00 1.09e+04 1.26e+12  -1.0 1.30e+05    -  1.00e+00 5.04e-04h 11
      38  0.0000000e+00 1.08e+04 1.73e+12  -1.0 1.30e+05    -  4.08e-01 5.04e-04h 11
      39  0.0000000e+00 1.38e+06 1.60e+12  -1.0 1.30e+05    -  1.00e+00 5.17e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      40  0.0000000e+00 6.48e+05 4.25e+12  -1.0 6.18e+03    -  6.48e-01 5.40e-01w  1
      41  0.0000000e+00 3.87e+05 3.24e+13  -1.0 3.57e+03    -  1.24e-01 4.07e-01w  1
      42  0.0000000e+00 1.08e+04 3.31e+12  -1.0 2.46e+03    -  1.00e+00 5.04e-04h 10
      43  0.0000000e+00 1.08e+04 4.54e+12  -1.0 1.30e+05    -  4.08e-01 5.04e-04h 11
      44  0.0000000e+00 1.08e+04 8.71e+12  -1.0 1.30e+05    -  1.00e+00 5.04e-04h 11
      45  0.0000000e+00 1.08e+04 1.20e+13  -1.0 1.30e+05    -  4.09e-01 5.04e-04h 11
      46  0.0000000e+00 1.08e+04 2.30e+13  -1.0 1.30e+05    -  1.00e+00 5.04e-04h 11
      47  0.0000000e+00 1.08e+04 3.16e+13  -1.0 1.30e+05    -  4.10e-01 5.04e-04h 11
      48  0.0000000e+00 1.08e+04 6.06e+13  -1.0 1.29e+05    -  1.00e+00 5.04e-04h 11
      49  0.0000000e+00 1.08e+04 8.34e+13  -1.0 1.29e+05    -  4.10e-01 5.04e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      50  0.0000000e+00 1.08e+04 1.60e+14  -1.0 1.29e+05    -  1.00e+00 5.04e-04h 11
      51  0.0000000e+00 1.08e+04 2.20e+14  -1.0 1.29e+05    -  4.11e-01 5.04e-04h 11
      52  0.0000000e+00 1.35e+06 2.04e+14  -1.0 1.29e+05    -  1.00e+00 5.16e-01w  1
      53  0.0000000e+00 6.35e+05 5.43e+14  -1.0 6.12e+03    -  6.49e-01 5.41e-01w  1
      54  0.0000000e+00 3.82e+05 4.02e+15  -1.0 3.50e+03    -  1.23e-01 4.02e-01w  1
      55  0.0000000e+00 1.08e+04 4.22e+14  -1.0 2.45e+03    -  1.00e+00 5.04e-04h 10
      56  0.0000000e+00 1.08e+04 5.82e+14  -1.0 1.29e+05    -  4.12e-01 5.04e-04h 11
      57  0.0000000e+00 1.08e+04 1.07e+15  -1.0 1.29e+05    -  1.00e+00 5.04e-04h 11
      58  0.0000000e+00 1.08e+04 1.07e+15  -1.0 1.23e+05    -  4.33e-01 5.03e-04h 11
      59  0.0000000e+00 1.08e+04 1.07e+15  -1.0 9.10e+04    -  1.00e+00 9.95e-04h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      60  0.0000000e+00 1.07e+04 1.07e+15  -1.0 6.69e+04    -  9.00e-01 9.86e-04h 10
      61  0.0000000e+00 1.07e+04 1.07e+15  -1.0 6.63e+04    -  1.00e+00 9.86e-04h 10
      62  0.0000000e+00 1.07e+04 1.07e+15  -1.0 6.62e+04    -  9.09e-01 9.86e-04h 10
      63  0.0000000e+00 1.07e+04 1.07e+15  -1.0 6.61e+04    -  1.00e+00 9.86e-04h 10
      64  0.0000000e+00 1.07e+04 1.07e+15  -1.0 6.60e+04    -  9.11e-01 9.86e-04h 10
      65  0.0000000e+00 5.45e+05 1.04e+15  -1.0 6.59e+04    -  1.00e+00 5.05e-01w  1
      66  0.0000000e+00 2.50e+05 2.69e+15  -1.0 1.96e+04    -  6.61e-01 5.41e-01w  1
      67  0.0000000e+00 6.21e+05 4.94e+16  -1.0 6.75e+04    -  1.24e-01 5.85e-01w  1
      68  0.0000000e+00 1.07e+04 1.07e+15  -1.0 3.50e+03    -  1.00e+00 9.86e-04h  9
      69  0.0000000e+00 1.07e+04 1.07e+15  -1.0 6.59e+04    -  9.13e-01 9.86e-04h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      70  0.0000000e+00 1.07e+04 1.07e+15  -1.0 6.58e+04    -  1.00e+00 9.86e-04h 10
      71  0.0000000e+00 1.07e+04 1.07e+15  -1.0 6.57e+04    -  9.15e-01 9.86e-04h 10
      72  0.0000000e+00 1.07e+04 1.07e+15  -1.0 6.56e+04    -  1.00e+00 9.86e-04h 10
      73  0.0000000e+00 1.06e+04 1.07e+15  -1.0 6.56e+04    -  9.17e-01 9.86e-04h 10
      74  0.0000000e+00 1.06e+04 1.07e+15  -1.0 6.55e+04    -  1.00e+00 9.86e-04h 10
      75  0.0000000e+00 1.06e+04 1.07e+15  -1.0 6.54e+04    -  9.19e-01 9.86e-04h 10
      76  0.0000000e+00 1.06e+04 1.07e+15  -1.0 6.53e+04    -  1.00e+00 9.86e-04h 10
      77  0.0000000e+00 1.06e+04 1.07e+15  -1.0 6.53e+04    -  9.21e-01 9.86e-04h 10
      78  0.0000000e+00 5.92e+05 1.05e+15  -1.0 6.52e+04    -  1.00e+00 5.05e-01w  1
      79  0.0000000e+00 2.72e+05 2.75e+15  -1.0 1.96e+04    -  6.59e-01 5.42e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      80  0.0000000e+00 5.82e+05 4.90e+16  -1.0 6.78e+04    -  1.23e-01 5.78e-01w  1
      81  0.0000000e+00 1.06e+04 1.07e+15  -1.0 3.53e+03    -  1.00e+00 9.86e-04h  9
      82  0.0000000e+00 1.06e+04 1.07e+15  -1.0 6.51e+04    -  9.23e-01 9.86e-04h 10
      83  0.0000000e+00 1.06e+04 1.07e+15  -1.0 6.51e+04    -  1.00e+00 9.85e-04h 10
      84  0.0000000e+00 1.06e+04 1.07e+15  -1.0 6.50e+04    -  9.25e-01 9.85e-04h 10
      85  0.0000000e+00 1.06e+04 1.07e+15  -1.0 6.49e+04    -  1.00e+00 9.85e-04h 10
      86  0.0000000e+00 1.05e+04 1.07e+15  -1.0 6.49e+04    -  9.27e-01 9.85e-04h 10
      87  0.0000000e+00 1.05e+04 1.07e+15  -1.0 6.48e+04    -  1.00e+00 9.85e-04h 10
      88  0.0000000e+00 1.05e+04 1.07e+15  -1.0 6.47e+04    -  9.29e-01 9.85e-04h 10
      89  0.0000000e+00 1.05e+04 1.07e+15  -1.0 6.46e+04    -  1.00e+00 9.85e-04h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      90  0.0000000e+00 1.05e+04 1.07e+15  -1.0 6.46e+04    -  9.31e-01 9.85e-04h 10
      91  0.0000000e+00 5.81e+05 1.07e+15  -1.0 6.45e+04    -  1.00e+00 5.05e-01w  1
      92  0.0000000e+00 2.66e+05 2.80e+15  -1.0 1.95e+04    -  6.58e-01 5.42e-01w  1
      93  0.0000000e+00 5.69e+05 4.87e+16  -1.0 6.81e+04    -  1.22e-01 5.71e-01w  1
      94  0.0000000e+00 1.05e+04 1.07e+15  -1.0 3.56e+03    -  1.00e+00 9.85e-04h  9
      95  0.0000000e+00 1.05e+04 1.07e+15  -1.0 6.44e+04    -  9.33e-01 9.85e-04h 10
      96  0.0000000e+00 1.05e+04 1.07e+15  -1.0 6.43e+04    -  1.00e+00 9.85e-04h 10
      97  0.0000000e+00 1.05e+04 1.07e+15  -1.0 6.43e+04    -  9.35e-01 9.85e-04h 10
      98  0.0000000e+00 1.04e+04 1.07e+15  -1.0 6.42e+04    -  1.00e+00 9.85e-04h 10
      99  0.0000000e+00 1.04e+04 1.07e+15  -1.0 6.41e+04    -  9.37e-01 9.85e-04h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     100  0.0000000e+00 1.04e+04 1.07e+15  -1.0 6.41e+04    -  1.00e+00 9.85e-04h 10
     101  0.0000000e+00 1.04e+04 1.07e+15  -1.0 6.40e+04    -  9.39e-01 9.85e-04h 10
     102  0.0000000e+00 1.04e+04 1.07e+15  -1.0 6.39e+04    -  1.00e+00 9.85e-04h 10
     103  0.0000000e+00 1.04e+04 1.07e+15  -1.0 6.39e+04    -  9.41e-01 9.85e-04h 10
     104  0.0000000e+00 5.68e+05 1.08e+15  -1.0 6.38e+04    -  1.00e+00 5.04e-01w  1
     105  0.0000000e+00 2.60e+05 2.86e+15  -1.0 1.95e+04    -  6.56e-01 5.43e-01w  1
     106  0.0000000e+00 5.58e+05 4.86e+16  -1.0 6.83e+04    -  1.21e-01 5.65e-01w  1
     107  0.0000000e+00 1.04e+04 1.07e+15  -1.0 3.58e+03    -  1.00e+00 9.85e-04h  9
     108  0.0000000e+00 1.04e+04 1.07e+15  -1.0 6.37e+04    -  9.43e-01 9.85e-04h 10
     109  0.0000000e+00 1.04e+04 1.07e+15  -1.0 6.36e+04    -  1.00e+00 9.85e-04h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     110  0.0000000e+00 1.04e+04 1.07e+15  -1.0 6.36e+04    -  9.45e-01 9.85e-04h 10
     111  0.0000000e+00 1.03e+04 1.07e+15  -1.0 6.35e+04    -  1.00e+00 9.85e-04h 10
     112  0.0000000e+00 1.03e+04 1.07e+15  -1.0 6.34e+04    -  9.47e-01 9.85e-04h 10
     113  0.0000000e+00 1.03e+04 1.07e+15  -1.0 6.34e+04    -  1.00e+00 9.85e-04h 10
     114  0.0000000e+00 1.03e+04 1.07e+15  -1.0 6.33e+04    -  9.49e-01 9.85e-04h 10
     115  0.0000000e+00 1.03e+04 1.07e+15  -1.0 6.32e+04    -  1.00e+00 9.85e-04h 10
     116  0.0000000e+00 1.03e+04 1.07e+15  -1.0 6.32e+04    -  9.51e-01 9.85e-04h 10
     117  0.0000000e+00 5.56e+05 1.09e+15  -1.0 6.31e+04    -  1.00e+00 5.04e-01w  1
     118  0.0000000e+00 2.54e+05 2.92e+15  -1.0 1.94e+04    -  6.54e-01 5.43e-01w  1
     119  0.0000000e+00 5.48e+05 4.92e+16  -1.0 6.82e+04    -  1.19e-01 5.62e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     120  0.0000000e+00 1.03e+04 1.07e+15  -1.0 3.62e+03    -  1.00e+00 9.85e-04h  9
     121  0.0000000e+00 1.03e+04 1.07e+15  -1.0 6.30e+04    -  9.53e-01 9.85e-04h 10
     122  0.0000000e+00 1.03e+04 1.07e+15  -1.0 6.30e+04    -  1.00e+00 9.85e-04h 10
     123  0.0000000e+00 1.03e+04 1.07e+15  -1.0 6.29e+04    -  9.55e-01 9.85e-04h 10
     124  0.0000000e+00 1.02e+04 1.07e+15  -1.0 6.28e+04    -  1.00e+00 9.85e-04h 10
     125  0.0000000e+00 1.02e+04 1.07e+15  -1.0 6.27e+04    -  9.57e-01 9.85e-04h 10
     126  0.0000000e+00 1.02e+04 1.07e+15  -1.0 6.27e+04    -  1.00e+00 9.85e-04h 10
     127  0.0000000e+00 1.02e+04 1.07e+15  -1.0 6.26e+04    -  9.59e-01 9.85e-04h 10
     128  0.0000000e+00 1.02e+04 1.07e+15  -1.0 6.25e+04    -  1.00e+00 9.85e-04h 10
     129  0.0000000e+00 1.02e+04 1.07e+15  -1.0 6.25e+04    -  9.61e-01 9.85e-04h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     130  0.0000000e+00 5.44e+05 1.10e+15  -1.0 6.24e+04    -  1.00e+00 5.04e-01w  1
     131  0.0000000e+00 2.48e+05 2.98e+15  -1.0 1.94e+04    -  6.52e-01 5.44e-01w  1
     132  0.0000000e+00 5.37e+05 5.19e+16  -1.0 6.71e+04    -  1.18e-01 5.67e-01w  1
     133  0.0000000e+00 1.02e+04 1.07e+15  -1.0 3.66e+03    -  1.00e+00 9.85e-04h  9
     134  0.0000000e+00 1.02e+04 1.07e+15  -1.0 6.23e+04    -  9.63e-01 9.85e-04h 10
     135  0.0000000e+00 1.02e+04 1.07e+15  -1.0 6.23e+04    -  1.00e+00 9.85e-04h 10
     136  0.0000000e+00 1.02e+04 1.07e+15  -1.0 6.22e+04    -  9.65e-01 9.85e-04h 10
     137  0.0000000e+00 1.01e+04 1.07e+15  -1.0 6.21e+04    -  1.00e+00 9.85e-04h 10
     138  0.0000000e+00 1.01e+04 1.07e+15  -1.0 6.21e+04    -  9.67e-01 9.85e-04h 10
     139  0.0000000e+00 1.01e+04 1.07e+15  -1.0 6.20e+04    -  1.00e+00 9.85e-04h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     140  0.0000000e+00 1.01e+04 1.07e+15  -1.0 6.19e+04    -  9.69e-01 9.85e-04h 10
     141  0.0000000e+00 1.01e+04 1.07e+15  -1.0 6.18e+04    -  1.00e+00 9.85e-04h 10
     142  0.0000000e+00 1.01e+04 1.07e+15  -1.0 6.18e+04    -  9.71e-01 9.85e-04h 10
     143  0.0000000e+00 5.32e+05 1.11e+15  -1.0 6.17e+04    -  1.00e+00 5.04e-01w  1
     144  0.0000000e+00 2.43e+05 3.04e+15  -1.0 1.94e+04    -  6.49e-01 5.45e-01w  1
     145  0.0000000e+00 5.25e+05 6.33e+16  -1.0 6.28e+04    -  1.17e-01 6.01e-01w  1
     146  0.0000000e+00 1.01e+04 1.07e+15  -1.0 3.74e+03    -  1.00e+00 9.85e-04h  9
     147  0.0000000e+00 1.01e+04 1.07e+15  -1.0 6.16e+04    -  9.73e-01 9.85e-04h 10
     148  0.0000000e+00 1.01e+04 1.07e+15  -1.0 6.16e+04    -  1.00e+00 9.85e-04h 10
     149  0.0000000e+00 1.01e+04 1.07e+15  -1.0 6.15e+04    -  9.75e-01 9.85e-04h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     150  0.0000000e+00 1.00e+04 1.07e+15  -1.0 6.14e+04    -  1.00e+00 9.85e-04h 10
     151  0.0000000e+00 1.00e+04 1.07e+15  -1.0 6.14e+04    -  9.77e-01 9.85e-04h 10
     152  0.0000000e+00 1.00e+04 1.07e+15  -1.0 6.13e+04    -  1.00e+00 9.85e-04h 10
     153  0.0000000e+00 1.00e+04 1.07e+15  -1.0 6.12e+04    -  9.79e-01 9.85e-04h 10
     154  0.0000000e+00 1.00e+04 1.07e+15  -1.0 6.11e+04    -  1.00e+00 9.85e-04h 10
     155  0.0000000e+00 9.99e+03 1.07e+15  -1.0 6.11e+04    -  9.81e-01 9.85e-04h 10
     156  0.0000000e+00 5.21e+05 1.12e+15  -1.0 6.10e+04    -  1.00e+00 5.04e-01w  1
     157  0.0000000e+00 2.37e+05 3.11e+15  -1.0 1.95e+04    -  6.46e-01 5.45e-01w  1
     158  0.0000000e+00 5.09e+05 1.25e+17  -1.0 5.13e+04    -  1.16e-01 7.29e-01w  1
     159  0.0000000e+00 9.98e+03 1.07e+15  -1.0 3.96e+03    -  1.00e+00 9.85e-04h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     160  0.0000000e+00 9.97e+03 1.07e+15  -1.0 6.09e+04    -  9.84e-01 9.85e-04h 10
     161  0.0000000e+00 9.97e+03 1.07e+15  -1.0 6.09e+04    -  1.00e+00 9.85e-04h 10
     162  0.0000000e+00 9.96e+03 1.07e+15  -1.0 6.08e+04    -  9.86e-01 9.85e-04h 10
     163  0.0000000e+00 9.95e+03 1.07e+15  -1.0 6.07e+04    -  1.00e+00 9.85e-04h 10
     164  0.0000000e+00 9.94e+03 1.07e+15  -1.0 6.06e+04    -  9.88e-01 9.85e-04h 10
     165  0.0000000e+00 9.93e+03 1.08e+15  -1.0 6.05e+04    -  1.00e+00 9.85e-04h 10
     166  0.0000000e+00 9.92e+03 1.08e+15  -1.0 6.05e+04    -  9.91e-01 9.85e-04h 10
     167  0.0000000e+00 9.91e+03 1.08e+15  -1.0 6.04e+04    -  1.00e+00 9.85e-04h 10
     168  0.0000000e+00 9.90e+03 1.08e+15  -1.0 6.03e+04    -  9.93e-01 9.85e-04h 10
     169  0.0000000e+00 5.08e+05 1.14e+15  -1.0 6.02e+04    -  1.00e+00 5.04e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     170  0.0000000e+00 2.31e+05 3.18e+15  -1.0 1.98e+04    -  6.39e-01 5.46e-01w  1
     171  0.0000000e+00 3.47e+05 4.69e+17  -1.0 3.21e+04    -  1.16e-01 9.14e-01w  1
     172  0.0000000e+00 9.89e+03 1.08e+15  -1.0 2.86e+05    -  1.00e+00 9.85e-04h  9
     173  0.0000000e+00 9.88e+03 1.08e+15  -1.0 6.01e+04    -  9.96e-01 9.85e-04h 10
     174  0.0000000e+00 9.87e+03 1.08e+15  -1.0 6.00e+04    -  1.00e+00 9.85e-04h 10
     175  0.0000000e+00 9.86e+03 1.08e+15  -1.0 5.99e+04    -  9.99e-01 9.85e-04h 10
     176  0.0000000e+00 9.85e+03 1.08e+15  -1.0 5.98e+04    -  1.00e+00 9.85e-04h 10
     177  0.0000000e+00 9.84e+03 1.08e+15  -1.0 5.97e+04    -  1.00e+00 9.85e-04h 10
     178  0.0000000e+00 9.83e+03 1.08e+15  -1.0 5.96e+04    -  1.00e+00 9.84e-04h 10
     179  0.0000000e+00 9.82e+03 1.08e+15  -1.0 5.95e+04    -  1.00e+00 9.84e-04h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     180  0.0000000e+00 9.81e+03 1.08e+15  -1.0 5.94e+04    -  1.00e+00 9.84e-04h 10
     181  0.0000000e+00 9.80e+03 1.08e+15  -1.0 5.93e+04    -  1.00e+00 9.84e-04h 10
     182  0.0000000e+00 4.93e+05 1.15e+15  -1.0 5.91e+04    -  1.00e+00 5.04e-01w  1
     183  0.0000000e+00 2.24e+05 3.29e+15  -1.0 2.11e+04    -  6.21e-01 5.46e-01w  1
     184  0.0000000e+00 1.40e+05 1.32e+18  -1.0 1.54e+04    -  1.15e-01 9.76e-01w  1
     185  0.0000000e+00 9.79e+03 1.08e+15  -1.0 2.69e+06    -  1.00e+00 9.84e-04h  9
     186  0.0000000e+00 9.78e+03 1.08e+15  -1.0 5.90e+04    -  1.00e+00 9.84e-04h 10
     187  0.0000000e+00 9.77e+03 1.08e+15  -1.0 5.88e+04    -  1.00e+00 9.84e-04h 10
     188  0.0000000e+00 9.76e+03 1.08e+15  -1.0 5.86e+04    -  1.00e+00 9.84e-04h 10
     189  0.0000000e+00 9.75e+03 1.08e+15  -1.0 5.84e+04    -  1.00e+00 9.84e-04h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     190  0.0000000e+00 9.73e+03 1.08e+15  -1.0 5.82e+04    -  1.00e+00 1.97e-03h  9
     191  0.0000000e+00 9.71e+03 1.08e+15  -1.0 5.78e+04    -  1.00e+00 1.97e-03h  9
     192  0.0000000e+00 9.69e+03 1.08e+15  -1.0 5.75e+04    -  1.00e+00 1.97e-03h  9
     193  0.0000000e+00 9.68e+03 1.08e+15  -1.0 5.71e+04    -  1.00e+00 1.97e-03h  9
     194  0.0000000e+00 9.66e+03 1.08e+15  -1.0 5.67e+04    -  1.00e+00 1.97e-03h  9
     195  0.0000000e+00 4.60e+05 1.17e+15  -1.0 5.63e+04    -  1.00e+00 5.04e-01w  1
     196  0.0000000e+00 2.10e+05 3.50e+15  -1.0 2.56e+04    -  5.72e-01 5.47e-01w  1
     197  0.0000000e+00 7.75e+04 1.48e+18  -1.0 8.15e+03    -  1.13e-01 9.76e-01w  1
     198  0.0000000e+00 9.64e+03 1.08e+15  -1.0 2.74e+06    -  1.00e+00 1.97e-03h  8
     199  0.0000000e+00 9.62e+03 1.08e+15  -1.0 5.57e+04    -  1.00e+00 1.97e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     200  0.0000000e+00 9.60e+03 1.08e+15  -1.0 5.52e+04    -  1.00e+00 1.97e-03h  9
     201  0.0000000e+00 9.58e+03 1.08e+15  -1.0 5.45e+04    -  1.00e+00 1.97e-03h  9
     202  0.0000000e+00 9.56e+03 1.08e+15  -1.0 5.37e+04    -  1.00e+00 1.97e-03h  9
     203  0.0000000e+00 9.54e+03 1.08e+15  -1.0 5.27e+04    -  1.00e+00 1.96e-03h  9
     204  0.0000000e+00 9.53e+03 1.08e+15  -1.0 5.15e+04    -  1.00e+00 1.96e-03h  9
     205  0.0000000e+00 9.51e+03 1.09e+15  -1.0 5.00e+04    -  1.00e+00 1.96e-03h  9
     206  0.0000000e+00 9.49e+03 1.09e+15  -1.0 4.81e+04    -  1.00e+00 1.96e-03h  9
     207  0.0000000e+00 9.47e+03 1.09e+15  -1.0 4.57e+04    -  1.00e+00 1.96e-03h  9
     208  0.0000000e+00 3.32e+05 1.21e+15  -1.0 4.24e+04    -  1.00e+00 5.01e-01w  1
     209  0.0000000e+00 1.68e+05 3.84e+15  -1.0 3.78e+04    -  5.05e-01 5.46e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     210  0.0000000e+00 9.10e+04 2.02e+18  -1.0 2.88e+03  -2.0 1.12e-01 9.81e-01h  1
     211  0.0000000e+00 9.10e+04 2.02e+18  -1.0 2.81e+06  -2.5 1.50e-03 4.32e-04f  2
     212  0.0000000e+00 9.10e+04 2.02e+18  -1.0 3.54e+06  -3.0 1.47e-03 2.20e-04f  3
     213  0.0000000e+00 9.10e+04 2.02e+18  -1.0 4.12e+06  -3.4 2.08e-03 2.30e-04f  3
     214  0.0000000e+00 9.10e+04 2.02e+18  -1.0 4.54e+06  -3.9 2.12e-03 1.30e-04f  4
     215  0.0000000e+00 9.10e+04 2.02e+18  -1.0 4.66e+06  -4.4 2.19e-03 1.41e-04f  4
     216  0.0000000e+00 9.10e+04 2.02e+18  -1.0 4.76e+06  -4.9 2.35e-03 1.55e-04f  4
     217  0.0000000e+00 9.10e+04 2.02e+18  -1.0 4.83e+06  -5.3 1.05e-03 3.41e-04f  3
     218  0.0000000e+00 9.10e+04 2.02e+18  -1.0 4.89e+06  -4.9 7.77e-03 4.29e-04f  3
     219  0.0000000e+00 9.10e+04 2.02e+18  -1.0 4.92e+06  -5.4 7.08e-03 6.75e-05h  6
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     220  0.0000000e+00 9.10e+04 2.02e+18  -1.0 4.68e+06  -5.0 1.23e-02 7.27e-05h  6
     221  0.0000000e+00 9.14e+04 2.01e+18  -1.0 4.92e+06  -5.4 4.90e-03 2.27e-03w  1
     222  0.0000000e+00 9.20e+04 2.01e+18  -1.0 4.73e+06  -5.9 1.48e-02 3.34e-03w  1
     223  0.0000000e+00 9.20e+04 2.01e+18  -1.0 4.86e+06  -6.4 1.60e-02 4.11e-05w  1
     224  0.0000000e+00 9.10e+04 2.02e+18  -1.0 4.82e+05  -6.0 4.90e-03 5.67e-04f  2
     225  0.0000000e+00 9.10e+04 2.01e+18  -1.0 3.07e+06  -5.5 2.01e-02 1.17e-03h  3
     226  0.0000000e+00 9.10e+04 2.01e+18  -1.0 4.91e+06  -6.0 5.86e-03 4.42e-04h  4
     227  0.0000000e+00 9.04e+04 2.00e+18  -1.0 6.56e+05  -5.6 1.15e-01 8.89e-03h  3
     228  0.0000000e+00 9.04e+04 2.00e+18  -1.0 4.86e+06  -6.1 6.24e-03 7.34e-05h  7
     229  0.0000000e+00 9.14e+04 1.99e+18  -1.0 4.85e+06  -6.5 4.75e-03 4.75e-03s 15
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     230  0.0000000e+00 9.13e+04 1.99e+18  -1.0 4.80e+06  -7.0 1.40e-03 1.40e-03s 15
     231r 0.0000000e+00 9.13e+04 1.00e+03   1.6 0.00e+00  -6.6 0.00e+00 0.00e+00R  1
     232r 0.0000000e+00 8.99e+04 1.75e+05   1.6 5.18e+04    -  1.03e-04 3.15e-04f  1
     233  0.0000000e+00 6.51e+04 5.74e+05  -1.0 3.65e+04    -  6.21e-01 1.86e-01h  1
     234  0.0000000e+00 2.69e+04 3.90e+05  -1.0 2.11e+04    -  3.49e-01 3.16e-01h  1
     235  0.0000000e+00 2.60e+04 9.15e+05  -1.0 5.23e+04    -  9.89e-01 3.19e-02h  1
     236  0.0000000e+00 4.51e+04 8.83e+05  -1.0 2.65e+05    -  1.29e-01 8.45e-02h  3
     237  0.0000000e+00 5.31e+04 1.45e+06  -1.0 2.43e+05    -  6.74e-01 7.12e-02h  3
     238  0.0000000e+00 5.23e+04 1.55e+06  -1.0 2.28e+05    -  1.73e-01 2.90e-02h  4
     239  0.0000000e+00 5.20e+04 2.40e+06  -1.0 2.22e+05    -  9.91e-01 2.63e-02h  4
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     240  0.0000000e+00 5.18e+04 2.48e+06  -1.0 2.17e+05    -  1.60e-01 1.19e-02h  5
     241  0.0000000e+00 5.17e+04 4.75e+06  -1.0 2.14e+05    -  1.00e+00 2.82e-03h  7
     242  0.0000000e+00 5.17e+04 6.97e+06  -1.0 2.13e+05    -  1.59e-01 1.40e-03h  8
     243  0.0000000e+00 5.16e+04 3.00e+07  -1.0 2.13e+05    -  1.00e+00 6.94e-04h  9
     244  0.0000000e+00 5.24e+04 4.55e+07  -1.0 2.12e+05    -  1.59e-01 2.22e-02h  4
     245  0.0000000e+00 5.26e+04 2.76e+08  -1.0 2.08e+05    -  1.00e+00 1.99e-02h  4
     246  0.0000000e+00 8.73e+04 2.33e+08  -1.0 2.04e+05    -  1.41e-01 1.43e-01w  1
     247  0.0000000e+00 8.71e+04 5.17e+10  -1.0 1.70e+05    -  1.63e-01 1.60e-03w  1
     248  0.0000000e+00 8.71e+04 4.01e+13  -1.0 1.80e+05    -  1.63e-01 1.84e-04w  1
     249  0.0000000e+00 5.34e+04 4.38e+08  -1.0 5.68e+05    -  1.41e-01 3.58e-02h  2
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     250  0.0000000e+00 5.34e+04 2.30e+09  -1.0 1.98e+05    -  4.78e-01 2.73e-02h  3
     251  0.0000000e+00 5.54e+04 4.26e+09  -1.0 1.93e+05    -  1.22e-01 4.33e-02h  2
     252  0.0000000e+00 5.68e+04 6.04e+10  -1.0 1.86e+05    -  5.11e-01 4.45e-02h  1
     253  0.0000000e+00 5.66e+04 1.37e+12  -1.0 1.81e+05    -  5.97e-02 2.90e-03h  1
     254  0.0000000e+00 5.66e+04 5.54e+15  -1.0 1.14e+05    -  5.47e-01 1.14e-04h  1
     255  0.0000000e+00 5.66e+04 5.42e+13  -1.0 4.14e-06  17.2 9.90e-01 1.00e+00h  1
     256  0.0000000e+00 5.66e+04 8.87e+11  -1.0 2.49e-05  16.8 1.00e+00 2.21e-01h  1
     257r 0.0000000e+00 5.66e+04 1.00e+03   1.1 0.00e+00  16.3 0.00e+00 3.34e-07R 13
     258r 0.0000000e+00 2.64e+04 1.89e+05   1.1 1.83e+04    -  3.84e-03 9.92e-04f  1
     259  0.0000000e+00 2.59e+04 2.98e+04  -1.0 6.87e+04    -  4.83e-01 1.15e-01f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     260  0.0000000e+00 2.58e+04 5.12e+05  -1.0 7.38e+04    -  3.53e-01 4.98e-03h  1
     261  0.0000000e+00 2.68e+04 5.25e+06  -1.0 1.45e+05    -  5.24e-01 5.23e-02h  1
     262  0.0000000e+00 2.68e+04 6.73e+08  -1.0 1.01e+05    -  9.06e-02 7.60e-04h  1
     263r 0.0000000e+00 2.68e+04 1.00e+03   0.5 0.00e+00    -  0.00e+00 3.01e-07R  6
     264r 0.0000000e+00 2.66e+04 1.94e+04   0.5 3.35e+03    -  6.06e-03 9.91e-04f  1
     265  0.0000000e+00 2.58e+04 5.30e+04  -1.0 6.92e+04    -  4.83e-01 3.60e-02f  1
     266  0.0000000e+00 2.57e+04 7.43e+05  -1.0 1.11e+05    -  1.03e-01 4.42e-03h  1
     267  0.0000000e+00 2.56e+04 1.80e+07  -1.0 1.28e+05    -  3.39e-01 1.59e-02h  1
     268  0.0000000e+00 2.56e+04 7.64e+09  -1.0 9.92e+04    -  8.87e-02 2.06e-04h  1
     269r 0.0000000e+00 2.56e+04 1.00e+03  -0.2 0.00e+00    -  0.00e+00 2.59e-07R  5
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     270r 0.0000000e+00 2.48e+04 3.26e+03  -0.2 1.79e+04    -  3.61e-03 9.91e-04f  1
     271  0.0000000e+00 2.46e+04 7.23e+04  -1.0 6.94e+04    -  4.83e-01 9.59e-03f  1
     272  0.0000000e+00 2.45e+04 2.62e+06  -1.0 1.12e+05    -  7.08e-02 3.25e-03h  1
     273  0.0000000e+00 2.45e+04 2.34e+08  -1.0 1.17e+05    -  2.23e-01 2.66e-03h  1
     274  0.0000000e+00 2.45e+04 6.07e+11  -1.0 1.08e+05    -  7.96e-02 2.89e-05h  1
     275  0.0000000e+00 2.45e+04 3.09e+16  -1.0 5.39e+05    -  9.90e-01 1.88e-05h  1
     276  0.0000000e+00 2.44e+04 3.09e+14  -1.0 2.39e-06  15.8 9.90e-01 7.13e-01h  1
     277  0.0000000e+00 2.44e+04 2.50e+12  -1.0 1.32e-06  16.2 9.92e-01 1.28e-02h  1
     278  0.0000000e+00 2.44e+04 1.14e+10  -1.0 6.32e-06  15.8 1.00e+00 4.23e-03h  1
     279  0.0000000e+00 2.43e+04 1.58e+10  -1.0 1.43e-06  16.2 1.00e+00 4.13e-01h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     280  0.0000000e+00 2.41e+04 3.22e+10  -1.0 8.70e-06  15.7 1.00e+00 5.58e-01h  1
     281  0.0000000e+00 2.39e+04 4.61e+10  -1.0 3.48e-06  16.1 1.00e+00 8.91e-01h  1
     282  0.0000000e+00 2.35e+04 9.77e+10  -1.0 2.95e-05  15.7 1.00e+00 5.46e-01h  1
     283  0.0000000e+00 2.31e+04 1.39e+11  -1.0 1.13e-05  16.1 1.00e+00 1.00e+00f  1
     284  0.0000000e+00 2.31e+04 6.60e+11  -1.0 3.66e-04  15.6 1.00e+00 1.25e-02h  4
     285  0.0000000e+00 2.26e+04 2.48e+11  -1.0 2.21e-05  16.0 1.00e+00 1.00e+00s 22
     286r 0.0000000e+00 2.26e+04 1.00e+03  -0.8 0.00e+00  16.5 0.00e+00 0.00e+00R  1
     287r 0.0000000e+00 2.10e+04 4.95e+05  -0.8 3.60e+04    -  2.04e-03 1.01e-03f  1
     288  0.0000000e+00 2.46e+04 2.45e+05  -1.0 6.94e+04    -  4.83e-01 5.71e-02f  1
     289  0.0000000e+00 2.38e+04 1.08e+06  -1.0 5.58e+04    -  3.63e-02 1.51e-01f  2
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     290  0.0000000e+00 2.31e+04 1.34e+06  -1.0 5.27e+04    -  7.92e-01 5.91e-02f  3
     291  0.0000000e+00 1.35e+04 2.11e+06  -1.0 5.66e+03    -  1.68e-01 5.02e-01h  1
     292  0.0000000e+00 1.68e+04 2.89e+08  -1.0 1.29e+04    -  1.57e-01 7.49e-01H  1
     293  0.0000000e+00 1.39e+04 1.64e+09  -1.0 1.68e+04    -  8.12e-02 1.76e-01H  1
     294  0.0000000e+00 1.38e+04 4.25e+11  -1.0 1.51e+04    -  3.19e-01 3.04e-03H  1
    In iteration 294, 1 Slack too small, adjusting variable bound
     295  0.0000000e+00 1.38e+04 8.06e+15  -1.0 1.57e+04    -  6.05e-01 3.23e-05F  1
     296  0.0000000e+00 1.30e+04 1.48e+18  -1.0 3.37e+04    -  1.17e-04 9.55e-02f  2
     297r 0.0000000e+00 1.30e+04 1.00e+03  -1.0 0.00e+00    -  0.00e+00 4.03e-07R  5
     298r 0.0000000e+00 1.26e+04 1.51e+05  -1.0 6.59e+04    -  5.36e-04 5.34e-04f  1
     299  0.0000000e+00 7.03e+03 8.67e+05  -1.0 1.65e+04    -  9.78e-01 5.10e-01h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     300  0.0000000e+00 6.79e+03 6.63e+08  -1.0 4.17e+03    -  9.90e-01 3.36e-02h  5
     301  0.0000000e+00 6.76e+03 1.95e+09  -1.0 3.32e+03    -  9.92e-01 4.21e-03h  8
     302  0.0000000e+00 6.75e+03 4.34e+09  -1.0 4.06e+03    -  1.00e+00 2.11e-03h  9
     303  0.0000000e+00 6.73e+03 8.70e+09  -1.0 3.47e+03    -  1.00e+00 2.11e-03h  9
     304  0.0000000e+00 6.72e+03 1.67e+10  -1.0 3.57e+03    -  1.00e+00 2.11e-03h  9
     305  0.0000000e+00 6.71e+03 3.12e+10  -1.0 2.93e+03    -  1.00e+00 2.11e-03h  9
     306  0.0000000e+00 6.69e+03 5.77e+10  -1.0 2.05e+03    -  1.00e+00 2.11e-03h  9
     307  0.0000000e+00 6.68e+03 1.06e+11  -1.0 2.10e+03    -  1.00e+00 2.12e-03h  9
     308  0.0000000e+00 6.66e+03 1.63e+11  -1.0 4.82e+03    -  6.61e-01 2.13e-03h  9
     309  0.0000000e+00 6.65e+03 2.08e+11  -1.0 1.10e+04    -  3.54e-01 2.15e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     310  0.0000000e+00 9.25e+03 1.50e+11  -1.0 1.38e+04    -  3.02e-01 5.52e-01w  1
     311  0.0000000e+00 7.29e+03 7.08e+12  -1.0 3.04e+03    -  2.54e-01 9.91e-01w  1
     312  0.0000000e+00 7.84e+03 6.96e+12  -1.0 2.07e+05    -  2.90e-02 3.48e-02w  1
     313  0.0000000e+00 6.63e+03 2.58e+11  -1.0 7.02e+03    -  3.02e-01 2.16e-03h  8
     314  0.0000000e+00 6.62e+03 3.17e+11  -1.0 1.62e+04    -  2.98e-01 2.16e-03h  9
     315  0.0000000e+00 6.61e+03 3.87e+11  -1.0 1.76e+04    -  2.81e-01 2.17e-03h  9
     316  0.0000000e+00 6.59e+03 4.73e+11  -1.0 1.86e+04    -  2.84e-01 2.17e-03h  9
     317  0.0000000e+00 6.58e+03 5.75e+11  -1.0 1.92e+04    -  2.79e-01 2.17e-03h  9
     318  0.0000000e+00 6.56e+03 7.02e+11  -1.0 1.96e+04    -  2.83e-01 2.16e-03h  9
     319  0.0000000e+00 6.55e+03 8.56e+11  -1.0 1.98e+04    -  2.83e-01 2.15e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     320  0.0000000e+00 6.53e+03 1.21e+12  -1.0 1.30e+04    -  5.18e-01 2.16e-03h  9
     321  0.0000000e+00 6.52e+03 1.37e+12  -1.0 2.35e+04    -  1.72e-01 1.90e-03h  9
     322  0.0000000e+00 6.51e+03 1.90e+12  -1.0 2.41e+04    -  5.03e-01 1.84e-03h  9
     323  0.0000000e+00 1.81e+04 1.25e+12  -1.0 2.09e+04    -  2.01e-01 5.38e-01w  1
     324  0.0000000e+00 7.70e+03 1.19e+13  -1.0 8.85e+02    -  9.74e-02 9.90e-01w  1
     325  0.0000000e+00 6.09e+03 1.78e+13  -1.0 1.05e+04    -  1.08e-01 4.40e-01w  1
     326  0.0000000e+00 6.50e+03 2.19e+12  -1.0 1.39e+04    -  2.01e-01 2.10e-03h  8
     327  0.0000000e+00 6.48e+03 3.05e+12  -1.0 2.10e+04    -  5.02e-01 2.09e-03h  9
     328  0.0000000e+00 6.47e+03 3.56e+12  -1.0 1.99e+04    -  2.17e-01 2.15e-03h  9
     329  0.0000000e+00 6.45e+03 4.96e+12  -1.0 2.01e+04    -  5.02e-01 2.14e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     330  0.0000000e+00 6.44e+03 5.83e+12  -1.0 1.96e+04    -  2.27e-01 2.16e-03h  9
     331  0.0000000e+00 6.43e+03 8.13e+12  -1.0 1.97e+04    -  5.07e-01 2.16e-03h  9
     332  0.0000000e+00 6.41e+03 9.60e+12  -1.0 1.95e+04    -  2.35e-01 2.17e-03h  9
     333  0.0000000e+00 6.40e+03 1.35e+13  -1.0 1.95e+04    -  5.16e-01 2.17e-03h  9
     334  0.0000000e+00 6.39e+03 1.60e+13  -1.0 1.94e+04    -  2.42e-01 2.17e-03h  9
     335  0.0000000e+00 6.37e+03 2.25e+13  -1.0 1.93e+04    -  5.27e-01 2.17e-03h  9
     336  0.0000000e+00 1.66e+04 1.70e+13  -1.0 1.93e+04    -  2.48e-01 5.57e-01w  1
     337  0.0000000e+00 7.53e+03 7.80e+14  -1.0 7.02e+02    -  8.57e-02 9.91e-01w  1
     338  0.0000000e+00 7.46e+03 7.57e+14  -1.0 5.98e+04    -  9.81e-02 1.12e-02w  1
     339  0.0000000e+00 6.36e+03 2.68e+13  -1.0 4.19e+02    -  2.48e-01 2.18e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     340  0.0000000e+00 6.34e+03 3.81e+13  -1.0 1.92e+04    -  5.39e-01 2.18e-03h  9
     341  0.0000000e+00 6.33e+03 4.55e+13  -1.0 1.92e+04    -  2.54e-01 2.18e-03h  9
     342  0.0000000e+00 6.32e+03 6.51e+13  -1.0 1.91e+04    -  5.52e-01 2.18e-03h  9
     343  0.0000000e+00 6.30e+03 7.81e+13  -1.0 1.91e+04    -  2.60e-01 2.18e-03h  9
     344  0.0000000e+00 6.29e+03 1.13e+14  -1.0 1.90e+04    -  5.66e-01 2.18e-03h  9
     345  0.0000000e+00 6.27e+03 1.35e+14  -1.0 1.90e+04    -  2.65e-01 2.18e-03h  9
     346  0.0000000e+00 6.26e+03 1.97e+14  -1.0 1.89e+04    -  5.80e-01 2.18e-03h  9
     347  0.0000000e+00 6.25e+03 2.37e+14  -1.0 1.89e+04    -  2.71e-01 2.18e-03h  9
     348  0.0000000e+00 6.23e+03 3.47e+14  -1.0 1.88e+04    -  5.94e-01 2.18e-03h  9
     349  0.0000000e+00 1.59e+04 2.54e+14  -1.0 1.88e+04    -  2.76e-01 5.58e-01w  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     350  0.0000000e+00 7.09e+03 1.18e+16  -1.0 6.83e+02    -  1.49e-01 9.91e-01w  1
     351  0.0000000e+00 7.06e+03 1.16e+16  -1.0 4.04e+04    -  1.51e-01 5.51e-03w  1
     352  0.0000000e+00 6.22e+03 4.21e+14  -1.0 4.21e+02    -  2.76e-01 2.18e-03h  8
     353  0.0000000e+00 6.21e+03 6.20e+14  -1.0 1.87e+04    -  6.09e-01 2.18e-03h  9
     354  0.0000000e+00 6.19e+03 7.53e+14  -1.0 1.88e+04    -  2.81e-01 2.17e-03h  9
     355  0.0000000e+00 6.18e+03 1.12e+15  -1.0 1.87e+04    -  6.23e-01 2.17e-03h  9
     356  0.0000000e+00 6.17e+03 1.36e+15  -1.0 1.87e+04    -  2.87e-01 2.16e-03h  9
     357  0.0000000e+00 6.15e+03 2.04e+15  -1.0 1.86e+04    -  6.38e-01 2.17e-03h  9
     358  0.0000000e+00 6.14e+03 2.49e+15  -1.0 1.86e+04    -  2.92e-01 2.16e-03h  9
     359  0.0000000e+00 6.13e+03 3.76e+15  -1.0 1.85e+04    -  6.53e-01 2.16e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     360  0.0000000e+00 6.11e+03 4.61e+15  -1.0 1.85e+04    -  2.97e-01 2.15e-03h  9
     361  0.0000000e+00 6.10e+03 6.99e+15  -1.0 1.84e+04    -  6.68e-01 2.16e-03h  9
     362  0.0000000e+00 1.47e+04 4.57e+15  -1.0 1.84e+04    -  3.02e-01 5.49e-01w  1
     363  0.0000000e+00 6.43e+03 8.73e+16  -1.0 6.80e+02    -  3.24e-01 9.90e-01w  1
     364  0.0000000e+00 6.28e+03 8.01e+16  -1.0 8.13e+03    -  7.01e-01 2.38e-02w  1
     365  0.0000000e+00 6.09e+03 8.61e+15  -1.0 4.14e+02    -  3.02e-01 2.15e-03h  8
     366  0.0000000e+00 6.07e+03 1.32e+16  -1.0 1.83e+04    -  6.83e-01 2.15e-03h  9
     367  0.0000000e+00 6.06e+03 1.62e+16  -1.0 1.83e+04    -  3.07e-01 2.14e-03h  9
     368  0.0000000e+00 6.05e+03 2.50e+16  -1.0 1.82e+04    -  6.98e-01 2.14e-03h  9
     369  0.0000000e+00 6.03e+03 3.10e+16  -1.0 1.82e+04    -  3.12e-01 2.14e-03h  9
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     370  0.0000000e+00 6.02e+03 4.81e+16  -1.0 1.81e+04    -  7.13e-01 2.14e-03h  9
     371  0.0000000e+00 6.01e+03 5.97e+16  -1.0 1.81e+04    -  3.17e-01 2.13e-03h  9
     372  0.0000000e+00 6.00e+03 9.33e+16  -1.0 1.80e+04    -  7.29e-01 2.14e-03h  9
     373  0.0000000e+00 5.98e+03 1.16e+17  -1.0 1.80e+04    -  3.22e-01 2.13e-03h  9
     374  0.0000000e+00 5.97e+03 1.83e+17  -1.0 1.79e+04    -  7.44e-01 2.13e-03h  9
     375  0.0000000e+00 1.36e+04 1.09e+17  -1.0 1.79e+04    -  3.28e-01 5.44e-01w  1
     376  0.0000000e+00 5.90e+03 1.47e+18  -1.0 6.74e+02    -  4.59e-01 9.90e-01w  1
     377  0.0000000e+00 4.63e+03 1.10e+18  -1.0 7.56e+02    -  1.00e+00 2.18e-01w  1
     378  0.0000000e+00 1.76e+02 3.94e+16  -1.0 3.32e+02    -  9.72e-01 1.00e+00h  1
     379  0.0000000e+00 1.76e-01 3.87e+15  -1.0 2.17e+00    -  9.38e-01 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     380  0.0000000e+00 2.36e-06 9.71e+12  -1.0 1.38e-01    -  9.90e-01 1.00e+00f  1
     381  0.0000000e+00 2.19e-02 2.03e+12  -1.0 1.34e+01    -  7.92e-01 1.00e+00f  1
     382  0.0000000e+00 5.23e-01 1.40e+12  -1.0 6.56e+01    -  3.12e-01 1.00e+00f  1
     383  0.0000000e+00 1.16e+00 5.38e+11  -1.0 9.76e+01    -  6.15e-01 1.00e+00f  1
     384  0.0000000e+00 8.83e+00 2.36e+11  -1.0 2.69e+02    -  5.60e-01 1.00e+00f  1
     385  0.0000000e+00 6.02e+01 1.09e+11  -1.0 7.01e+02    -  5.38e-01 1.00e+00f  1
     386  0.0000000e+00 2.51e+00 5.54e+10  -1.0 2.01e+03    -  4.92e-01 1.00e+00F  1
     387  0.0000000e+00 2.24e+02 5.54e+10  -1.0 1.31e+06    -  2.38e-03 1.01e-03f  5
     388  0.0000000e+00 2.68e+02 1.33e+10  -1.0 2.34e+04    -  7.86e-01 2.64e-02h  6
     389  0.0000000e+00 1.72e+02 3.65e+09  -1.0 7.70e+03    -  7.18e-01 1.00e+00H  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     390  0.0000000e+00 3.44e+01 2.12e+06  -1.0 4.58e+03    -  1.00e+00 1.00e+00H  1
     391  0.0000000e+00 8.52e+00 3.96e+07  -1.7 2.36e+02    -  1.00e+00 1.00e+00h  1
     392  0.0000000e+00 1.44e+02 7.44e+03  -1.7 9.71e+02    -  1.00e+00 1.00e+00h  1
     393  0.0000000e+00 1.17e-04 3.09e+04  -2.5 8.28e-01    -  1.00e+00 1.00e+00h  1
     394  0.0000000e+00 6.17e-04 1.03e-02  -2.5 2.21e+00    -  1.00e+00 1.00e+00h  1
     395  0.0000000e+00 6.71e-08 5.81e+00  -7.0 2.47e-06    -  1.00e+00 1.00e+00h  1
     396  0.0000000e+00 8.20e-08 8.81e-10  -7.0 1.90e-04    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 396
    
                                       (scaled)                 (unscaled)
    Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    Dual infeasibility......:   8.8079467308911410e-10    8.8079467308911410e-10
    Constraint violation....:   2.9103830456733704e-11    8.1956386566162109e-08
    Complementarity.........:   9.0909090909091298e-08    9.0909090909091298e-08
    Overall NLP error.......:   7.3586953589023978e-09    9.0909090909091298e-08
    
    
    Number of objective function evaluations             = 3210
    Number of objective gradient evaluations             = 397
    Number of equality constraint evaluations            = 3210
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 403
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 396
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.621
    Total CPU secs in NLP function evaluations           =      0.092
    
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
