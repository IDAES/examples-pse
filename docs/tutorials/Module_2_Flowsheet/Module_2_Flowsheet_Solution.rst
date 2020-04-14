Learning outcomes
-----------------

-  Construct a steady-state flowsheet using the IDAES unit model library
-  Connecting unit models in a flowsheet using Arcs
-  Using the SequentialDecomposition tool to initialize a flowsheet with
   recycle
-  Fomulate and solve an optimization problem

   -  Defining an objective function
   -  Setting variable bounds
   -  Adding additional constraints

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

-  hda\_ideal\_VLE.py
-  hda\_reaction.py

The state variables chosen for the property package are **flows of
component by phase, temperature and pressure**. The components
considered are: **toluene, hydrogen, benzene and methane**. Therefore,
every stream has 8 flow variables, 1 temperature and 1 pressure
variable.

.. figure:: module_2_flowsheet.png
   :alt: 

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
models: - Mixer - Heater - StoichiometricReactor - **Flash** - Separator
(splitter) - PressureChanger

.. code:: ipython3

    from idaes.core import FlowsheetBlock

.. code:: ipython3

    from idaes.generic_models.unit_models import (PressureChanger,
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

    from idaes.generic_models.unit_models import Flash

We will also be needing some utility tools to put together the flowsheet
and calculate the degrees of freedom.

.. code:: ipython3

    from idaes.generic_models.unit_models.pressure_changer import ThermodynamicAssumption
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

::

         <li>hda_ideal_VLE as thermo_props</li>
         <li>hda_reaction as reaction_props </li>
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

Adding Unit Models
------------------

Let us start adding the unit models we have imported to the flowsheet.
Here, we are adding the Mixer (assigned a name M101) and a Heater
(assigned a name H101). Note that, all unit models need to be given a
property package argument. In addition to that, there are several
arguments depending on the unit model, please refer to the documentation
for more details
(https://idaes-pse.readthedocs.io/en/latest/model\_libraries/core\_lib/unit\_models/index.html).
For example, the Mixer unit model here is given a ``list`` consisting of
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

::

         <li>"property_package": m.fs.thermo_params</li>
         <li>"reaction_package": m.fs.reaction_params </li>
         <li>"has_heat_of_reaction": True </li>
         <li>"has_heat_transfer": True</li>
         <li>"has_pressure_change": False</li>
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

::

         <li>"property_package": m.fs.thermo_params</li>
         <li>"has_heat_transfer": True</li>
         <li>"has_pressure_change": False</li>
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

.. figure:: module_2_flowsheet.png
   :alt: 

.. raw:: html

   <div class="alert alert-block alert-info">

Inline Exercise: Now, connect the H101 outlet to the R101 inlet using
the cell above as a guide.

.. raw:: html

   </div>

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
https://pyomo.readthedocs.io/en/latest/pyomo\_modeling\_components/Expressions.html

For this flowsheet, we are interested in computing the purity of the
product Benzene stream (i.e. the mole fraction) and the operating cost
which is a sum of the cooling and heating cost.

Let us first add an Expression to compute the mole fraction of benzene
in the ``vap_outlet`` of F102 which is our product stream. Please note
that the var flow\_mol\_phase\_comp has the index - [time, phase,
component]. As this is a steady-state flowsheet, the time index by
default is 0. The valid phases are ["Liq", "Vap"]. Similarly the valid
component list is ["benzene", "toluene", "hydrogen", "methane"].

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

::

         <li>2.2E-4 dollars/kW for H101</li>
         <li>1.9E-4 dollars/kW for F102</li>
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

    29


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

::

         <li>F<sub>H2</sub> = 0.30 mol/s</li>
         <li>F<sub>CH4</sub> = 0.02 mol/s</li>
         <li>Remaining components = 1e-5 mol/s</li>
         <li>T = 303.2 K</li>
         <li>P = 350000 Pa</li>
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
adiabatic i.e. Q = 0.

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

::

         <li>T = 375 K</li>
         <li>deltaP = -200000</li>
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
degrees of freedom i.e. should be a square problem. Please check that
the degrees of freedom is 0.

Use Shift+Enter to run the cell once you have typed in your code.

.. raw:: html

   </div>

.. code:: ipython3

    print(degrees_of_freedom(m))


.. parsed-literal::

    0


Initialization
--------------

This section will demonstrate how to use the built-in sequential
decomposition tool to initialize our flowsheet.

.. figure:: module_2_flowsheet.png
   :alt: 

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


.. figure:: module_2_tear_stream.png
   :alt: 

The SequentialDecomposition tool has determined that the tear stream is
the mixer outlet. We will need to provide a reasonable guess for this.

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
will be writing a python function which takes in a "unit" and calls the
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

    2020-04-14 00:01:57 [INFO] idaes.init.fs.H101.control_volume: Initialization Complete
    2020-04-14 00:01:57 [INFO] idaes.init.fs.H101: Initialization Step 1 Complete.
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: ******************************************************************************
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: contain the following acknowledgement:
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: ******************************************************************************
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: Number of nonzeros in equality constraint Jacobian...:      124
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: Number of nonzeros in Lagrangian Hessian.............:      112
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: Total number of variables............................:       41
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: variables with only lower bounds:        0
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: variables with lower and upper bounds:        9
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: variables with only upper bounds:        0
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: Total number of equality constraints.................:       41
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: Total number of inequality constraints...............:        0
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: inequality constraints with only lower bounds:        0
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: inequality constraints with only upper bounds:        0
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 0  0.0000000e+00 1.44e+05 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 1  0.0000000e+00 8.53e+04 1.03e+01  -1.0 3.65e+04    -  1.44e-01 5.96e-01h  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 2  0.0000000e+00 5.59e+04 4.56e+02  -1.0 1.46e+04    -  9.90e-01 3.84e-01h  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 3  0.0000000e+00 5.46e+04 2.28e+04  -1.0 9.01e+03    -  9.64e-01 2.49e-02h  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 4  0.0000000e+00 5.45e+04 8.50e+07  -1.0 8.79e+03    -  9.91e-01 2.77e-04h  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 5r 0.0000000e+00 5.45e+04 1.00e+03   0.7 0.00e+00    -  0.00e+00 3.46e-07R  4
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 6r 0.0000000e+00 4.36e+04 3.24e+03   0.7 2.01e+04    -  9.16e-02 2.64e-03f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 7r 0.0000000e+00 3.79e+04 5.91e+03   0.7 2.19e+04    -  4.65e-02 8.63e-02f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 8r 0.0000000e+00 3.17e+04 5.52e+03   0.7 2.00e+04    -  1.58e-01 1.21e-01f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 9r 0.0000000e+00 2.24e+04 4.07e+03   0.7 1.69e+04    -  2.62e-01 2.16e-01f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 10r 0.0000000e+00 2.81e+03 6.09e+03   0.7 1.32e+04    -  1.00e+00 7.31e-01f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 11r 0.0000000e+00 1.14e+03 5.67e+02   0.7 1.61e+03    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 12r 0.0000000e+00 1.94e+02 4.92e+02  -0.0 7.79e+02    -  1.00e+00 9.01e-01f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 13r 0.0000000e+00 2.70e+01 8.63e+03  -0.0 4.68e+02    -  8.57e-01 3.24e-01f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 14r 0.0000000e+00 2.21e+01 2.68e+04  -0.0 1.19e+03    -  1.00e+00 4.58e-01f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 15r 0.0000000e+00 1.50e+01 2.79e+02  -0.0 2.96e+02    -  1.00e+00 1.00e+00f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 16r 0.0000000e+00 9.88e+00 1.66e+01  -0.0 1.08e+02    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 17r 0.0000000e+00 3.29e+00 1.77e+02  -1.4 5.95e+01    -  7.72e-01 9.68e-01f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 18r 0.0000000e+00 5.49e+02 4.70e+03  -1.4 8.66e+03    -  7.57e-01 3.13e-01f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 19r 0.0000000e+00 3.25e+03 1.47e+04  -1.4 2.57e+03    -  1.00e+00 9.97e-01f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 20r 0.0000000e+00 7.66e+01 9.65e+02  -1.4 1.06e+03    -  4.61e-01 1.00e+00h  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 21r 0.0000000e+00 2.47e+00 2.28e+01  -1.4 5.47e+01    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 22r 0.0000000e+00 8.50e-02 3.54e-02  -1.4 1.24e+00    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 23r 0.0000000e+00 4.43e-01 6.74e+01  -4.7 1.56e+03    -  8.51e-01 8.92e-01f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 24r 0.0000000e+00 1.18e+01 6.96e+03  -4.7 1.12e+02  -4.0 7.45e-01 8.89e-01f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 25r 0.0000000e+00 2.53e+03 6.99e+03  -4.7 9.70e+05    -  1.34e-02 5.52e-03f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 26r 0.0000000e+00 2.53e+03 1.30e+04  -4.7 2.82e+05    -  1.56e-01 1.94e-05f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 27r 0.0000000e+00 2.53e+03 1.41e+04  -4.7 3.91e+05    -  2.35e-01 5.39e-02f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 28r 0.0000000e+00 2.48e+03 9.30e+04  -4.7 3.89e+05    -  9.70e-01 1.83e-02f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 29r 0.0000000e+00 2.48e+03 9.80e+04  -4.7 1.06e+05    -  1.00e+00 3.51e-03f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 30r 0.0000000e+00 6.31e+02 1.61e+05  -4.7 4.03e+03    -  1.00e+00 7.63e-01f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 31r 0.0000000e+00 6.03e+02 8.31e+05  -4.7 9.55e+02    -  1.00e+00 8.77e-02f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 32r 0.0000000e+00 1.66e+02 7.40e+05  -4.7 8.71e+02    -  1.00e+00 7.25e-01f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 33r 0.0000000e+00 1.13e+02 8.82e+05  -4.7 2.39e+02    -  1.00e+00 3.18e-01f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 34r 0.0000000e+00 1.13e+01 7.38e+07  -4.7 1.63e+02    -  1.00e+00 9.70e-01f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 35r 0.0000000e+00 9.91e+00 6.88e+07  -4.7 1.16e+00  -1.8 9.77e-01 6.84e-01h  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 36r 0.0000000e+00 9.90e+00 3.24e+08  -4.7 2.52e+02    -  1.00e+00 5.62e-04h  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 37r 0.0000000e+00 8.52e+00 2.84e+08  -4.7 4.92e+00    -  1.17e-01 1.39e-01f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 38r 0.0000000e+00 1.25e+00 2.52e+08  -4.7 4.24e+00    -  2.12e-01 1.00e+00f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 39r 0.0000000e+00 8.16e-02 2.46e+07  -4.7 7.06e-01    -  9.33e-01 1.00e+00h  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 40r 0.0000000e+00 3.37e-01 2.12e+06  -4.7 8.24e+00    -  9.11e-01 5.50e-01f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 41r 0.0000000e+00 3.16e-01 1.68e+07  -4.7 1.05e+02    -  8.97e-01 1.40e-01f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 42r 0.0000000e+00 1.59e-01 4.30e+06  -4.7 4.26e+01    -  1.00e+00 5.00e-01f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 43r 0.0000000e+00 1.31e-01 1.03e+08  -4.7 3.05e+00    -  1.00e+00 2.23e-01f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 44r 0.0000000e+00 7.23e-02 7.07e+08  -4.7 9.03e-01    -  9.95e-02 4.32e-01f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 45r 0.0000000e+00 2.99e-02 3.02e+08  -4.7 7.49e-02    -  1.00e+00 5.76e-01f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 46r 0.0000000e+00 3.96e-03 7.36e+05  -4.7 4.38e-02    -  1.00e+00 1.00e+00f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 47r 0.0000000e+00 3.96e-03 1.08e+04  -4.7 3.69e-02    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 48r 0.0000000e+00 3.96e-03 1.35e-01  -4.7 1.37e-04    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 49r 0.0000000e+00 3.96e-03 7.39e+05  -7.0 8.44e-01    -  1.00e+00 9.45e-01f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: 50r 0.0000000e+00 3.44e-05 9.17e+05  -7.0 6.03e+04    -  1.00e+00 6.57e-02f  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: Number of Iterations....: 50
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: (scaled)                 (unscaled)
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: Constraint violation....:   2.7355712726718819e-09    3.4409734794849101e-05
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: Overall NLP error.......:   2.7355712726718819e-09    3.4409734794849101e-05
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: Number of objective function evaluations             = 55
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: Number of objective gradient evaluations             = 7
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: Number of equality constraint evaluations            = 55
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: Number of equality constraint Jacobian evaluations   = 52
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: Number of Lagrangian Hessian evaluations             = 50
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.050
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: Total CPU secs in NLP function evaluations           =      0.002
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.H101: EXIT: Optimal Solution Found.
    2020-04-14 00:01:57 [INFO] idaes.init.fs.H101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:01:57 [INFO] idaes.init.fs.H101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:01:57 [INFO] idaes.init.fs.R101.control_volume: Initialization Complete
    2020-04-14 00:01:57 [INFO] idaes.init.fs.R101: Initialization Step 1 Complete.
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: ******************************************************************************
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: contain the following acknowledgement:
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: ******************************************************************************
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: Number of nonzeros in equality constraint Jacobian...:       93
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: Number of nonzeros in Lagrangian Hessian.............:       61
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: Total number of variables............................:       39
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: variables with only lower bounds:        0
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: variables with lower and upper bounds:       10
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: variables with only upper bounds:        0
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: Total number of equality constraints.................:       39
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: Total number of inequality constraints...............:        0
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: inequality constraints with only lower bounds:        0
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: inequality constraints with only upper bounds:        0
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: 0  0.0000000e+00 3.38e+04 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: 1  0.0000000e+00 2.99e+06 2.44e+02  -1.0 1.10e+05    -  7.92e-02 2.48e-01f  3
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: 2  0.0000000e+00 2.00e+06 9.88e+04  -1.0 9.86e+04    -  6.04e-01 9.90e-01H  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: 3  0.0000000e+00 3.15e+04 1.42e+05  -1.0 2.84e+03    -  9.90e-01 9.90e-01h  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: 4  0.0000000e+00 5.77e+04 4.90e+04  -1.0 5.46e+04    -  9.90e-01 1.00e+00h  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: 5  0.0000000e+00 7.12e-04 1.85e-01  -1.0 5.77e+04    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: 6  0.0000000e+00 2.24e-08 2.37e-09  -7.0 5.26e-07    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: Number of Iterations....: 6
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: (scaled)                 (unscaled)
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: Constraint violation....:   3.4711595092880110e-11    2.2351741790771484e-08
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: Overall NLP error.......:   3.4711595092880110e-11    2.2351741790771484e-08
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: Number of objective function evaluations             = 11
    2020-04-14 00:01:57 [DEBUG] idaes.solve.fs.R101: Number of objective gradient evaluations             = 7
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.R101: Number of equality constraint evaluations            = 14
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.R101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.R101: Number of equality constraint Jacobian evaluations   = 7
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.R101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.R101: Number of Lagrangian Hessian evaluations             = 6
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.R101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.002
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.R101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.R101: EXIT: Optimal Solution Found.
    2020-04-14 00:01:58 [INFO] idaes.init.fs.R101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:01:58 [INFO] idaes.init.fs.R101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:01:58 [INFO] idaes.init.fs.F101.control_volume: Initialization Complete
    2020-04-14 00:01:58 [INFO] idaes.init.fs.F101: Initialization Step 1 Complete.
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: ******************************************************************************
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: contain the following acknowledgement:
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: ******************************************************************************
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: Number of nonzeros in equality constraint Jacobian...:      124
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: Number of nonzeros in Lagrangian Hessian.............:      112
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: Total number of variables............................:       41
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: variables with only lower bounds:        0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: variables with lower and upper bounds:        9
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: variables with only upper bounds:        0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: Total number of equality constraints.................:       41
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: Total number of inequality constraints...............:        0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: inequality constraints with only lower bounds:        0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: inequality constraints with only upper bounds:        0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: 0  0.0000000e+00 1.17e+05 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: 1  0.0000000e+00 1.28e+05 8.02e+01  -1.0 8.03e+04    -  1.39e-02 1.19e-01f  1
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: 2  0.0000000e+00 1.16e+05 2.15e+02  -1.0 2.81e+04    -  3.47e-03 1.32e-01f  2
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: 3  0.0000000e+00 1.07e+05 1.99e+02  -1.0 4.03e+04    -  9.89e-01 7.41e-02h  4
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: 4  0.0000000e+00 1.06e+05 1.98e+02  -1.0 3.57e+04    -  5.06e-01 5.95e-03h  8
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: 5  0.0000000e+00 1.06e+05 1.96e+02  -1.0 3.54e+04    -  9.90e-01 6.03e-03h  8
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: 6  0.0000000e+00 1.05e+05 1.94e+02  -1.0 3.51e+04    -  5.89e-01 6.11e-03h  8
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: 7  0.0000000e+00 1.04e+05 1.93e+02  -1.0 3.48e+04    -  1.00e+00 6.18e-03h  8
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: 8  0.0000000e+00 1.04e+05 1.92e+02  -1.7 3.45e+04    -  4.11e-01 6.26e-03h  8
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: 9  0.0000000e+00 1.03e+05 1.92e+02  -1.7 3.42e+04    -  1.00e+00 6.34e-03h  8
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: 10  0.0000000e+00 1.02e+05 1.91e+02  -1.7 3.40e+04    -  6.84e-01 6.41e-03h  8
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: 11  0.0000000e+00 1.02e+05 1.90e+02  -1.7 3.37e+04    -  1.00e+00 6.49e-03h  8
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: 12  0.0000000e+00 1.01e+05 1.89e+02  -1.7 3.34e+04    -  6.89e-01 6.56e-03h  8
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: 13  0.0000000e+00 1.00e+05 1.88e+02  -1.7 3.31e+04    -  1.00e+00 6.63e-03h  8
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: 14  0.0000000e+00 9.96e+04 1.87e+02  -1.7 3.29e+04    -  7.27e-01 6.70e-03h  8
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: 15  0.0000000e+00 9.89e+04 1.86e+02  -1.7 3.26e+04    -  1.00e+00 6.77e-03h  8
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: 16  0.0000000e+00 9.82e+04 1.84e+02  -1.7 3.23e+04    -  7.65e-01 6.84e-03h  8
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: 17  0.0000000e+00 9.75e+04 1.83e+02  -1.7 3.21e+04    -  1.00e+00 6.91e-03h  8
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: 18  0.0000000e+00 1.20e+04 2.74e+02  -1.7 3.18e+04    -  8.05e-01 8.94e-01w  1
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: 19  0.0000000e+00 4.15e+03 1.32e+03  -1.7 2.18e+03    -  6.55e-02 8.35e-01w  1
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: 20  0.0000000e+00 1.84e+02 5.93e+01  -1.7 4.30e+02    -  8.96e-01 9.70e-01h  1
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: 21  0.0000000e+00 1.24e-01 1.12e+01  -1.7 1.39e+01    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: 22  0.0000000e+00 5.63e-08 6.77e-03  -3.8 1.52e-03    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: Number of Iterations....: 22
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: (scaled)                 (unscaled)
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: Constraint violation....:   2.0565647774824511e-11    5.6330463849008083e-08
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: Overall NLP error.......:   2.0565647774824511e-11    5.6330463849008083e-08
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: Number of objective function evaluations             = 157
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: Number of objective gradient evaluations             = 23
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: Number of equality constraint evaluations            = 157
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: Number of equality constraint Jacobian evaluations   = 23
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: Number of Lagrangian Hessian evaluations             = 22
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.007
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: Total CPU secs in NLP function evaluations           =      0.001
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F101: EXIT: Optimal Solution Found.
    2020-04-14 00:01:58 [INFO] idaes.init.fs.F101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:01:58 [INFO] idaes.init.fs.F101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: Ipopt 3.13.2:
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: ******************************************************************************
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: contain the following acknowledgement:
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: ******************************************************************************
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: Number of nonzeros in equality constraint Jacobian...:       29
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: Number of nonzeros in Lagrangian Hessian.............:        0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: Total number of variables............................:       21
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: variables with only lower bounds:        0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: variables with lower and upper bounds:       20
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: variables with only upper bounds:        0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: Total number of equality constraints.................:       21
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: Total number of inequality constraints...............:        0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: inequality constraints with only lower bounds:        0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: inequality constraints with only upper bounds:        0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: 0  0.0000000e+00 3.00e-01 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: 1  0.0000000e+00 3.00e-03 2.20e-06  -1.0 3.00e-01    -  9.90e-01 9.90e-01h  1
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: 2  0.0000000e+00 2.99e-05 2.00e-02  -1.0 3.00e-03    -  9.90e-01 9.90e-01h  1
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: 3  0.0000000e+00 2.40e-07 2.01e+02  -1.0 2.99e-05    -  9.94e-01 9.92e-01h  1
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: 4  0.0000000e+00 4.14e-25 1.05e-10  -1.7 2.40e-07    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: Number of Iterations....: 4
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: (scaled)                 (unscaled)
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: Constraint violation....:   4.1359030627651384e-25    4.1359030627651384e-25
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: Overall NLP error.......:   4.1359030627651384e-25    4.1359030627651384e-25
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: Number of objective function evaluations             = 5
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: Number of objective gradient evaluations             = 5
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: Number of equality constraint evaluations            = 5
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: Number of equality constraint Jacobian evaluations   = 5
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: Number of Lagrangian Hessian evaluations             = 4
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.001
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: EXIT: Optimal Solution Found.
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.S101: 
    2020-04-14 00:01:58 [INFO] idaes.init.fs.S101: Initialization Step 2 Complete: optimal - Optimal Solution Found
    2020-04-14 00:01:58 [INFO] idaes.init.fs.F102.control_volume: Initialization Complete
    2020-04-14 00:01:58 [INFO] idaes.init.fs.F102: Initialization Step 1 Complete.
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: ******************************************************************************
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: This version of Ipopt was compiled from source code available at
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: contain the following acknowledgement:
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: ******************************************************************************
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: Number of nonzeros in equality constraint Jacobian...:      124
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: Number of nonzeros in Lagrangian Hessian.............:      112
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: Total number of variables............................:       41
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: variables with only lower bounds:        0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: variables with lower and upper bounds:        9
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: variables with only upper bounds:        0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: Total number of equality constraints.................:       41
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: Total number of inequality constraints...............:        0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: inequality constraints with only lower bounds:        0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: inequality constraints with only upper bounds:        0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: 0  0.0000000e+00 6.83e+04 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: 1  0.0000000e+00 2.82e+04 5.69e+00  -1.0 2.00e+05    -  1.98e-01 8.56e-01h  1
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: 2  0.0000000e+00 8.44e+03 1.15e+02  -1.0 2.88e+04    -  8.76e-01 9.57e-01h  1
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: 3  0.0000000e+00 1.46e+03 8.84e+04  -1.0 1.98e+03    -  4.14e-01 9.65e-01h  1
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: 4  0.0000000e+00 1.42e+03 1.94e+06  -1.0 9.55e+02    -  9.91e-01 1.00e+00h  1
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: 5  0.0000000e+00 9.65e-01 8.14e+04  -1.0 1.42e+01    -  9.91e-01 1.00e+00h  1
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: 6  0.0000000e+00 2.25e-04 9.29e+01  -1.0 3.62e-01    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: 7  0.0000000e+00 1.46e-11 1.72e-04  -2.5 8.18e-06    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: Number of Iterations....: 7
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: (scaled)                 (unscaled)
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: Constraint violation....:   2.2737367544323206e-13    1.4551915228366852e-11
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: Overall NLP error.......:   2.2737367544323206e-13    1.4551915228366852e-11
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: Number of objective function evaluations             = 8
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: Number of objective gradient evaluations             = 8
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: Number of equality constraint evaluations            = 8
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: Number of inequality constraint evaluations          = 0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: Number of equality constraint Jacobian evaluations   = 8
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: Number of Lagrangian Hessian evaluations             = 7
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: Total CPU secs in IPOPT (w/o function evaluations)   =      0.022
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:01:58 [DEBUG] idaes.solve.fs.F102: EXIT: Optimal Solution Found.
    2020-04-14 00:01:58 [INFO] idaes.init.fs.F102: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:01:58 [INFO] idaes.init.fs.F102: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:01:59 [INFO] idaes.init.fs.C101.control_volume: Initialization Complete
    2020-04-14 00:01:59 [INFO] idaes.init.fs.C101: Initialization Step 1 Complete.
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: ******************************************************************************
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: contain the following acknowledgement:
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: ******************************************************************************
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: Number of nonzeros in equality constraint Jacobian...:       74
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: Number of nonzeros in Lagrangian Hessian.............:       61
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: Total number of variables............................:       32
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: variables with only lower bounds:        0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: variables with lower and upper bounds:        9
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: variables with only upper bounds:        0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: Total number of equality constraints.................:       32
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: Total number of inequality constraints...............:        0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: inequality constraints with only lower bounds:        0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: inequality constraints with only upper bounds:        0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: 0  0.0000000e+00 9.18e+02 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: 1  0.0000000e+00 1.31e+01 1.60e-03  -1.0 2.28e+01    -  9.90e-01 9.90e-01h  1
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: 2  0.0000000e+00 1.32e-01 9.90e+00  -1.0 2.25e+01    -  1.00e+00 9.90e-01h  1
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: 3  0.0000000e+00 2.35e-07 2.83e-11  -1.0 2.26e-01    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: Number of Iterations....: 3
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: (scaled)                 (unscaled)
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: Constraint violation....:   1.7853365888371418e-10    2.3471693566534668e-07
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: Overall NLP error.......:   1.7853365888371418e-10    2.3471693566534668e-07
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: Number of objective function evaluations             = 4
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: Number of objective gradient evaluations             = 4
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: Number of equality constraint evaluations            = 4
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: Number of equality constraint Jacobian evaluations   = 4
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: Number of Lagrangian Hessian evaluations             = 3
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.002
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.C101: EXIT: Optimal Solution Found.
    2020-04-14 00:01:59 [INFO] idaes.init.fs.C101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:01:59 [INFO] idaes.init.fs.C101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: Ipopt 3.13.2:
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: ******************************************************************************
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: contain the following acknowledgement:
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: ******************************************************************************
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: Number of nonzeros in equality constraint Jacobian...:      117
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: Number of nonzeros in Lagrangian Hessian.............:       63
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: Total number of variables............................:       53
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: variables with only lower bounds:        0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: variables with lower and upper bounds:       10
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: variables with only upper bounds:        0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: Total number of equality constraints.................:       53
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: Total number of inequality constraints...............:        0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: inequality constraints with only lower bounds:        0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: inequality constraints with only upper bounds:        0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: 0  0.0000000e+00 3.50e+05 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: 1  0.0000000e+00 2.81e+03 2.79e+00  -1.0 3.50e+05    -  9.85e-01 9.92e-01h  1
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: 2  0.0000000e+00 3.49e+00 1.01e+01  -1.0 2.81e+03    -  9.90e-01 1.00e+00h  1
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: 3  0.0000000e+00 5.90e-05 2.93e+01  -1.0 3.24e-02    -  9.93e-01 1.00e+00h  1
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: 4  0.0000000e+00 7.45e-09 7.21e-13  -1.0 5.90e-05    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: Number of Iterations....: 4
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: (scaled)                 (unscaled)
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: Constraint violation....:   2.9103830456733704e-11    7.4505805969238281e-09
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: Overall NLP error.......:   2.9103830456733704e-11    7.4505805969238281e-09
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: Number of objective function evaluations             = 5
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: Number of objective gradient evaluations             = 5
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: Number of equality constraint evaluations            = 5
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: Number of equality constraint Jacobian evaluations   = 5
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: Number of Lagrangian Hessian evaluations             = 4
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.002
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: EXIT: Optimal Solution Found.
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.M101: 
    2020-04-14 00:01:59 [INFO] idaes.init.fs.M101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:01:59 [INFO] idaes.init.fs.H101.control_volume: Initialization Complete
    2020-04-14 00:01:59 [INFO] idaes.init.fs.H101: Initialization Step 1 Complete.
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: ******************************************************************************
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: contain the following acknowledgement:
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: ******************************************************************************
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: Number of nonzeros in equality constraint Jacobian...:      124
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: Number of nonzeros in Lagrangian Hessian.............:      112
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: Total number of variables............................:       41
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: variables with only lower bounds:        0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: variables with lower and upper bounds:        9
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: variables with only upper bounds:        0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: Total number of equality constraints.................:       41
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: Total number of inequality constraints...............:        0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: inequality constraints with only lower bounds:        0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: inequality constraints with only upper bounds:        0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: 0  0.0000000e+00 1.05e+05 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: 1  0.0000000e+00 5.93e+04 1.13e+01  -1.0 2.36e+04    -  9.10e-02 5.64e-01h  1
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: 2  0.0000000e+00 2.98e+04 4.03e+00  -1.0 1.03e+04    -  9.90e-01 5.47e-01h  1
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: 3  0.0000000e+00 1.41e+04 1.06e+04  -1.0 4.72e+03    -  9.90e-01 5.55e-01h  1
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: 4  0.0000000e+00 1.46e+03 1.62e+08  -1.0 2.11e+03    -  9.94e-01 1.00e+00h  1
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: 5  0.0000000e+00 3.85e+02 7.63e+06  -1.0 2.44e+01    -  1.00e+00 5.12e-01h  1
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: 6  0.0000000e+00 1.80e+02 3.34e+06  -1.0 1.19e+01    -  1.00e+00 4.96e-01h  2
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: 7  0.0000000e+00 1.30e+02 1.36e+10  -1.0 6.01e+00    -  1.00e+00 9.92e-01h  1
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: 8  0.0000000e+00 1.01e+02 5.67e+09  -1.0 4.72e-02    -  1.00e+00 6.99e-01h  1
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: 9  0.0000000e+00 2.57e+01 1.29e+11  -1.0 1.42e-02    -  1.00e+00 5.00e-01f  2
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: 10  0.0000000e+00 2.48e+01 2.15e+11  -1.0 7.10e-03    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: 11  0.0000000e+00 2.00e-03 1.38e+07  -1.0 2.77e-06    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: 12  0.0000000e+00 1.46e-11 4.94e+02  -1.7 1.35e-07    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: Number of Iterations....: 12
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: (scaled)                 (unscaled)
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: Constraint violation....:   7.1054273576010019e-14    1.4551915228366852e-11
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: Overall NLP error.......:   7.1054273576010019e-14    1.4551915228366852e-11
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: Number of objective function evaluations             = 19
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: Number of objective gradient evaluations             = 13
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: Number of equality constraint evaluations            = 19
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: Number of equality constraint Jacobian evaluations   = 13
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: Number of Lagrangian Hessian evaluations             = 12
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.004
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:01:59 [DEBUG] idaes.solve.fs.H101: EXIT: Optimal Solution Found.
    2020-04-14 00:01:59 [INFO] idaes.init.fs.H101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:01:59 [INFO] idaes.init.fs.H101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:00 [INFO] idaes.init.fs.R101.control_volume: Initialization Complete
    2020-04-14 00:02:00 [INFO] idaes.init.fs.R101: Initialization Step 1 Complete.
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: ******************************************************************************
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: contain the following acknowledgement:
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: ******************************************************************************
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: Number of nonzeros in equality constraint Jacobian...:       93
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: Number of nonzeros in Lagrangian Hessian.............:       61
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: Total number of variables............................:       39
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: variables with only lower bounds:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: variables with lower and upper bounds:       10
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: variables with only upper bounds:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: Total number of equality constraints.................:       39
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: 0  0.0000000e+00 3.34e+04 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: 1  0.0000000e+00 5.16e+05 3.30e+01  -1.0 9.28e+04    -  1.51e-01 1.24e-01h  4
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: 2  0.0000000e+00 5.14e+05 3.93e+01  -1.0 8.65e+04    -  5.58e-01 1.55e-02h  7
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: 3  0.0000000e+00 5.12e+05 4.27e+01  -1.0 8.58e+04    -  3.40e-01 1.55e-02h  7
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: 4  0.0000000e+00 5.10e+05 4.92e+01  -1.0 8.50e+04    -  6.07e-01 1.55e-02h  7
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: 5  0.0000000e+00 5.08e+05 5.29e+01  -1.0 8.42e+04    -  3.85e-01 1.55e-02h  7
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: 6  0.0000000e+00 5.05e+05 6.40e+01  -1.0 8.34e+04    -  9.90e-01 1.55e-02h  7
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: 7  0.0000000e+00 5.03e+05 6.72e+01  -1.0 8.26e+04    -  3.52e-01 1.55e-02h  7
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: 8  0.0000000e+00 5.00e+05 7.84e+01  -1.0 8.18e+04    -  9.91e-01 1.55e-02h  7
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: 9  0.0000000e+00 4.98e+05 8.20e+01  -1.0 8.10e+04    -  3.96e-01 1.55e-02h  7
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: 10  0.0000000e+00 4.95e+05 9.34e+01  -1.0 8.03e+04    -  1.00e+00 1.55e-02h  7
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: 11  0.0000000e+00 2.18e+07 1.75e+04  -1.0 7.95e+04    -  4.38e-01 9.90e-01w  1
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: 12  0.0000000e+00 2.41e+05 5.85e+03  -1.0 1.79e+04    -  9.35e-01 9.90e-01w  1
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: 13  0.0000000e+00 5.29e+04 2.08e+04  -1.0 4.07e+04    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: 14  0.0000000e+00 2.48e-04 2.22e-02  -1.0 5.29e+04    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: 15  0.0000000e+00 5.22e-08 9.54e-10  -7.0 1.76e-07    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: Number of Iterations....: 15
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: (scaled)                 (unscaled)
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: Constraint violation....:   2.0826957055728064e-11    5.2154064178466797e-08
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: Overall NLP error.......:   2.0826957055728064e-11    5.2154064178466797e-08
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: Number of objective function evaluations             = 103
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: Number of objective gradient evaluations             = 16
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: Number of equality constraint evaluations            = 113
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: Number of equality constraint Jacobian evaluations   = 16
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: Number of Lagrangian Hessian evaluations             = 15
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.008
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: Total CPU secs in NLP function evaluations           =      0.001
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.R101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:00 [INFO] idaes.init.fs.R101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:02:00 [INFO] idaes.init.fs.R101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:00 [INFO] idaes.init.fs.F101.control_volume: Initialization Complete
    2020-04-14 00:02:00 [INFO] idaes.init.fs.F101: Initialization Step 1 Complete.
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: ******************************************************************************
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: contain the following acknowledgement:
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: ******************************************************************************
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: Number of nonzeros in equality constraint Jacobian...:      124
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: Number of nonzeros in Lagrangian Hessian.............:      112
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: Total number of variables............................:       41
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: variables with only lower bounds:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: variables with lower and upper bounds:        9
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: variables with only upper bounds:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: Total number of equality constraints.................:       41
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: 0  0.0000000e+00 8.74e+04 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: 1  0.0000000e+00 7.30e+04 3.10e+01  -1.0 1.09e+04    -  4.52e-02 3.85e-01f  1
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: 2  0.0000000e+00 3.99e+04 3.26e+01  -1.0 7.60e+03    -  2.91e-01 4.83e-01h  1
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: 3  0.0000000e+00 6.21e+02 1.64e+01  -1.0 3.50e+03    -  9.90e-01 1.00e+00H  1
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: 4  0.0000000e+00 9.72e-01 8.42e+00  -1.0 1.57e+01    -  9.90e-01 1.00e+00h  1
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: 5  0.0000000e+00 6.12e-06 6.99e+02  -1.0 3.15e-02    -  9.96e-01 1.00e+00h  1
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: Number of Iterations....: 5
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: (scaled)                 (unscaled)
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: Constraint violation....:   2.4734746446097313e-09    6.1192513385321945e-06
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: Overall NLP error.......:   2.4734746446097313e-09    6.1192513385321945e-06
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: Number of objective function evaluations             = 7
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: Number of objective gradient evaluations             = 6
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: Number of equality constraint evaluations            = 7
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: Number of equality constraint Jacobian evaluations   = 6
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: Number of Lagrangian Hessian evaluations             = 5
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.002
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.F101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:00 [INFO] idaes.init.fs.F101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:02:00 [INFO] idaes.init.fs.F101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: Ipopt 3.13.2:
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: ******************************************************************************
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: contain the following acknowledgement:
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: ******************************************************************************
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: Number of nonzeros in equality constraint Jacobian...:       29
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: Number of nonzeros in Lagrangian Hessian.............:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: Total number of variables............................:       21
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: variables with only lower bounds:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: variables with lower and upper bounds:       20
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: variables with only upper bounds:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: Total number of equality constraints.................:       21
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: 0  0.0000000e+00 1.00e-02 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: 1  0.0000000e+00 1.00e-04 2.20e-06  -1.0 1.00e-02    -  9.90e-01 9.90e-01h  1
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: 2  0.0000000e+00 9.98e-07 2.00e-02  -1.0 1.00e-04    -  9.90e-01 9.90e-01h  1
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: 3  0.0000000e+00 8.00e-09 2.01e+02  -1.0 9.98e-07    -  9.94e-01 9.92e-01h  1
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: Number of Iterations....: 3
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: (scaled)                 (unscaled)
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: Constraint violation....:   7.9999999999999045e-09    7.9999999999999045e-09
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: Overall NLP error.......:   7.9999999999999045e-09    7.9999999999999045e-09
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: Number of objective function evaluations             = 4
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: Number of objective gradient evaluations             = 4
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: Number of equality constraint evaluations            = 4
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: Number of equality constraint Jacobian evaluations   = 4
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: Number of Lagrangian Hessian evaluations             = 3
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.001
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.S101: 
    2020-04-14 00:02:00 [INFO] idaes.init.fs.S101: Initialization Step 2 Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:00 [INFO] idaes.init.fs.C101.control_volume: Initialization Complete
    2020-04-14 00:02:00 [INFO] idaes.init.fs.C101: Initialization Step 1 Complete.
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: ******************************************************************************
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: contain the following acknowledgement:
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: ******************************************************************************
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: Number of nonzeros in equality constraint Jacobian...:       74
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: Number of nonzeros in Lagrangian Hessian.............:       61
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: Total number of variables............................:       32
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: variables with only lower bounds:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: variables with lower and upper bounds:        9
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: variables with only upper bounds:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: Total number of equality constraints.................:       32
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: 0  0.0000000e+00 3.63e+02 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: 1  0.0000000e+00 1.41e-01 1.60e-03  -1.0 3.56e+00    -  9.90e-01 9.90e-01h  1
    2020-04-14 00:02:00 [DEBUG] idaes.solve.fs.C101: 2  0.0000000e+00 1.03e-03 9.84e+00  -1.0 3.52e+00    -  1.00e+00 9.90e-01h  1
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.C101: 3  0.0000000e+00 3.51e-08 8.70e-13  -1.0 3.50e-02    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.C101: Number of Iterations....: 3
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.C101: (scaled)                 (unscaled)
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.C101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.C101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.C101: Constraint violation....:   4.8173230335450958e-11    3.5146513255313039e-08
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.C101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.C101: Overall NLP error.......:   4.8173230335450958e-11    3.5146513255313039e-08
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.C101: Number of objective function evaluations             = 4
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.C101: Number of objective gradient evaluations             = 4
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.C101: Number of equality constraint evaluations            = 4
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.C101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.C101: Number of equality constraint Jacobian evaluations   = 4
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.C101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.C101: Number of Lagrangian Hessian evaluations             = 3
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.C101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.002
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.C101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.C101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:01 [INFO] idaes.init.fs.C101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:02:01 [INFO] idaes.init.fs.C101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: Ipopt 3.13.2:
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: ******************************************************************************
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: contain the following acknowledgement:
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: ******************************************************************************
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: Number of nonzeros in equality constraint Jacobian...:      117
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: Number of nonzeros in Lagrangian Hessian.............:       63
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: Total number of variables............................:       53
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: variables with only lower bounds:        0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: variables with lower and upper bounds:       10
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: variables with only upper bounds:        0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: Total number of equality constraints.................:       53
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: 0  0.0000000e+00 2.24e+02 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: 1  0.0000000e+00 2.98e+02 5.73e-01  -1.0 2.35e+02    -  9.88e-01 9.92e-01h  1
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: 2  0.0000000e+00 9.91e-01 9.97e+00  -1.0 1.69e+01    -  9.90e-01 1.00e+00h  1
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: 3  0.0000000e+00 1.79e-07 2.66e+01  -1.0 1.15e-02    -  9.93e-01 1.00e+00h  1
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: Number of Iterations....: 3
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: (scaled)                 (unscaled)
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: Constraint violation....:   2.0915153662155086e-10    1.7881393432617188e-07
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: Overall NLP error.......:   2.0915153662155086e-10    1.7881393432617188e-07
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: Number of objective function evaluations             = 4
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: Number of objective gradient evaluations             = 4
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: Number of equality constraint evaluations            = 4
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: Number of equality constraint Jacobian evaluations   = 4
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: Number of Lagrangian Hessian evaluations             = 3
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.001
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.M101: 
    2020-04-14 00:02:01 [INFO] idaes.init.fs.M101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:01 [INFO] idaes.init.fs.H101.control_volume: Initialization Complete
    2020-04-14 00:02:01 [INFO] idaes.init.fs.H101: Initialization Step 1 Complete.
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: ******************************************************************************
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: contain the following acknowledgement:
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: ******************************************************************************
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: Number of nonzeros in equality constraint Jacobian...:      124
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: Number of nonzeros in Lagrangian Hessian.............:      112
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: Total number of variables............................:       41
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: variables with only lower bounds:        0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: variables with lower and upper bounds:        9
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: variables with only upper bounds:        0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: Total number of equality constraints.................:       41
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: 0  0.0000000e+00 5.41e+04 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: 1  0.0000000e+00 2.71e+04 1.18e+01  -1.0 1.40e+04    -  6.37e-02 5.29e-01h  1
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: 2  0.0000000e+00 1.27e+04 5.28e+00  -1.0 8.31e+03    -  9.90e-01 5.21e-01h  1
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: 3  0.0000000e+00 5.86e+03 9.25e+03  -1.0 4.22e+03    -  9.90e-01 5.33e-01h  1
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: 4  0.0000000e+00 1.95e+03 4.39e+08  -1.0 2.04e+03    -  9.94e-01 1.00e+00h  1
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: 5  0.0000000e+00 6.46e+01 3.12e+06  -1.0 1.93e+01    -  1.00e+00 5.00e-01h  1
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: 6  0.0000000e+00 3.02e+01 4.75e+06  -1.0 9.65e+00    -  1.00e+00 4.96e-01h  2
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: 7  0.0000000e+00 2.11e+01 3.98e+10  -1.0 4.86e+00    -  1.00e+00 9.94e-01h  1
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: 8  0.0000000e+00 5.93e+00 3.45e+09  -1.0 2.80e-02    -  1.00e+00 5.00e-01h  2
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: 9  0.0000000e+00 4.59e+00 2.70e+10  -1.0 1.40e-02    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: 10  0.0000000e+00 1.11e-04 3.02e+05  -1.0 4.56e-07    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: 11  0.0000000e+00 1.09e-11 2.18e+01  -3.8 3.46e-08    -  1.00e+00 1.00e+00   0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: Number of Iterations....: 11
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: (scaled)                 (unscaled)
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: Constraint violation....:   2.2737367544323206e-13    1.0913936421275139e-11
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: Overall NLP error.......:   2.2737367544323206e-13    1.0913936421275139e-11
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: Number of objective function evaluations             = 17
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: Number of objective gradient evaluations             = 12
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: Number of equality constraint evaluations            = 17
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: Number of equality constraint Jacobian evaluations   = 12
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: Number of Lagrangian Hessian evaluations             = 11
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.011
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.H101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:01 [INFO] idaes.init.fs.H101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:02:01 [INFO] idaes.init.fs.H101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:01 [INFO] idaes.init.fs.R101.control_volume: Initialization Complete
    2020-04-14 00:02:01 [INFO] idaes.init.fs.R101: Initialization Step 1 Complete.
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: ******************************************************************************
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: contain the following acknowledgement:
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: ******************************************************************************
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: Number of nonzeros in equality constraint Jacobian...:       93
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: Number of nonzeros in Lagrangian Hessian.............:       61
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: Total number of variables............................:       39
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: variables with only lower bounds:        0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: variables with lower and upper bounds:       10
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: variables with only upper bounds:        0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: Total number of equality constraints.................:       39
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: 0  0.0000000e+00 3.25e+04 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: 1  0.0000000e+00 3.24e+04 5.27e+00  -1.0 5.45e+04    -  4.79e-01 6.04e-05h 15
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: 2  0.0000000e+00 3.24e+04 1.18e+01  -1.0 5.45e+04    -  6.18e-01 6.04e-05h 15
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: 3  0.0000000e+00 3.24e+04 1.87e+01  -1.0 5.45e+04    -  6.78e-01 6.04e-05h 15
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: 4  0.0000000e+00 3.24e+04 2.86e+01  -1.0 5.45e+04    -  9.90e-01 6.04e-05h 15
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: 5  0.0000000e+00 3.24e+04 3.48e+01  -1.0 5.45e+04    -  6.15e-01 6.04e-05h 15
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: 6  0.0000000e+00 3.24e+04 4.47e+01  -1.0 5.45e+04    -  9.92e-01 6.04e-05h 15
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: 7  0.0000000e+00 3.24e+04 5.01e+01  -1.0 5.45e+04    -  5.38e-01 6.04e-05h 15
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: 8  0.0000000e+00 3.24e+04 6.01e+01  -1.0 5.45e+04    -  1.00e+00 6.04e-05h 15
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: 9  0.0000000e+00 3.24e+04 6.55e+01  -1.0 5.45e+04    -  5.37e-01 6.04e-05h 15
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: 10  0.0000000e+00 3.24e+04 7.55e+01  -1.0 5.45e+04    -  1.00e+00 6.04e-05h 15
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: 11  0.0000000e+00 1.24e+07 8.23e+01  -1.0 5.45e+04    -  5.37e-01 9.90e-01w  1
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: 12  0.0000000e+00 1.71e+05 9.11e+01  -1.0 7.89e+03    -  1.00e+00 9.90e-01w  1
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: 13  0.0000000e+00 4.02e+04 1.64e+00  -1.0 2.26e+04    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:01 [DEBUG] idaes.solve.fs.R101: 14  0.0000000e+00 3.56e-05 5.20e-02  -3.8 4.02e+04    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.R101: Number of Iterations....: 14
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.R101: (scaled)                 (unscaled)
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.R101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.R101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.R101: Constraint violation....:   1.2334621097335475e-08    3.5621225833892822e-05
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.R101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.R101: Overall NLP error.......:   1.2334621097335475e-08    3.5621225833892822e-05
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.R101: Number of objective function evaluations             = 185
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.R101: Number of objective gradient evaluations             = 15
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.R101: Number of equality constraint evaluations            = 195
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.R101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.R101: Number of equality constraint Jacobian evaluations   = 15
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.R101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.R101: Number of Lagrangian Hessian evaluations             = 14
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.R101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.005
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.R101: Total CPU secs in NLP function evaluations           =      0.001
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.R101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:02 [INFO] idaes.init.fs.R101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:02:02 [INFO] idaes.init.fs.R101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:02 [INFO] idaes.init.fs.F101.control_volume: Initialization Complete
    2020-04-14 00:02:02 [INFO] idaes.init.fs.F101: Initialization Step 1 Complete.
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: ******************************************************************************
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: contain the following acknowledgement:
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: ******************************************************************************
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: Number of nonzeros in equality constraint Jacobian...:      124
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: Number of nonzeros in Lagrangian Hessian.............:      112
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: Total number of variables............................:       41
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: variables with only lower bounds:        0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: variables with lower and upper bounds:        9
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: variables with only upper bounds:        0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: Total number of equality constraints.................:       41
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: 0  0.0000000e+00 5.18e+04 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: 1  0.0000000e+00 3.81e+04 6.10e+00  -1.0 1.67e+04    -  1.31e-01 6.60e-01f  1
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: 2  0.0000000e+00 9.85e+03 2.82e+01  -1.0 5.27e+03    -  8.83e-01 6.24e-01h  1
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: 3  0.0000000e+00 1.72e+02 1.97e+01  -1.0 2.29e+03    -  8.30e-01 1.00e+00h  1
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: 4  0.0000000e+00 1.41e+00 2.45e+00  -1.0 1.54e+01    -  9.90e-01 9.92e-01h  1
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: 5  0.0000000e+00 2.66e-06 6.13e+02  -1.0 1.28e-01    -  9.92e-01 1.00e+00h  1
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: Number of Iterations....: 5
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: (scaled)                 (unscaled)
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: Constraint violation....:   3.3424899880925952e-11    2.6555062504485250e-06
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: Overall NLP error.......:   3.3424899880925952e-11    2.6555062504485250e-06
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: Number of objective function evaluations             = 6
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: Number of objective gradient evaluations             = 6
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: Number of equality constraint evaluations            = 6
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: Number of equality constraint Jacobian evaluations   = 6
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: Number of Lagrangian Hessian evaluations             = 5
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.001
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.F101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:02 [INFO] idaes.init.fs.F101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:02:02 [INFO] idaes.init.fs.F101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: Ipopt 3.13.2:
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: ******************************************************************************
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: contain the following acknowledgement:
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: ******************************************************************************
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: Number of nonzeros in equality constraint Jacobian...:       29
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: Number of nonzeros in Lagrangian Hessian.............:        0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: Total number of variables............................:       21
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: variables with only lower bounds:        0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: variables with lower and upper bounds:       20
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: variables with only upper bounds:        0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: Total number of equality constraints.................:       21
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: 0  0.0000000e+00 1.00e-02 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: 1  0.0000000e+00 1.00e-04 2.20e-06  -1.0 1.00e-02    -  9.90e-01 9.90e-01h  1
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: 2  0.0000000e+00 9.98e-07 2.00e-02  -1.0 1.00e-04    -  9.90e-01 9.90e-01h  1
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: 3  0.0000000e+00 8.00e-09 2.01e+02  -1.0 9.98e-07    -  9.94e-01 9.92e-01h  1
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: Number of Iterations....: 3
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: (scaled)                 (unscaled)
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: Constraint violation....:   7.9999999999999045e-09    7.9999999999999045e-09
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: Overall NLP error.......:   7.9999999999999045e-09    7.9999999999999045e-09
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: Number of objective function evaluations             = 4
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: Number of objective gradient evaluations             = 4
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: Number of equality constraint evaluations            = 4
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: Number of equality constraint Jacobian evaluations   = 4
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: Number of Lagrangian Hessian evaluations             = 3
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.001
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.S101: 
    2020-04-14 00:02:02 [INFO] idaes.init.fs.S101: Initialization Step 2 Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:02 [INFO] idaes.init.fs.C101.control_volume: Initialization Complete
    2020-04-14 00:02:02 [INFO] idaes.init.fs.C101: Initialization Step 1 Complete.
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: ******************************************************************************
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: contain the following acknowledgement:
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: ******************************************************************************
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: Number of nonzeros in equality constraint Jacobian...:       74
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: Number of nonzeros in Lagrangian Hessian.............:       61
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: Total number of variables............................:       32
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: variables with only lower bounds:        0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: variables with lower and upper bounds:        9
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: variables with only upper bounds:        0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: Total number of equality constraints.................:       32
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: 0  0.0000000e+00 1.00e-02 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: 1  0.0000000e+00 1.00e-04 1.60e-03  -1.0 1.00e-02    -  9.90e-01 9.90e-01h  1
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: 2  0.0000000e+00 9.84e-07 9.84e+00  -1.0 1.00e-04    -  1.00e+00 9.90e-01h  1
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: Number of Iterations....: 2
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: (scaled)                 (unscaled)
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: Constraint violation....:   9.8400000480001186e-07    9.8400000480001186e-07
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: Overall NLP error.......:   9.8400000480001186e-07    9.8400000480001186e-07
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: Number of objective function evaluations             = 3
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: Number of objective gradient evaluations             = 3
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: Number of equality constraint evaluations            = 3
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: Number of equality constraint Jacobian evaluations   = 3
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: Number of Lagrangian Hessian evaluations             = 2
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.001
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:02 [DEBUG] idaes.solve.fs.C101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:02 [INFO] idaes.init.fs.C101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:02:02 [INFO] idaes.init.fs.C101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: Ipopt 3.13.2:
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: ******************************************************************************
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: contain the following acknowledgement:
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: ******************************************************************************
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: Number of nonzeros in equality constraint Jacobian...:      117
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: Number of nonzeros in Lagrangian Hessian.............:       63
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: Total number of variables............................:       53
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: variables with only lower bounds:        0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: variables with lower and upper bounds:       10
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: variables with only upper bounds:        0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: Total number of equality constraints.................:       53
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: 0  0.0000000e+00 1.26e+02 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: 1  0.0000000e+00 2.35e+00 3.39e+00  -1.0 5.83e+02    -  7.69e-01 9.92e-01H  1
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: 2  0.0000000e+00 2.43e-01 9.99e+00  -1.0 7.44e+00    -  9.90e-01 1.00e+00h  1
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: 3  0.0000000e+00 7.45e-09 3.70e+01  -1.0 5.30e-03    -  9.90e-01 1.00e+00h  1
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: Number of Iterations....: 3
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: (scaled)                 (unscaled)
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: Constraint violation....:   2.9103830456733704e-11    7.4505805969238281e-09
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: Overall NLP error.......:   2.9103830456733704e-11    7.4505805969238281e-09
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: Number of objective function evaluations             = 5
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: Number of objective gradient evaluations             = 4
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: Number of equality constraint evaluations            = 5
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: Number of equality constraint Jacobian evaluations   = 4
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: Number of Lagrangian Hessian evaluations             = 3
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.002
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.M101: 
    2020-04-14 00:02:03 [INFO] idaes.init.fs.M101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:03 [INFO] idaes.init.fs.H101.control_volume: Initialization Complete
    2020-04-14 00:02:03 [INFO] idaes.init.fs.H101: Initialization Step 1 Complete.
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: ******************************************************************************
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: contain the following acknowledgement:
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: ******************************************************************************
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: Number of nonzeros in equality constraint Jacobian...:      124
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: Number of nonzeros in Lagrangian Hessian.............:      112
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: Total number of variables............................:       41
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: variables with only lower bounds:        0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: variables with lower and upper bounds:        9
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: variables with only upper bounds:        0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: Total number of equality constraints.................:       41
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: 0  0.0000000e+00 5.31e+04 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: 1  0.0000000e+00 2.73e+04 1.13e+01  -1.0 9.78e+03    -  9.31e-02 5.26e-01h  1
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: 2  0.0000000e+00 1.30e+04 5.37e+00  -1.0 4.93e+03    -  9.90e-01 5.20e-01h  1
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: 3  0.0000000e+00 6.03e+03 9.17e+03  -1.0 2.43e+03    -  9.90e-01 5.32e-01h  1
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: 4  0.0000000e+00 1.71e+03 4.38e+08  -1.0 1.15e+03    -  9.94e-01 1.00e+00h  1
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: 5  0.0000000e+00 6.61e+01 3.17e+06  -1.0 1.78e+01    -  1.00e+00 5.00e-01h  1
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: 6  0.0000000e+00 3.12e+01 4.80e+06  -1.0 8.89e+00    -  1.00e+00 4.96e-01h  2
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: 7  0.0000000e+00 1.89e+01 3.99e+10  -1.0 4.48e+00    -  1.00e+00 9.94e-01h  1
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: 8  0.0000000e+00 5.23e+00 3.45e+09  -1.0 2.60e-02    -  1.00e+00 5.00e-01h  2
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: 9  0.0000000e+00 4.19e+00 2.86e+10  -1.0 1.30e-02    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: 10  0.0000000e+00 9.04e-05 2.85e+05  -1.0 4.05e-07    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: Number of Iterations....: 10
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: (scaled)                 (unscaled)
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: Constraint violation....:   1.7768289534621779e-08    9.0412653662497178e-05
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: Overall NLP error.......:   1.7768289534621779e-08    9.0412653662497178e-05
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: Number of objective function evaluations             = 16
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: Number of objective gradient evaluations             = 11
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: Number of equality constraint evaluations            = 16
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: Number of equality constraint Jacobian evaluations   = 11
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: Number of Lagrangian Hessian evaluations             = 10
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.003
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.H101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:03 [INFO] idaes.init.fs.H101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:02:03 [INFO] idaes.init.fs.H101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:03 [INFO] idaes.init.fs.R101.control_volume: Initialization Complete
    2020-04-14 00:02:03 [INFO] idaes.init.fs.R101: Initialization Step 1 Complete.
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: ******************************************************************************
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: contain the following acknowledgement:
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: ******************************************************************************
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: Number of nonzeros in equality constraint Jacobian...:       93
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: Number of nonzeros in Lagrangian Hessian.............:       61
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: Total number of variables............................:       39
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: variables with only lower bounds:        0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: variables with lower and upper bounds:       10
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: variables with only upper bounds:        0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: Total number of equality constraints.................:       39
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: 0  0.0000000e+00 3.26e+04 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: 1  0.0000000e+00 3.26e+04 4.69e+00  -1.0 5.64e+04    -  4.27e-01 6.04e-05h 15
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: 2  0.0000000e+00 3.26e+04 1.11e+01  -1.0 5.64e+04    -  6.10e-01 6.04e-05h 15
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: 3  0.0000000e+00 3.26e+04 1.74e+01  -1.0 5.64e+04    -  6.16e-01 6.04e-05h 15
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: 4  0.0000000e+00 3.26e+04 2.74e+01  -1.0 5.64e+04    -  9.90e-01 6.04e-05h 15
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: 5  0.0000000e+00 3.26e+04 3.18e+01  -1.0 5.64e+04    -  4.40e-01 6.04e-05h 15
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: 6  0.0000000e+00 3.26e+04 4.17e+01  -1.0 5.64e+04    -  9.91e-01 6.04e-05h 15
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: 7  0.0000000e+00 3.26e+04 4.61e+01  -1.0 5.64e+04    -  4.42e-01 6.04e-05h 15
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: 8  0.0000000e+00 3.26e+04 5.61e+01  -1.0 5.64e+04    -  1.00e+00 6.04e-05h 15
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: 9  0.0000000e+00 3.26e+04 6.05e+01  -1.0 5.64e+04    -  4.40e-01 6.04e-05h 15
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: 10  0.0000000e+00 3.26e+04 7.05e+01  -1.0 5.64e+04    -  1.00e+00 6.04e-05h 15
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: 11  0.0000000e+00 1.34e+07 1.55e+02  -1.0 5.64e+04    -  4.40e-01 9.90e-01w  1
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: 12  0.0000000e+00 1.85e+05 1.65e+02  -1.0 8.57e+03    -  1.00e+00 9.90e-01w  1
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: 13  0.0000000e+00 4.16e+04 6.67e+00  -1.0 2.28e+04    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: 14  0.0000000e+00 3.73e-05 5.65e-02  -2.5 4.16e+04    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: Number of Iterations....: 14
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: (scaled)                 (unscaled)
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: Constraint violation....:   1.2935751389848925e-08    3.7349760532379150e-05
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: Overall NLP error.......:   1.2935751389848925e-08    3.7349760532379150e-05
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: Number of objective function evaluations             = 185
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: Number of objective gradient evaluations             = 15
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: Number of equality constraint evaluations            = 195
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: Number of equality constraint Jacobian evaluations   = 15
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: Number of Lagrangian Hessian evaluations             = 14
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.008
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: Total CPU secs in NLP function evaluations           =      0.001
    2020-04-14 00:02:03 [DEBUG] idaes.solve.fs.R101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:03 [INFO] idaes.init.fs.R101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:02:03 [INFO] idaes.init.fs.R101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:03 [INFO] idaes.init.fs.F101.control_volume: Initialization Complete
    2020-04-14 00:02:03 [INFO] idaes.init.fs.F101: Initialization Step 1 Complete.
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: ******************************************************************************
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: contain the following acknowledgement:
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: ******************************************************************************
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: Number of nonzeros in equality constraint Jacobian...:      124
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: Number of nonzeros in Lagrangian Hessian.............:      112
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: Total number of variables............................:       41
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: variables with only lower bounds:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: variables with lower and upper bounds:        9
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: variables with only upper bounds:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: Total number of equality constraints.................:       41
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: 0  0.0000000e+00 4.96e+04 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: 1  0.0000000e+00 3.54e+04 5.11e+00  -1.0 3.23e+03    -  1.46e-01 6.92e-01f  1
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: 2  0.0000000e+00 9.50e+03 2.82e+01  -1.0 2.15e+03    -  8.67e-01 6.00e-01h  1
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: 3  0.0000000e+00 1.87e+02 2.35e+01  -1.0 7.00e+02    -  7.91e-01 1.00e+00h  1
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: 4  0.0000000e+00 1.47e+00 2.76e+00  -1.0 1.47e+01    -  9.90e-01 9.92e-01h  1
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: 5  0.0000000e+00 4.89e-07 7.23e+02  -1.0 1.16e-01    -  9.92e-01 1.00e+00h  1
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: Number of Iterations....: 5
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: (scaled)                 (unscaled)
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: Constraint violation....:   2.0246222515506914e-11    4.8884612624533474e-07
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: Overall NLP error.......:   2.0246222515506914e-11    4.8884612624533474e-07
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: Number of objective function evaluations             = 6
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: Number of objective gradient evaluations             = 6
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: Number of equality constraint evaluations            = 6
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: Number of equality constraint Jacobian evaluations   = 6
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: Number of Lagrangian Hessian evaluations             = 5
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.001
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.F101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:04 [INFO] idaes.init.fs.F101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:02:04 [INFO] idaes.init.fs.F101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: Ipopt 3.13.2:
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: ******************************************************************************
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: contain the following acknowledgement:
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: ******************************************************************************
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: Number of nonzeros in equality constraint Jacobian...:       29
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: Number of nonzeros in Lagrangian Hessian.............:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: Total number of variables............................:       21
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: variables with only lower bounds:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: variables with lower and upper bounds:       20
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: variables with only upper bounds:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: Total number of equality constraints.................:       21
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: 0  0.0000000e+00 1.00e-02 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: 1  0.0000000e+00 1.00e-04 2.20e-06  -1.0 1.00e-02    -  9.90e-01 9.90e-01h  1
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: 2  0.0000000e+00 9.98e-07 2.00e-02  -1.0 1.00e-04    -  9.90e-01 9.90e-01h  1
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: 3  0.0000000e+00 8.00e-09 2.01e+02  -1.0 9.98e-07    -  9.94e-01 9.92e-01h  1
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: Number of Iterations....: 3
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: (scaled)                 (unscaled)
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: Constraint violation....:   7.9999999999999045e-09    7.9999999999999045e-09
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: Overall NLP error.......:   7.9999999999999045e-09    7.9999999999999045e-09
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: Number of objective function evaluations             = 4
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: Number of objective gradient evaluations             = 4
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: Number of equality constraint evaluations            = 4
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: Number of equality constraint Jacobian evaluations   = 4
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: Number of Lagrangian Hessian evaluations             = 3
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.001
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.S101: 
    2020-04-14 00:02:04 [INFO] idaes.init.fs.S101: Initialization Step 2 Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:04 [INFO] idaes.init.fs.C101.control_volume: Initialization Complete
    2020-04-14 00:02:04 [INFO] idaes.init.fs.C101: Initialization Step 1 Complete.
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: ******************************************************************************
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: contain the following acknowledgement:
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: ******************************************************************************
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: Number of nonzeros in equality constraint Jacobian...:       74
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: Number of nonzeros in Lagrangian Hessian.............:       61
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: Total number of variables............................:       32
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: variables with only lower bounds:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: variables with lower and upper bounds:        9
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: variables with only upper bounds:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: Total number of equality constraints.................:       32
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: 0  0.0000000e+00 1.00e-02 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: 1  0.0000000e+00 1.00e-04 1.60e-03  -1.0 1.00e-02    -  9.90e-01 9.90e-01h  1
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: 2  0.0000000e+00 9.84e-07 9.84e+00  -1.0 1.00e-04    -  1.00e+00 9.90e-01h  1
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: Number of Iterations....: 2
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: (scaled)                 (unscaled)
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: Constraint violation....:   9.8400000480001186e-07    9.8400000480001186e-07
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: Overall NLP error.......:   9.8400000480001186e-07    9.8400000480001186e-07
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: Number of objective function evaluations             = 3
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: Number of objective gradient evaluations             = 3
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: Number of equality constraint evaluations            = 3
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: Number of equality constraint Jacobian evaluations   = 3
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: Number of Lagrangian Hessian evaluations             = 2
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.001
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.C101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:04 [INFO] idaes.init.fs.C101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:02:04 [INFO] idaes.init.fs.C101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: Ipopt 3.13.2:
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: ******************************************************************************
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: contain the following acknowledgement:
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: ******************************************************************************
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: Number of nonzeros in equality constraint Jacobian...:      117
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: Number of nonzeros in Lagrangian Hessian.............:       63
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: Total number of variables............................:       53
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: variables with only lower bounds:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: variables with lower and upper bounds:       10
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: variables with only upper bounds:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: Total number of equality constraints.................:       53
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: 0  0.0000000e+00 1.26e+02 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: 1  0.0000000e+00 2.04e+00 3.36e+00  -1.0 5.77e+02    -  7.71e-01 9.92e-01H  1
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: 2  0.0000000e+00 2.39e-01 9.99e+00  -1.0 7.37e+00    -  9.90e-01 1.00e+00h  1
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: 3  0.0000000e+00 7.45e-09 3.70e+01  -1.0 5.26e-03    -  9.90e-01 1.00e+00h  1
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: Number of Iterations....: 3
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: (scaled)                 (unscaled)
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: Constraint violation....:   2.9103830456733704e-11    7.4505805969238281e-09
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: Overall NLP error.......:   2.9103830456733704e-11    7.4505805969238281e-09
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: Number of objective function evaluations             = 5
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: Number of objective gradient evaluations             = 4
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: Number of equality constraint evaluations            = 5
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: Number of equality constraint Jacobian evaluations   = 4
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: Number of Lagrangian Hessian evaluations             = 3
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.002
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:04 [DEBUG] idaes.solve.fs.M101: 
    2020-04-14 00:02:04 [INFO] idaes.init.fs.M101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:05 [INFO] idaes.init.fs.H101.control_volume: Initialization Complete
    2020-04-14 00:02:05 [INFO] idaes.init.fs.H101: Initialization Step 1 Complete.
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: ******************************************************************************
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: contain the following acknowledgement:
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: ******************************************************************************
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: Number of nonzeros in equality constraint Jacobian...:      124
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: Number of nonzeros in Lagrangian Hessian.............:      112
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: Total number of variables............................:       41
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: variables with only lower bounds:        0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: variables with lower and upper bounds:        9
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: variables with only upper bounds:        0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: Total number of equality constraints.................:       41
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: 0  0.0000000e+00 5.39e+04 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: 1  0.0000000e+00 2.72e+04 1.17e+01  -1.0 9.60e+03    -  6.94e-02 5.28e-01h  1
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: 2  0.0000000e+00 1.28e+04 5.29e+00  -1.0 4.88e+03    -  9.90e-01 5.21e-01h  1
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: 3  0.0000000e+00 5.91e+03 9.23e+03  -1.0 2.41e+03    -  9.90e-01 5.33e-01h  1
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: 4  0.0000000e+00 1.90e+03 4.37e+08  -1.0 1.14e+03    -  9.94e-01 1.00e+00h  1
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: 5  0.0000000e+00 6.52e+01 3.15e+06  -1.0 1.90e+01    -  1.00e+00 5.00e-01h  1
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: 6  0.0000000e+00 3.06e+01 4.77e+06  -1.0 9.46e+00    -  1.00e+00 4.96e-01h  2
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: 7  0.0000000e+00 2.08e+01 3.98e+10  -1.0 4.77e+00    -  1.00e+00 9.94e-01h  1
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: 8  0.0000000e+00 5.81e+00 3.45e+09  -1.0 2.75e-02    -  1.00e+00 5.00e-01h  2
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: 9  0.0000000e+00 4.56e+00 2.76e+10  -1.0 1.38e-02    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: 10  0.0000000e+00 1.09e-04 3.05e+05  -1.0 4.50e-07    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: 11  0.0000000e+00 2.18e-11 4.40e+01  -3.8 3.66e-08    -  1.00e+00 1.00e+00   0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: Number of Iterations....: 11
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: (scaled)                 (unscaled)
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: Constraint violation....:   2.2737367544323206e-13    2.1827872842550274e-11
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: Overall NLP error.......:   2.2737367544323206e-13    2.1827872842550274e-11
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: Number of objective function evaluations             = 17
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: Number of objective gradient evaluations             = 12
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: Number of equality constraint evaluations            = 17
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: Number of equality constraint Jacobian evaluations   = 12
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: Number of Lagrangian Hessian evaluations             = 11
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.004
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.H101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:05 [INFO] idaes.init.fs.H101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:02:05 [INFO] idaes.init.fs.H101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:05 [INFO] idaes.init.fs.R101.control_volume: Initialization Complete
    2020-04-14 00:02:05 [INFO] idaes.init.fs.R101: Initialization Step 1 Complete.
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: ******************************************************************************
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: contain the following acknowledgement:
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: ******************************************************************************
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: Number of nonzeros in equality constraint Jacobian...:       93
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: Number of nonzeros in Lagrangian Hessian.............:       61
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: Total number of variables............................:       39
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: variables with only lower bounds:        0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: variables with lower and upper bounds:       10
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: variables with only upper bounds:        0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: Total number of equality constraints.................:       39
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: 0  0.0000000e+00 3.25e+04 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: 1  0.0000000e+00 3.25e+04 5.13e+00  -1.0 5.50e+04    -  4.66e-01 6.04e-05h 15
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: 2  0.0000000e+00 3.25e+04 1.16e+01  -1.0 5.50e+04    -  6.16e-01 6.04e-05h 15
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: 3  0.0000000e+00 3.25e+04 1.85e+01  -1.0 5.50e+04    -  6.71e-01 6.04e-05h 15
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: 4  0.0000000e+00 3.25e+04 2.84e+01  -1.0 5.50e+04    -  9.90e-01 6.04e-05h 15
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: 5  0.0000000e+00 3.25e+04 3.40e+01  -1.0 5.50e+04    -  5.58e-01 6.04e-05h 15
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: 6  0.0000000e+00 3.25e+04 4.39e+01  -1.0 5.50e+04    -  9.91e-01 6.04e-05h 15
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: 7  0.0000000e+00 3.25e+04 4.91e+01  -1.0 5.50e+04    -  5.13e-01 6.04e-05h 15
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: 8  0.0000000e+00 3.25e+04 5.91e+01  -1.0 5.50e+04    -  1.00e+00 6.04e-05h 15
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: 9  0.0000000e+00 3.25e+04 6.42e+01  -1.0 5.49e+04    -  5.12e-01 6.04e-05h 15
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: 10  0.0000000e+00 3.25e+04 7.42e+01  -1.0 5.49e+04    -  1.00e+00 6.04e-05h 15
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: 11  0.0000000e+00 1.27e+07 9.78e+01  -1.0 5.49e+04    -  5.12e-01 9.90e-01w  1
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: 12  0.0000000e+00 1.74e+05 1.07e+02  -1.0 8.06e+03    -  1.00e+00 9.90e-01w  1
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: 13  0.0000000e+00 4.06e+04 2.63e+00  -1.0 2.26e+04    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: 14  0.0000000e+00 3.61e-05 5.36e-02  -3.8 4.06e+04    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: Number of Iterations....: 14
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: (scaled)                 (unscaled)
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: Constraint violation....:   1.2500647559077285e-08    3.6127865314483643e-05
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: Overall NLP error.......:   1.2500647559077285e-08    3.6127865314483643e-05
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: Number of objective function evaluations             = 185
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: Number of objective gradient evaluations             = 15
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: Number of equality constraint evaluations            = 195
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: Number of equality constraint Jacobian evaluations   = 15
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: Number of Lagrangian Hessian evaluations             = 14
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.027
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: Total CPU secs in NLP function evaluations           =      0.001
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.R101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:05 [INFO] idaes.init.fs.R101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:02:05 [INFO] idaes.init.fs.R101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:05 [INFO] idaes.init.fs.F101.control_volume: Initialization Complete
    2020-04-14 00:02:05 [INFO] idaes.init.fs.F101: Initialization Step 1 Complete.
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: ******************************************************************************
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: contain the following acknowledgement:
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: ******************************************************************************
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: Number of nonzeros in equality constraint Jacobian...:      124
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: Number of nonzeros in Lagrangian Hessian.............:      112
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: Total number of variables............................:       41
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: variables with only lower bounds:        0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: variables with lower and upper bounds:        9
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: variables with only upper bounds:        0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: Total number of equality constraints.................:       41
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: 0  0.0000000e+00 5.12e+04 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: 1  0.0000000e+00 3.74e+04 5.83e+00  -1.0 3.26e+03    -  1.35e-01 6.68e-01f  1
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: 2  0.0000000e+00 9.79e+03 2.83e+01  -1.0 2.31e+03    -  8.78e-01 6.18e-01h  1
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: 3  0.0000000e+00 1.77e+02 2.07e+01  -1.0 7.29e+02    -  8.20e-01 1.00e+00h  1
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: 4  0.0000000e+00 1.44e+00 2.51e+00  -1.0 1.53e+01    -  9.90e-01 9.92e-01h  1
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: 5  0.0000000e+00 2.07e-06 6.41e+02  -1.0 1.26e-01    -  9.92e-01 1.00e+00h  1
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: Number of Iterations....: 5
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: (scaled)                 (unscaled)
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: Constraint violation....:   4.2890462401563167e-11    2.0655534171964973e-06
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: Overall NLP error.......:   4.2890462401563167e-11    2.0655534171964973e-06
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: Number of objective function evaluations             = 6
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: Number of objective gradient evaluations             = 6
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: Number of equality constraint evaluations            = 6
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: Number of equality constraint Jacobian evaluations   = 6
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: Number of Lagrangian Hessian evaluations             = 5
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.012
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.F101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:05 [INFO] idaes.init.fs.F101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:02:05 [INFO] idaes.init.fs.F101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.S101: Ipopt 3.13.2:
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.S101: ******************************************************************************
    2020-04-14 00:02:05 [DEBUG] idaes.solve.fs.S101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: contain the following acknowledgement:
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: ******************************************************************************
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: Number of nonzeros in equality constraint Jacobian...:       29
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: Number of nonzeros in Lagrangian Hessian.............:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: Total number of variables............................:       21
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: variables with only lower bounds:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: variables with lower and upper bounds:       20
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: variables with only upper bounds:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: Total number of equality constraints.................:       21
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: 0  0.0000000e+00 1.00e-02 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: 1  0.0000000e+00 1.00e-04 2.20e-06  -1.0 1.00e-02    -  9.90e-01 9.90e-01h  1
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: 2  0.0000000e+00 9.98e-07 2.00e-02  -1.0 1.00e-04    -  9.90e-01 9.90e-01h  1
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: 3  0.0000000e+00 8.00e-09 2.01e+02  -1.0 9.98e-07    -  9.94e-01 9.92e-01h  1
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: Number of Iterations....: 3
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: (scaled)                 (unscaled)
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: Constraint violation....:   7.9999999999999045e-09    7.9999999999999045e-09
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: Overall NLP error.......:   7.9999999999999045e-09    7.9999999999999045e-09
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: Number of objective function evaluations             = 4
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: Number of objective gradient evaluations             = 4
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: Number of equality constraint evaluations            = 4
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: Number of equality constraint Jacobian evaluations   = 4
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: Number of Lagrangian Hessian evaluations             = 3
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.000
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.S101: 
    2020-04-14 00:02:06 [INFO] idaes.init.fs.S101: Initialization Step 2 Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:06 [INFO] idaes.init.fs.C101.control_volume: Initialization Complete
    2020-04-14 00:02:06 [INFO] idaes.init.fs.C101: Initialization Step 1 Complete.
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: ******************************************************************************
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: contain the following acknowledgement:
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: ******************************************************************************
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: Number of nonzeros in equality constraint Jacobian...:       74
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: Number of nonzeros in Lagrangian Hessian.............:       61
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: Total number of variables............................:       32
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: variables with only lower bounds:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: variables with lower and upper bounds:        9
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: variables with only upper bounds:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: Total number of equality constraints.................:       32
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: 0  0.0000000e+00 1.00e-02 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: 1  0.0000000e+00 1.00e-04 1.60e-03  -1.0 1.00e-02    -  9.90e-01 9.90e-01h  1
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: 2  0.0000000e+00 9.84e-07 9.84e+00  -1.0 1.00e-04    -  1.00e+00 9.90e-01h  1
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: Number of Iterations....: 2
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: (scaled)                 (unscaled)
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: Constraint violation....:   9.8400000480001186e-07    9.8400000480001186e-07
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: Overall NLP error.......:   9.8400000480001186e-07    9.8400000480001186e-07
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: Number of objective function evaluations             = 3
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: Number of objective gradient evaluations             = 3
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: Number of equality constraint evaluations            = 3
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: Number of equality constraint Jacobian evaluations   = 3
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: Number of Lagrangian Hessian evaluations             = 2
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.001
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.C101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:06 [INFO] idaes.init.fs.C101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:02:06 [INFO] idaes.init.fs.C101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: Ipopt 3.13.2:
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: ******************************************************************************
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: contain the following acknowledgement:
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: ******************************************************************************
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: Number of nonzeros in equality constraint Jacobian...:      117
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: Number of nonzeros in Lagrangian Hessian.............:       63
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: Total number of variables............................:       53
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: variables with only lower bounds:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: variables with lower and upper bounds:       10
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: variables with only upper bounds:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: Total number of equality constraints.................:       53
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: 0  0.0000000e+00 1.26e+02 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: 1  0.0000000e+00 2.20e+00 3.38e+00  -1.0 5.80e+02    -  7.70e-01 9.92e-01H  1
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: 2  0.0000000e+00 2.41e-01 9.99e+00  -1.0 7.41e+00    -  9.90e-01 1.00e+00h  1
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: 3  0.0000000e+00 7.45e-09 3.70e+01  -1.0 5.28e-03    -  9.90e-01 1.00e+00h  1
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: Number of Iterations....: 3
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: (scaled)                 (unscaled)
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: Constraint violation....:   2.9103830456733704e-11    7.4505805969238281e-09
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: Overall NLP error.......:   2.9103830456733704e-11    7.4505805969238281e-09
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: Number of objective function evaluations             = 5
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: Number of objective gradient evaluations             = 4
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: Number of equality constraint evaluations            = 5
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: Number of equality constraint Jacobian evaluations   = 4
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: Number of Lagrangian Hessian evaluations             = 3
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.001
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.M101: 
    2020-04-14 00:02:06 [INFO] idaes.init.fs.M101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:06 [INFO] idaes.init.fs.H101.control_volume: Initialization Complete
    2020-04-14 00:02:06 [INFO] idaes.init.fs.H101: Initialization Step 1 Complete.
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: ******************************************************************************
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: contain the following acknowledgement:
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: ******************************************************************************
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: Number of nonzeros in equality constraint Jacobian...:      124
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: Number of nonzeros in Lagrangian Hessian.............:      112
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: Total number of variables............................:       41
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: variables with only lower bounds:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: variables with lower and upper bounds:        9
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: variables with only upper bounds:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: Total number of equality constraints.................:       41
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: 0  0.0000000e+00 5.39e+04 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: 1  0.0000000e+00 2.72e+04 1.17e+01  -1.0 9.61e+03    -  6.95e-02 5.28e-01h  1
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: 2  0.0000000e+00 1.28e+04 5.29e+00  -1.0 4.88e+03    -  9.90e-01 5.21e-01h  1
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: 3  0.0000000e+00 5.91e+03 9.23e+03  -1.0 2.41e+03    -  9.90e-01 5.33e-01h  1
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: 4  0.0000000e+00 1.90e+03 4.37e+08  -1.0 1.14e+03    -  9.94e-01 1.00e+00h  1
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: 5  0.0000000e+00 6.53e+01 3.15e+06  -1.0 1.90e+01    -  1.00e+00 5.00e-01h  1
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: 6  0.0000000e+00 3.06e+01 4.77e+06  -1.0 9.46e+00    -  1.00e+00 4.96e-01h  2
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: 7  0.0000000e+00 2.08e+01 3.97e+10  -1.0 4.76e+00    -  1.00e+00 9.94e-01h  1
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: 8  0.0000000e+00 5.82e+00 3.45e+09  -1.0 2.75e-02    -  1.00e+00 5.00e-01h  2
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: 9  0.0000000e+00 4.58e+00 2.76e+10  -1.0 1.38e-02    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: 10  0.0000000e+00 1.09e-04 3.07e+05  -1.0 4.51e-07    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: 11  0.0000000e+00 7.28e-12 4.32e+01  -3.8 2.89e-08    -  1.00e+00 1.00e+00   0
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: Number of Iterations....: 11
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: (scaled)                 (unscaled)
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: Constraint violation....:   2.2737367544323206e-13    7.2759576141834259e-12
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: Overall NLP error.......:   2.2737367544323206e-13    7.2759576141834259e-12
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: Number of objective function evaluations             = 17
    2020-04-14 00:02:06 [DEBUG] idaes.solve.fs.H101: Number of objective gradient evaluations             = 12
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.H101: Number of equality constraint evaluations            = 17
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.H101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.H101: Number of equality constraint Jacobian evaluations   = 12
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.H101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.H101: Number of Lagrangian Hessian evaluations             = 11
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.H101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.003
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.H101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.H101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:07 [INFO] idaes.init.fs.H101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:02:07 [INFO] idaes.init.fs.H101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:07 [INFO] idaes.init.fs.R101.control_volume: Initialization Complete
    2020-04-14 00:02:07 [INFO] idaes.init.fs.R101: Initialization Step 1 Complete.
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: ******************************************************************************
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: contain the following acknowledgement:
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: ******************************************************************************
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: Number of nonzeros in equality constraint Jacobian...:       93
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: Number of nonzeros in Lagrangian Hessian.............:       61
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: Total number of variables............................:       39
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: variables with only lower bounds:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: variables with lower and upper bounds:       10
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: variables with only upper bounds:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: Total number of equality constraints.................:       39
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: 0  0.0000000e+00 3.25e+04 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: 1  0.0000000e+00 3.25e+04 5.13e+00  -1.0 5.50e+04    -  4.66e-01 6.04e-05h 15
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: 2  0.0000000e+00 3.25e+04 1.16e+01  -1.0 5.50e+04    -  6.16e-01 6.04e-05h 15
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: 3  0.0000000e+00 3.25e+04 1.85e+01  -1.0 5.50e+04    -  6.70e-01 6.04e-05h 15
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: 4  0.0000000e+00 3.25e+04 2.84e+01  -1.0 5.50e+04    -  9.90e-01 6.04e-05h 15
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: 5  0.0000000e+00 3.25e+04 3.40e+01  -1.0 5.50e+04    -  5.57e-01 6.04e-05h 15
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: 6  0.0000000e+00 3.25e+04 4.39e+01  -1.0 5.50e+04    -  9.91e-01 6.04e-05h 15
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: 7  0.0000000e+00 3.25e+04 4.90e+01  -1.0 5.50e+04    -  5.12e-01 6.04e-05h 15
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: 8  0.0000000e+00 3.25e+04 5.90e+01  -1.0 5.50e+04    -  1.00e+00 6.04e-05h 15
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: 9  0.0000000e+00 3.25e+04 6.41e+01  -1.0 5.50e+04    -  5.11e-01 6.04e-05h 15
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: 10  0.0000000e+00 3.25e+04 7.41e+01  -1.0 5.50e+04    -  1.00e+00 6.04e-05h 15
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: 11  0.0000000e+00 1.27e+07 9.86e+01  -1.0 5.50e+04    -  5.11e-01 9.90e-01w  1
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: 12  0.0000000e+00 1.75e+05 1.08e+02  -1.0 8.07e+03    -  1.00e+00 9.90e-01w  1
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: 13  0.0000000e+00 4.06e+04 2.68e+00  -1.0 2.26e+04    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: 14  0.0000000e+00 3.61e-05 5.37e-02  -3.8 4.06e+04    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: Number of Iterations....: 14
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: (scaled)                 (unscaled)
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: Constraint violation....:   1.2517822710291955e-08    3.6142766475677490e-05
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: Overall NLP error.......:   1.2517822710291955e-08    3.6142766475677490e-05
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: Number of objective function evaluations             = 185
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: Number of objective gradient evaluations             = 15
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: Number of equality constraint evaluations            = 195
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: Number of equality constraint Jacobian evaluations   = 15
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: Number of Lagrangian Hessian evaluations             = 14
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.008
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: Total CPU secs in NLP function evaluations           =      0.002
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.R101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:07 [INFO] idaes.init.fs.R101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:02:07 [INFO] idaes.init.fs.R101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:07 [INFO] idaes.init.fs.F101.control_volume: Initialization Complete
    2020-04-14 00:02:07 [INFO] idaes.init.fs.F101: Initialization Step 1 Complete.
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: ******************************************************************************
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: contain the following acknowledgement:
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: ******************************************************************************
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: Number of nonzeros in equality constraint Jacobian...:      124
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: Number of nonzeros in Lagrangian Hessian.............:      112
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: Total number of variables............................:       41
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: variables with only lower bounds:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: variables with lower and upper bounds:        9
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: variables with only upper bounds:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: Total number of equality constraints.................:       41
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: 0  0.0000000e+00 5.12e+04 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: 1  0.0000000e+00 3.74e+04 5.83e+00  -1.0 3.26e+03    -  1.35e-01 6.68e-01f  1
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: 2  0.0000000e+00 9.80e+03 2.83e+01  -1.0 2.31e+03    -  8.78e-01 6.17e-01h  1
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: 3  0.0000000e+00 1.77e+02 2.07e+01  -1.0 7.30e+02    -  8.20e-01 1.00e+00h  1
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: 4  0.0000000e+00 1.44e+00 2.51e+00  -1.0 1.53e+01    -  9.90e-01 9.92e-01h  1
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: 5  0.0000000e+00 2.07e-06 6.41e+02  -1.0 1.26e-01    -  9.92e-01 1.00e+00h  1
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: Number of Iterations....: 5
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: (scaled)                 (unscaled)
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: Constraint violation....:   4.2902110935282590e-11    2.0713669073302299e-06
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: Overall NLP error.......:   4.2902110935282590e-11    2.0713669073302299e-06
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: Number of objective function evaluations             = 6
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: Number of objective gradient evaluations             = 6
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: Number of equality constraint evaluations            = 6
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: Number of equality constraint Jacobian evaluations   = 6
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: Number of Lagrangian Hessian evaluations             = 5
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.017
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.F101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:07 [INFO] idaes.init.fs.F101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:02:07 [INFO] idaes.init.fs.F101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: Ipopt 3.13.2:
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: ******************************************************************************
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: contain the following acknowledgement:
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: ******************************************************************************
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: Number of nonzeros in equality constraint Jacobian...:       29
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: Number of nonzeros in Lagrangian Hessian.............:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: Total number of variables............................:       21
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: variables with only lower bounds:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: variables with lower and upper bounds:       20
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: variables with only upper bounds:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: Total number of equality constraints.................:       21
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: 0  0.0000000e+00 1.00e-02 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: 1  0.0000000e+00 1.00e-04 2.20e-06  -1.0 1.00e-02    -  9.90e-01 9.90e-01h  1
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: 2  0.0000000e+00 9.98e-07 2.00e-02  -1.0 1.00e-04    -  9.90e-01 9.90e-01h  1
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: 3  0.0000000e+00 8.00e-09 2.01e+02  -1.0 9.98e-07    -  9.94e-01 9.92e-01h  1
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: Number of Iterations....: 3
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: (scaled)                 (unscaled)
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: Constraint violation....:   7.9999999999999045e-09    7.9999999999999045e-09
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: Overall NLP error.......:   7.9999999999999045e-09    7.9999999999999045e-09
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: Number of objective function evaluations             = 4
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: Number of objective gradient evaluations             = 4
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: Number of equality constraint evaluations            = 4
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: Number of equality constraint Jacobian evaluations   = 4
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: Number of Lagrangian Hessian evaluations             = 3
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.001
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.S101: 
    2020-04-14 00:02:07 [INFO] idaes.init.fs.S101: Initialization Step 2 Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:07 [INFO] idaes.init.fs.C101.control_volume: Initialization Complete
    2020-04-14 00:02:07 [INFO] idaes.init.fs.C101: Initialization Step 1 Complete.
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: ******************************************************************************
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: contain the following acknowledgement:
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: ******************************************************************************
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: Number of nonzeros in equality constraint Jacobian...:       74
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: Number of nonzeros in Lagrangian Hessian.............:       61
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: Total number of variables............................:       32
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: variables with only lower bounds:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: variables with lower and upper bounds:        9
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: variables with only upper bounds:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: Total number of equality constraints.................:       32
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:07 [DEBUG] idaes.solve.fs.C101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.C101: 0  0.0000000e+00 1.00e-02 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.C101: 1  0.0000000e+00 1.00e-04 1.60e-03  -1.0 1.00e-02    -  9.90e-01 9.90e-01h  1
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.C101: 2  0.0000000e+00 9.84e-07 9.84e+00  -1.0 1.00e-04    -  1.00e+00 9.90e-01h  1
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.C101: Number of Iterations....: 2
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.C101: (scaled)                 (unscaled)
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.C101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.C101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.C101: Constraint violation....:   9.8400000480001186e-07    9.8400000480001186e-07
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.C101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.C101: Overall NLP error.......:   9.8400000480001186e-07    9.8400000480001186e-07
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.C101: Number of objective function evaluations             = 3
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.C101: Number of objective gradient evaluations             = 3
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.C101: Number of equality constraint evaluations            = 3
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.C101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.C101: Number of equality constraint Jacobian evaluations   = 3
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.C101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.C101: Number of Lagrangian Hessian evaluations             = 2
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.C101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.001
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.C101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.C101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:08 [INFO] idaes.init.fs.C101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:02:08 [INFO] idaes.init.fs.C101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: Ipopt 3.13.2:
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: ******************************************************************************
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: contain the following acknowledgement:
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: ******************************************************************************
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: Number of nonzeros in equality constraint Jacobian...:      117
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: Number of nonzeros in Lagrangian Hessian.............:       63
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: Total number of variables............................:       53
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: variables with only lower bounds:        0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: variables with lower and upper bounds:       10
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: variables with only upper bounds:        0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: Total number of equality constraints.................:       53
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: 0  0.0000000e+00 1.26e+02 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: 1  0.0000000e+00 2.17e+00 3.38e+00  -1.0 5.79e+02    -  7.70e-01 9.92e-01H  1
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: 2  0.0000000e+00 2.41e-01 9.99e+00  -1.0 7.40e+00    -  9.90e-01 1.00e+00h  1
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: 3  0.0000000e+00 7.45e-09 3.70e+01  -1.0 5.28e-03    -  9.90e-01 1.00e+00h  1
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: Number of Iterations....: 3
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: (scaled)                 (unscaled)
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: Constraint violation....:   2.9103830456733704e-11    7.4505805969238281e-09
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: Overall NLP error.......:   2.9103830456733704e-11    7.4505805969238281e-09
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: Number of objective function evaluations             = 5
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: Number of objective gradient evaluations             = 4
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: Number of equality constraint evaluations            = 5
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: Number of equality constraint Jacobian evaluations   = 4
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: Number of Lagrangian Hessian evaluations             = 3
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.001
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.M101: 
    2020-04-14 00:02:08 [INFO] idaes.init.fs.M101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:08 [INFO] idaes.init.fs.H101.control_volume: Initialization Complete
    2020-04-14 00:02:08 [INFO] idaes.init.fs.H101: Initialization Step 1 Complete.
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: ******************************************************************************
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: contain the following acknowledgement:
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: ******************************************************************************
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: Number of nonzeros in equality constraint Jacobian...:      124
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: Number of nonzeros in Lagrangian Hessian.............:      112
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: Total number of variables............................:       41
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: variables with only lower bounds:        0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: variables with lower and upper bounds:        9
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: variables with only upper bounds:        0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: Total number of equality constraints.................:       41
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: 0  0.0000000e+00 5.40e+04 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: 1  0.0000000e+00 2.72e+04 1.17e+01  -1.0 9.62e+03    -  6.95e-02 5.28e-01h  1
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: 2  0.0000000e+00 1.28e+04 5.29e+00  -1.0 4.89e+03    -  9.90e-01 5.21e-01h  1
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: 3  0.0000000e+00 5.92e+03 9.23e+03  -1.0 2.42e+03    -  9.90e-01 5.33e-01h  1
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: 4  0.0000000e+00 1.90e+03 4.37e+08  -1.0 1.14e+03    -  9.94e-01 1.00e+00h  1
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: 5  0.0000000e+00 6.55e+01 3.16e+06  -1.0 1.90e+01    -  1.00e+00 5.00e-01h  1
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: 6  0.0000000e+00 3.07e+01 4.77e+06  -1.0 9.46e+00    -  1.00e+00 4.96e-01h  2
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: 7  0.0000000e+00 2.09e+01 3.97e+10  -1.0 4.76e+00    -  1.00e+00 9.94e-01h  1
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: 8  0.0000000e+00 5.84e+00 3.45e+09  -1.0 2.76e-02    -  1.00e+00 5.00e-01h  2
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: 9  0.0000000e+00 4.60e+00 2.77e+10  -1.0 1.38e-02    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: 10  0.0000000e+00 1.10e-04 3.09e+05  -1.0 4.54e-07    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: 11  0.0000000e+00 7.28e-12 6.04e+01  -3.8 4.05e-08    -  1.00e+00 1.00e+00   0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: Number of Iterations....: 11
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: (scaled)                 (unscaled)
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: Constraint violation....:   2.2737367544323206e-13    7.2759576141834259e-12
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: Overall NLP error.......:   2.2737367544323206e-13    7.2759576141834259e-12
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: Number of objective function evaluations             = 17
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: Number of objective gradient evaluations             = 12
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: Number of equality constraint evaluations            = 17
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: Number of equality constraint Jacobian evaluations   = 12
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: Number of Lagrangian Hessian evaluations             = 11
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.002
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.H101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:08 [INFO] idaes.init.fs.H101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:02:08 [INFO] idaes.init.fs.H101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:08 [INFO] idaes.init.fs.R101.control_volume: Initialization Complete
    2020-04-14 00:02:08 [INFO] idaes.init.fs.R101: Initialization Step 1 Complete.
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: ******************************************************************************
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: contain the following acknowledgement:
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: ******************************************************************************
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: Number of nonzeros in equality constraint Jacobian...:       93
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: Number of nonzeros in Lagrangian Hessian.............:       61
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: Total number of variables............................:       39
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: variables with only lower bounds:        0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: variables with lower and upper bounds:       10
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: variables with only upper bounds:        0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: Total number of equality constraints.................:       39
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: 0  0.0000000e+00 3.25e+04 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: 1  0.0000000e+00 3.25e+04 5.12e+00  -1.0 5.50e+04    -  4.66e-01 6.04e-05h 15
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: 2  0.0000000e+00 3.25e+04 1.16e+01  -1.0 5.50e+04    -  6.15e-01 6.04e-05h 15
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: 3  0.0000000e+00 3.25e+04 1.84e+01  -1.0 5.50e+04    -  6.70e-01 6.04e-05h 15
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: 4  0.0000000e+00 3.25e+04 2.84e+01  -1.0 5.50e+04    -  9.90e-01 6.04e-05h 15
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: 5  0.0000000e+00 3.25e+04 3.40e+01  -1.0 5.50e+04    -  5.56e-01 6.04e-05h 15
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: 6  0.0000000e+00 3.25e+04 4.39e+01  -1.0 5.50e+04    -  9.91e-01 6.04e-05h 15
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: 7  0.0000000e+00 3.25e+04 4.90e+01  -1.0 5.50e+04    -  5.12e-01 6.04e-05h 15
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: 8  0.0000000e+00 3.25e+04 5.90e+01  -1.0 5.50e+04    -  1.00e+00 6.04e-05h 15
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: 9  0.0000000e+00 3.25e+04 6.41e+01  -1.0 5.50e+04    -  5.10e-01 6.04e-05h 15
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: 10  0.0000000e+00 3.25e+04 7.41e+01  -1.0 5.50e+04    -  1.00e+00 6.04e-05h 15
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: 11  0.0000000e+00 1.27e+07 9.94e+01  -1.0 5.50e+04    -  5.10e-01 9.90e-01w  1
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: 12  0.0000000e+00 1.75e+05 1.08e+02  -1.0 8.09e+03    -  1.00e+00 9.90e-01w  1
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: 13  0.0000000e+00 4.06e+04 2.73e+00  -1.0 2.26e+04    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: 14  0.0000000e+00 3.63e-05 5.38e-02  -3.8 4.06e+04    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: Number of Iterations....: 14
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: (scaled)                 (unscaled)
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: Constraint violation....:   1.2552173012721294e-08    3.6276876926422119e-05
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: Overall NLP error.......:   1.2552173012721294e-08    3.6276876926422119e-05
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: Number of objective function evaluations             = 185
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: Number of objective gradient evaluations             = 15
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: Number of equality constraint evaluations            = 195
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: Number of equality constraint Jacobian evaluations   = 15
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: Number of Lagrangian Hessian evaluations             = 14
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.032
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:08 [DEBUG] idaes.solve.fs.R101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:09 [INFO] idaes.init.fs.R101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:02:09 [INFO] idaes.init.fs.R101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:09 [INFO] idaes.init.fs.F101.control_volume: Initialization Complete
    2020-04-14 00:02:09 [INFO] idaes.init.fs.F101: Initialization Step 1 Complete.
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: ******************************************************************************
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: contain the following acknowledgement:
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: ******************************************************************************
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: Number of nonzeros in equality constraint Jacobian...:      124
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: Number of nonzeros in Lagrangian Hessian.............:      112
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: Total number of variables............................:       41
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: variables with only lower bounds:        0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: variables with lower and upper bounds:        9
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: variables with only upper bounds:        0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: Total number of equality constraints.................:       41
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: 0  0.0000000e+00 5.13e+04 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: 1  0.0000000e+00 3.75e+04 5.85e+00  -1.0 3.26e+03    -  1.34e-01 6.67e-01f  1
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: 2  0.0000000e+00 9.81e+03 2.83e+01  -1.0 2.31e+03    -  8.78e-01 6.17e-01h  1
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: 3  0.0000000e+00 1.77e+02 2.07e+01  -1.0 7.31e+02    -  8.20e-01 1.00e+00h  1
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: 4  0.0000000e+00 1.44e+00 2.51e+00  -1.0 1.54e+01    -  9.90e-01 9.92e-01h  1
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: 5  0.0000000e+00 2.11e-06 6.40e+02  -1.0 1.27e-01    -  9.92e-01 1.00e+00h  1
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: Number of Iterations....: 5
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: (scaled)                 (unscaled)
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: Constraint violation....:   4.2680665616048879e-11    2.1109808585606515e-06
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: Overall NLP error.......:   4.2680665616048879e-11    2.1109808585606515e-06
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: Number of objective function evaluations             = 6
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: Number of objective gradient evaluations             = 6
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: Number of equality constraint evaluations            = 6
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: Number of equality constraint Jacobian evaluations   = 6
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: Number of Lagrangian Hessian evaluations             = 5
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.002
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.F101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:09 [INFO] idaes.init.fs.F101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:02:09 [INFO] idaes.init.fs.F101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: Ipopt 3.13.2:
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: ******************************************************************************
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: contain the following acknowledgement:
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: ******************************************************************************
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: Number of nonzeros in equality constraint Jacobian...:       29
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: Number of nonzeros in Lagrangian Hessian.............:        0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: Total number of variables............................:       21
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: variables with only lower bounds:        0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: variables with lower and upper bounds:       20
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: variables with only upper bounds:        0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: Total number of equality constraints.................:       21
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: 0  0.0000000e+00 1.00e-02 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: 1  0.0000000e+00 1.00e-04 2.20e-06  -1.0 1.00e-02    -  9.90e-01 9.90e-01h  1
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: 2  0.0000000e+00 9.98e-07 2.00e-02  -1.0 1.00e-04    -  9.90e-01 9.90e-01h  1
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: 3  0.0000000e+00 8.00e-09 2.01e+02  -1.0 9.98e-07    -  9.94e-01 9.92e-01h  1
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: Number of Iterations....: 3
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: (scaled)                 (unscaled)
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: Constraint violation....:   7.9999999999999045e-09    7.9999999999999045e-09
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: Overall NLP error.......:   7.9999999999999045e-09    7.9999999999999045e-09
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: Number of objective function evaluations             = 4
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: Number of objective gradient evaluations             = 4
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: Number of equality constraint evaluations            = 4
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: Number of equality constraint Jacobian evaluations   = 4
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: Number of Lagrangian Hessian evaluations             = 3
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.001
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.S101: 
    2020-04-14 00:02:09 [INFO] idaes.init.fs.S101: Initialization Step 2 Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:09 [INFO] idaes.init.fs.C101.control_volume: Initialization Complete
    2020-04-14 00:02:09 [INFO] idaes.init.fs.C101: Initialization Step 1 Complete.
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: ******************************************************************************
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: contain the following acknowledgement:
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: ******************************************************************************
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: Number of nonzeros in equality constraint Jacobian...:       74
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: Number of nonzeros in Lagrangian Hessian.............:       61
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: Total number of variables............................:       32
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: variables with only lower bounds:        0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: variables with lower and upper bounds:        9
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: variables with only upper bounds:        0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: Total number of equality constraints.................:       32
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: 0  0.0000000e+00 1.00e-02 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: 1  0.0000000e+00 1.00e-04 1.60e-03  -1.0 1.00e-02    -  9.90e-01 9.90e-01h  1
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: 2  0.0000000e+00 9.84e-07 9.84e+00  -1.0 1.00e-04    -  1.00e+00 9.90e-01h  1
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: Number of Iterations....: 2
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: (scaled)                 (unscaled)
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: Constraint violation....:   9.8400000480001186e-07    9.8400000480001186e-07
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: Overall NLP error.......:   9.8400000480001186e-07    9.8400000480001186e-07
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: Number of objective function evaluations             = 3
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: Number of objective gradient evaluations             = 3
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: Number of equality constraint evaluations            = 3
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: Number of equality constraint Jacobian evaluations   = 3
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: Number of Lagrangian Hessian evaluations             = 2
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.001
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:09 [DEBUG] idaes.solve.fs.C101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:09 [INFO] idaes.init.fs.C101: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:02:09 [INFO] idaes.init.fs.C101: Initialization Complete: optimal - Optimal Solution Found
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: Ipopt 3.13.2:
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: ******************************************************************************
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: contain the following acknowledgement:
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: ******************************************************************************
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: Number of nonzeros in equality constraint Jacobian...:      117
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: Number of nonzeros in Lagrangian Hessian.............:       63
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: Total number of variables............................:       53
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: variables with only lower bounds:        0
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: variables with lower and upper bounds:       10
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: variables with only upper bounds:        0
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: Total number of equality constraints.................:       53
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: Total number of inequality constraints...............:        0
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: 0  0.0000000e+00 1.26e+02 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: 1  0.0000000e+00 2.12e+00 3.37e+00  -1.0 5.78e+02    -  7.71e-01 9.92e-01H  1
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: 2  0.0000000e+00 2.40e-01 9.99e+00  -1.0 7.39e+00    -  9.90e-01 1.00e+00h  1
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: 3  0.0000000e+00 7.45e-09 3.70e+01  -1.0 5.27e-03    -  9.90e-01 1.00e+00h  1
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: Number of Iterations....: 3
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: (scaled)                 (unscaled)
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: Constraint violation....:   2.9103830456733704e-11    7.4505805969238281e-09
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: Overall NLP error.......:   2.9103830456733704e-11    7.4505805969238281e-09
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: Number of objective function evaluations             = 5
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: Number of objective gradient evaluations             = 4
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: Number of equality constraint evaluations            = 5
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: Number of equality constraint Jacobian evaluations   = 4
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: Number of Lagrangian Hessian evaluations             = 3
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: Total CPU secs in IPOPT (w/o function evaluations)   =      0.001
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: EXIT: Optimal Solution Found.
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.M101: 
    2020-04-14 00:02:10 [INFO] idaes.init.fs.M101: Initialization Complete: optimal - Optimal Solution Found
    WARNING: Wegstein failed to converge in 5 iterations
    2020-04-14 00:02:10 [INFO] idaes.init.fs.F102.control_volume: Initialization Complete
    2020-04-14 00:02:10 [INFO] idaes.init.fs.F102: Initialization Step 1 Complete.
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: Ipopt 3.13.2: tol=1e-06
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: ******************************************************************************
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: This program contains Ipopt, a library for large-scale nonlinear optimization.
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: Ipopt is released as open source code under the Eclipse Public License (EPL).
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: For more information visit http://projects.coin-or.org/Ipopt
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: This version of Ipopt was compiled from source code available at
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: https://github.com/IDAES/Ipopt as part of the Institute for the Design of
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: Advanced Energy Systems Process Systems Engineering Framework (IDAES PSE
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: Framework) Copyright (c) 2018-2019. See https://github.com/IDAES/idaes-pse.
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: This version of Ipopt was compiled using HSL, a collection of Fortran codes
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: for large-scale scientific computation.  All technical papers, sales and
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: publicity material resulting from use of the HSL codes within IPOPT must
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: contain the following acknowledgement:
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: HSL, a collection of Fortran codes for large-scale scientific
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: computation. See http://www.hsl.rl.ac.uk.
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: ******************************************************************************
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: This is Ipopt version 3.13.2, running with linear solver ma27.
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: Number of nonzeros in equality constraint Jacobian...:      124
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: Number of nonzeros in inequality constraint Jacobian.:        0
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: Number of nonzeros in Lagrangian Hessian.............:      112
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: Total number of variables............................:       41
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: variables with only lower bounds:        0
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: variables with lower and upper bounds:        9
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: variables with only upper bounds:        0
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: Total number of equality constraints.................:       41
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: Total number of inequality constraints...............:        0
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: inequality constraints with only lower bounds:        0
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: inequality constraints with lower and upper bounds:        0
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: inequality constraints with only upper bounds:        0
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: 0  0.0000000e+00 6.83e+04 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: 1  0.0000000e+00 2.82e+04 5.69e+00  -1.0 2.00e+05    -  1.98e-01 8.56e-01h  1
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: 2  0.0000000e+00 8.44e+03 1.15e+02  -1.0 2.88e+04    -  8.76e-01 9.57e-01h  1
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: 3  0.0000000e+00 1.46e+03 8.84e+04  -1.0 1.94e+03    -  4.14e-01 9.65e-01h  1
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: 4  0.0000000e+00 1.42e+03 1.94e+06  -1.0 9.53e+02    -  9.91e-01 1.00e+00h  1
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: 5  0.0000000e+00 9.65e-01 8.14e+04  -1.0 1.42e+01    -  9.91e-01 1.00e+00h  1
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: 6  0.0000000e+00 2.25e-04 9.29e+01  -1.0 3.62e-01    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: 7  0.0000000e+00 1.46e-11 1.72e-04  -2.5 8.18e-06    -  1.00e+00 1.00e+00h  1
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: Number of Iterations....: 7
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: (scaled)                 (unscaled)
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: Constraint violation....:   2.2737367544323206e-13    1.4551915228366852e-11
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: Overall NLP error.......:   2.2737367544323206e-13    1.4551915228366852e-11
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: Number of objective function evaluations             = 8
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: Number of objective gradient evaluations             = 8
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: Number of equality constraint evaluations            = 8
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: Number of inequality constraint evaluations          = 0
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: Number of equality constraint Jacobian evaluations   = 8
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: Number of inequality constraint Jacobian evaluations = 0
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: Number of Lagrangian Hessian evaluations             = 7
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: Total CPU secs in IPOPT (w/o function evaluations)   =      0.014
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: Total CPU secs in NLP function evaluations           =      0.000
    2020-04-14 00:02:10 [DEBUG] idaes.solve.fs.F102: EXIT: Optimal Solution Found.
    2020-04-14 00:02:10 [INFO] idaes.init.fs.F102: Initialization Step 2 optimal - Optimal Solution Found.
    2020-04-14 00:02:10 [INFO] idaes.init.fs.F102: Initialization Complete: optimal - Optimal Solution Found


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
    results = solver.solve(m, tee=False)

.. code:: ipython3

    # For testing purposes
    from pyomo.environ import TerminationCondition
    assert results.solver.termination_condition == TerminationCondition.optimal

Analyze the results of the square problem
-----------------------------------------

What is the total operating cost?

.. code:: ipython3

    print('operating cost = $', value(m.fs.operating_cost))


.. parsed-literal::

    operating cost = $ 419122.33876779396


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
              Heat Duty :      7352.5 : False : (None, None)
        Pressure Change : -2.0000e+05 :  True : (None, None)
    
    ------------------------------------------------------------------------------------
        Stream Table
                                                   Inlet    Vapor Outlet  Liquid Outlet
        flow_mol_phase_comp ('Liq', 'benzene')     0.20460   1.0000e-08      0.062620  
        flow_mol_phase_comp ('Liq', 'toluene')    0.062520   1.0000e-08      0.032257  
        flow_mol_phase_comp ('Liq', 'hydrogen') 2.6712e-07   1.0000e-08    9.4877e-08  
        flow_mol_phase_comp ('Liq', 'methane')  2.6712e-07   1.0000e-08    9.4877e-08  
        flow_mol_phase_comp ('Vap', 'benzene')  1.0000e-08      0.14198    1.0000e-08  
        flow_mol_phase_comp ('Vap', 'toluene')  1.0000e-08     0.030264    1.0000e-08  
        flow_mol_phase_comp ('Vap', 'hydrogen') 1.0000e-08   1.8224e-07    1.0000e-08  
        flow_mol_phase_comp ('Vap', 'methane')  1.0000e-08   1.8224e-07    1.0000e-08  
        temperature                                 325.00       375.00        375.00  
        pressure                                3.5000e+05   1.5000e+05    1.5000e+05  
    ====================================================================================
    
    benzene purity =  0.8242962943918918


Next, let's look at how much benzene we are loosing with the light gases
out of F101. IDAES has tools for creating stream tables based on the
``Arcs`` and/or ``Ports`` in a flowsheet. Let us create and print a
simple stream table showing the stream leaving the reactor and the vapor
stream from F101.

.. raw:: html

   <div class="alert alert-block alert-info">

Inline Exercise: How much benzene are we loosing in the F101 vapor
outlet stream?

.. raw:: html

   </div>

.. code:: ipython3

    from idaes.core.util.tables import create_stream_table_dataframe, stream_table_dataframe_to_string
    
    st = create_stream_table_dataframe({"Reactor": m.fs.s05, "Light Gases": m.fs.s06})
    print(stream_table_dataframe_to_string(st))


.. parsed-literal::

                                              Reactor   Light Gases
    flow_mol_phase_comp ('Liq', 'benzene')  1.2993e-07  1.0000e-08 
    flow_mol_phase_comp ('Liq', 'toluene')  8.4147e-07  1.0000e-08 
    flow_mol_phase_comp ('Liq', 'hydrogen') 1.0000e-08  1.0000e-08 
    flow_mol_phase_comp ('Liq', 'methane')  1.0000e-08  1.0000e-08 
    flow_mol_phase_comp ('Vap', 'benzene')     0.35374     0.14915 
    flow_mol_phase_comp ('Vap', 'toluene')    0.078129    0.015610 
    flow_mol_phase_comp ('Vap', 'hydrogen')    0.32821     0.32821 
    flow_mol_phase_comp ('Vap', 'methane')      1.2721      1.2721 
    temperature                                 771.85      325.00 
    pressure                                3.5000e+05  3.5000e+05 


.. raw:: html

   <div class="alert alert-block alert-info">

Inline Exercise: You can querry additional variables here if you like.

Use Shift+Enter to run the cell once you have typed in your code.

.. raw:: html

   </div>

Optimization
------------

We saw from the results above that the total operating cost for the base
case was $419,122 per year. We are producing 0.142 mol/s of benzene at a
purity of 82%. However, we are losing around 42% of benzene in F101
vapor outlet stream.

Let us try to minimize this cost such that: - we are producing at least
0.15 mol/s of benzene in F102 vapor outlet i.e. our product stream -
purity of benzne i.e. the mole fraction of benzene in F102 vapor outlet
is at least 80% - restricting the benzene loss in F101 vapor outlet to
less than 20%

For this problem, our decision variables are as follows: - H101 outlet
temperature - R101 cooling duty provided - F101 outlet temperature -
F102 outlet temperature - F102 deltaP in the flash tank

Let us declare our objective function for this problem.

.. code:: ipython3

    m.fs.objective = Objective(expr=m.fs.operating_cost)

Now, we need to unfix the decision variables as we had solved a square
problem (degrees of freedom = 0) until now.

.. code:: ipython3

    m.fs.H101.outlet.temperature.unfix()
    m.fs.R101.heat_duty.unfix()
    m.fs.F101.vap_outlet.temperature.unfix()
    m.fs.F102.vap_outlet.temperature.unfix()

.. raw:: html

   <div class="alert alert-block alert-info">

Inline Exercise: Let us now unfix the remaining variable which is F102
pressure drop (F102.deltaP)

Use Shift+Enter to run the cell once you have typed in your code.

.. raw:: html

   </div>

.. code:: ipython3

    m.fs.F102.deltaP.unfix()

Next, we need to set bounds on these decision variables to values shown
below:

-  H101 outlet temperature [500, 600] K
-  R101 outlet temperature [600, 800] K
-  F101 outlet temperature [298, 450] K
-  F102 outlet temperature [298, 450] K
-  F102 outlet pressure [105000, 110000] Pa

Let us first set the variable bound for the H101 outlet temperature as
shown below:

.. code:: ipython3

    m.fs.H101.outlet.temperature[0].setlb(500)
    m.fs.H101.outlet.temperature[0].setub(600)

.. raw:: html

   <div class="alert alert-block alert-info">

Inline Exercise: Now, set the variable bound for the R101 outlet
temperature.

Use Shift+Enter to run the cell once you have typed in your code.

.. raw:: html

   </div>

.. code:: ipython3

    m.fs.R101.outlet.temperature[0].setlb(600)
    m.fs.R101.outlet.temperature[0].setub(800)

Let us fix the bounds for the rest of the decision variables.

.. code:: ipython3

    m.fs.F101.vap_outlet.temperature[0].setlb(298.0)
    m.fs.F101.vap_outlet.temperature[0].setub(450.0)
    m.fs.F102.vap_outlet.temperature[0].setlb(298.0)
    m.fs.F102.vap_outlet.temperature[0].setub(450.0)
    m.fs.F102.vap_outlet.pressure[0].setlb(105000)
    m.fs.F102.vap_outlet.pressure[0].setub(110000)

Now, the only things left to define are our constraints on overhead loss
in F101, product flow rate and purity in F102. Let us first look at
defining a constraint for the overhead loss in F101 where we are
restricting the benzene leaving the vapor stream to less than 20 % of
the benzene available in the reactor outlet.

.. code:: ipython3

    m.fs.overhead_loss = Constraint(
            expr=m.fs.F101.vap_outlet.flow_mol_phase_comp[0, "Vap", "benzene"] <=
            0.20 * m.fs.R101.outlet.flow_mol_phase_comp[0, "Vap", "benzene"])

.. raw:: html

   <div class="alert alert-block alert-info">

Inline Exercise: Now, add the constraint such that we are producing at
least 0.15 mol/s of benzene in the product stream which is the vapor
outlet of F102. Let us name this constraint as m.fs.product\_flow.

Use Shift+Enter to run the cell once you have typed in your code.

.. raw:: html

   </div>

.. code:: ipython3

    m.fs.product_flow = Constraint(
            expr=m.fs.F102.vap_outlet.flow_mol_phase_comp[0, "Vap", "benzene"] >=
            0.15)

Let us add the final constraint on product purity or the mole fraction
of benzene in the product stream such that it is at least greater than
80%.

.. code:: ipython3

    m.fs.product_purity = Constraint(expr=m.fs.purity >= 0.80)

We have now defined the optimization problem and we are now ready to
solve this problem.

.. code:: ipython3

    results = solver.solve(m, tee=True)


.. parsed-literal::

    Ipopt 3.13.2: tol=1e-06
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
    
    This is Ipopt version 3.13.2, running with linear solver ma27.
    
    Number of nonzeros in equality constraint Jacobian...:     1048
    Number of nonzeros in inequality constraint Jacobian.:        5
    Number of nonzeros in Lagrangian Hessian.............:      901
    
    Total number of variables............................:      343
                         variables with only lower bounds:        0
                    variables with lower and upper bounds:      149
                         variables with only upper bounds:        0
    Total number of equality constraints.................:      338
    Total number of inequality constraints...............:        3
            inequality constraints with only lower bounds:        2
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        1
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  4.1912234e+05 2.99e+05 6.94e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  4.1628396e+05 2.99e+05 6.94e+00  -1.0 4.82e+09    -  1.80e-05 5.83e-06f  1
       2  4.1616769e+05 2.99e+05 1.60e+02  -1.0 1.45e+09    -  5.86e-04 1.47e-05f  1
       3  4.0783429e+05 2.94e+05 4.85e+02  -1.0 1.35e+09    -  2.61e-04 9.35e-04f  1
       4  2.9670827e+05 2.83e+06 6.94e+02  -1.0 4.75e+08    -  7.35e-05 1.50e-03f  1
       5  2.9557701e+05 2.82e+06 4.95e+04  -1.0 1.90e+08    -  1.88e-01 1.04e-03f  1
       6  2.9452502e+05 2.72e+06 4.63e+05  -1.0 4.40e+07    -  1.88e-01 3.46e-02f  1
       7  2.9632753e+05 2.13e+06 4.47e+05  -1.0 1.47e+07    -  7.61e-02 2.19e-01h  1
       8  2.9636923e+05 2.12e+06 4.45e+05  -1.0 5.86e+06    -  6.38e-01 3.38e-03h  1
       9  2.9647019e+05 2.10e+06 4.42e+05  -1.0 6.53e+06    -  7.25e-01 7.18e-03h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  2.9958737e+05 1.63e+06 4.17e+05  -1.0 6.56e+06    -  3.55e-02 2.24e-01h  1
      11  3.0436334e+05 9.49e+05 6.96e+05  -1.0 5.55e+06    -  9.46e-01 4.19e-01h  1
      12  3.0792618e+05 5.00e+05 4.56e+06  -1.0 4.03e+06    -  9.90e-01 4.73e-01h  1
      13  3.0931998e+05 3.46e+05 1.59e+08  -1.0 2.67e+06    -  1.00e+00 3.08e-01h  2
      14  3.1261432e+05 5.80e+05 1.20e+11  -1.0 2.00e+06    -  1.00e+00 9.78e-01H  1
      15  3.1271509e+05 2.42e+05 8.71e+08  -1.0 1.43e+05    -  1.00e+00 5.84e-01h  1
      16  3.1278603e+05 2.73e+03 3.26e+11  -1.0 5.97e+04    -  1.00e+00 9.90e-01h  1
      17  3.1278674e+05 1.79e-01 3.96e+09  -1.0 6.18e+02    -  1.00e+00 1.00e+00h  1
      18  3.1278674e+05 1.91e-06 3.15e+05  -1.0 5.18e-03    -  1.00e+00 1.00e+00h  1
      19  3.1278634e+05 3.47e-05 1.62e+06  -3.8 2.02e+02    -  1.00e+00 1.00e+00f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      20  3.1278634e+05 1.49e-08 1.31e-03  -3.8 1.71e-01    -  1.00e+00 1.00e+00h  1
      21  3.1278634e+05 7.45e-09 3.81e+00  -7.0 3.04e-01    -  1.00e+00 1.00e+00f  1
      22  3.1278634e+05 7.45e-09 5.70e-05  -7.0 3.91e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 22
    
                                       (scaled)                 (unscaled)
    Objective...............:   3.1278633834102668e+05    3.1278633834102668e+05
    Dual infeasibility......:   5.7004195292406851e-05    5.7004195292406851e-05
    Constraint violation....:   2.9103830456733704e-11    7.4505805969238281e-09
    Complementarity.........:   9.0926527280252916e-08    9.0926527280252916e-08
    Overall NLP error.......:   2.0798568665482607e-09    5.7004195292406851e-05
    
    
    Number of objective function evaluations             = 26
    Number of objective gradient evaluations             = 23
    Number of equality constraint evaluations            = 26
    Number of inequality constraint evaluations          = 26
    Number of equality constraint Jacobian evaluations   = 23
    Number of inequality constraint Jacobian evaluations = 23
    Number of Lagrangian Hessian evaluations             = 22
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.018
    Total CPU secs in NLP function evaluations           =      0.003
    
    EXIT: Optimal Solution Found.


.. code:: ipython3

    # For testing purposes
    from pyomo.environ import TerminationCondition
    assert results.solver.termination_condition == TerminationCondition.optimal

Optimization Results
--------------------

Display the results and product specifications

.. code:: ipython3

    print('operating cost = $', value(m.fs.operating_cost))
    
    print()
    print('Product flow rate and purity in F102')
    
    m.fs.F102.report()
    
    print()
    print('benzene purity = ', value(m.fs.purity))
    
    print()
    print('Overhead loss in F101')
    m.fs.F101.report()


.. parsed-literal::

    operating cost = $ 312786.3383410267
    
    Product flow rate and purity in F102
    
    ====================================================================================
    Unit : fs.F102                                                             Time: 0.0
    ------------------------------------------------------------------------------------
        Unit Performance
    
        Variables: 
    
        Key             : Value       : Fixed : Bounds
              Heat Duty :      8377.0 : False : (None, None)
        Pressure Change : -2.4500e+05 : False : (None, None)
    
    ------------------------------------------------------------------------------------
        Stream Table
                                                   Inlet    Vapor Outlet  Liquid Outlet
        flow_mol_phase_comp ('Liq', 'benzene')     0.21743   1.0000e-08      0.067425  
        flow_mol_phase_comp ('Liq', 'toluene')    0.070695   1.0000e-08      0.037507  
        flow_mol_phase_comp ('Liq', 'hydrogen') 2.8812e-07   1.0000e-08    1.0493e-07  
        flow_mol_phase_comp ('Liq', 'methane')  2.8812e-07   1.0000e-08    1.0493e-07  
        flow_mol_phase_comp ('Vap', 'benzene')  1.0000e-08      0.15000    1.0000e-08  
        flow_mol_phase_comp ('Vap', 'toluene')  1.0000e-08     0.033189    1.0000e-08  
        flow_mol_phase_comp ('Vap', 'hydrogen') 1.0000e-08   1.9319e-07    1.0000e-08  
        flow_mol_phase_comp ('Vap', 'methane')  1.0000e-08   1.9319e-07    1.0000e-08  
        temperature                                 301.88       362.93        362.93  
        pressure                                3.5000e+05   1.0500e+05    1.0500e+05  
    ====================================================================================
    
    benzene purity =  0.8188276578112267
    
    Overhead loss in F101
    
    ====================================================================================
    Unit : fs.F101                                                             Time: 0.0
    ------------------------------------------------------------------------------------
        Unit Performance
    
        Variables: 
    
        Key             : Value   : Fixed : Bounds
              Heat Duty : -56354. : False : (None, None)
        Pressure Change :  0.0000 :  True : (None, None)
    
    ------------------------------------------------------------------------------------
        Stream Table
                                                   Inlet    Vapor Outlet  Liquid Outlet
        flow_mol_phase_comp ('Liq', 'benzene')  4.3534e-08   1.0000e-08       0.21743  
        flow_mol_phase_comp ('Liq', 'toluene')  7.5866e-07   1.0000e-08      0.070695  
        flow_mol_phase_comp ('Liq', 'hydrogen') 1.0000e-08   1.0000e-08    2.8812e-07  
        flow_mol_phase_comp ('Liq', 'methane')  1.0000e-08   1.0000e-08    2.8812e-07  
        flow_mol_phase_comp ('Vap', 'benzene')     0.27178     0.054356    1.0000e-08  
        flow_mol_phase_comp ('Vap', 'toluene')    0.076085    0.0053908    1.0000e-08  
        flow_mol_phase_comp ('Vap', 'hydrogen')    0.35887      0.35887    1.0000e-08  
        flow_mol_phase_comp ('Vap', 'methane')      1.2414       1.2414    1.0000e-08  
        temperature                                 696.12       301.88        301.88  
        pressure                                3.5000e+05   3.5000e+05    3.5000e+05  
    ====================================================================================


Display optimal values for the decision variables

.. code:: ipython3

    print('Optimal Values')
    print()
    
    print('H101 outlet temperature = ', value(m.fs.H101.outlet.temperature[0]), 'K')
    
    print()
    print('R101 outlet temperature = ', value(m.fs.R101.outlet.temperature[0]), 'K')
    
    print()
    print('F101 outlet temperature = ', value(m.fs.F101.vap_outlet.temperature[0]), 'K')
    
    print()
    print('F102 outlet temperature = ', value(m.fs.F102.vap_outlet.temperature[0]), 'K')
    print('F102 outlet pressure = ', value(m.fs.F102.vap_outlet.pressure[0]), 'Pa')


.. parsed-literal::

    Optimal Values
    
    H101 outlet temperature =  500.0 K
    
    R101 outlet temperature =  696.1161004637527 K
    
    F101 outlet temperature =  301.87847605692815 K
    
    F102 outlet temperature =  362.93476830548985 K
    F102 outlet pressure =  105000.0 Pa



