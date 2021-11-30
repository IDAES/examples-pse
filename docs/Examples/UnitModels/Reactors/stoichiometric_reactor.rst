Flowsheet Stoichiometric Reactor Simulation and Optimization of Ethylene Glycol Production
==========================================================================================

Learning outcomes
-----------------

-  Call and implement the IDAES StochiometricReactor unit model
-  Construct a steady-state flowsheet using the IDAES unit model library
-  Connecting unit models in a flowsheet using Arcs
-  Fomulate and solve an optimization problem

   -  Defining an objective function
   -  Setting variable bounds
   -  Adding additional constraints

Problem Statement
-----------------

Following the previous example implementing a CSTR unit model, we can
alter the flowsheet to use a stochiometric (or yield) reactor. As
before, this example is adapted from Fogler, H.S., Elements of Chemical
Reaction Engineering 5th ed., 2016, Prentice Hall, p. 157-160 with the
following chemical reaction, property packages and flowsheet. Unlike the
previous two reactors which apply performance equations to calculate
reaction extent, this simplified reactor model neglects all geometric
properties and allows the user to specify a yield per reaction. The
state variables chosen for the property package are **molar flows of
each component by phase in each stream, temperature of each stream and
pressure of each stream**. The components considered are: **ethylene
oxide, water, sulfuric acid and ethylene glycol** and the process occurs
in liquid phase only. Therefore, every stream has 4 flow variables, 1
temperature and 1 pressure variable.

Chemical reaction:

**C2H4O + H2O + H2SO4 → C2H6O2 + H2SO4**

Property Packages:

-  egprod_ideal.py
-  egprod_reaction.py

Flowsheet:

.. image:: egprod_flowsheet.png

Importing required Pyomo and IDAES components
---------------------------------------------

To construct a flowsheet, we will need several components from the pyomo
and idaes package. Let us first import the following components from
Pyomo: - Constraint (to write constraints) - Var (to declare variables)
- ConcreteModel (to create the concrete model object) - Expression (to
evaluate values as a function of variables defined in the model) -
Objective (to define an objective function for optimization) -
TransformationFactory (to apply certain transformations) - Arc (to
connect two unit models)

For further details on these components, please refer to the pyomo
documentation: https://pyomo.readthedocs.io/en/latest/

From idaes, we will be needing the FlowsheetBlock and the following unit
models: - Mixer - Heater - StoichiometricReactor

We will also be needing some utility tools to put together the flowsheet
and calculate the degrees of freedom, tools for model expressions and
calling variable values, and built-in functions to define property
packages, add unit containers to objects and define our initialization
scheme.

.. code:: ipython3

    from pyomo.environ import (Constraint,
                               exp,
                               Var,
                               ConcreteModel,
                               Expression,
                               Objective,
                               TransformationFactory,
                               value,
                               units as pyunits)
    from pyomo.network import Arc
    
    from idaes.core import FlowsheetBlock
    from idaes.generic_models.properties.core.generic.generic_property import (
            GenericParameterBlock)
    from idaes.generic_models.properties.core.generic.generic_reaction import (
            GenericReactionParameterBlock)
    from idaes.generic_models.unit_models import (Mixer,
                                                  Heater,
                                                  StoichiometricReactor)
    
    from idaes.core.util.constants import Constants
    from idaes.core.util import get_solver
    from idaes.core.util.model_statistics import degrees_of_freedom
    from idaes.core.util.initialization import propagate_state
    
    import idaes.logger as idaeslog

Importing required thermo and reaction package
----------------------------------------------

The final set of imports are to import the thermo and reaction package.
We have created a custom thermo package that support ideal vapor and
liquid behavior for this system, and in this case we will restrict it to
ideal liquid behavior only.

The reaction package here assumes Arrhenius kinetic behavior for the
StoichiometricReactor, for which :math:`k_0` and :math:`E_a` are known
*a priori* (if unknown, they may be obtained using one of the parameter
estimation tools within IDAES).

$ r = -kVC_{EO} $, $ k = k_0 e^{(-E_a/RT)}$, with the variables as
follows:

| :math:`r` - reaction rate extent in moles of ethylene oxide consumed
  per second; note that the traditional reaction rate would be given by
  :math:`rate = r/V` in moles per :math:`m^3` per second
| :math:`k` - reaction rate constant per second
| :math:`V` - volume of StoichiometricReactor in :math:`m^3`, note that
  this is *liquid volume* and not the *total volume* of the reactor
  itself
| :math:`C_{EO}` - bulk concentration of ethylene oxide in moles per
  :math:`m^3` (the limiting reagent, since we assume excess catalyst and
  water)
| :math:`k_0` - pre-exponential Arrhenius factor per second
| :math:`E_a` - reaction activation energy in kJ per mole of ethylene
  oxide consumed
| :math:`R` - gas constant in J/mol-K
| :math:`T` - reactor temperature in K

These calculations are contained within the property, reaction and unit
model packages, and do not need to be entered into the flowsheet. More
information on property estimation may be found below:

| ParamEst parameter estimation:
  https://idaes-pse.readthedocs.io/en/stable/user_guide/workflow/data_rec_parmest.html?highlight=paramest
| HELMET thermodynamic estimation:
  https://idaes-pse.readthedocs.io/en/stable/user_guide/modeling_extensions/surrogate/helmet/index.html
| RIPE reaction estimation:
  https://idaes-pse.readthedocs.io/en/stable/user_guide/modeling_extensions/surrogate/ripe/index.html

Let us import the following modules from the same directory as this
Jupyter notebook: - egprod_ideal as thermo_props - egprod_reaction as
reaction_props

.. code:: ipython3

    import egprod_ideal as thermo_props
    import egprod_reaction as reaction_props

Constructing the Flowsheet
--------------------------

We have now imported all the components, unit models, and property
modules we need to construct a flowsheet. Let us create a ConcreteModel
and add the flowsheet block.

.. code:: ipython3

    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

We now need to add the property packages to the flowsheet. Unlike Module
1, where we only had a thermo property package, for this flowsheet we
will also need to add a reaction property package. We will use the
Generic Property and Generic Reaction Frameworks; more information may
be found on these methods at
https://idaes-pse.readthedocs.io/en/1.8.0/user_guide/components/property_package/index.html.

.. code:: ipython3

    m.fs.thermo_params = GenericParameterBlock(default=thermo_props.config_dict)
    m.fs.reaction_params = GenericReactionParameterBlock(default={"property_package": m.fs.thermo_params,
                                                                  **reaction_props.config_dict})

Adding Unit Models
------------------

Let us start adding the unit models we have imported to the flowsheet.
Here, we are adding a Mixer (assigned a name M101), a Heater (assigned a
name H101) and a StoichiometricReactor (assigned a name R101). Note that
all unit models need to be given a property package argument. In
addition to that, there are several arguments depending on the unit
model, please refer to the documentation for more details
(https://idaes-pse.readthedocs.io/en/latest/model_libraries/core_lib/unit_models/index.html).
For example, the Mixer unit model here is given a ``list`` consisting of
names to the two inlets.

.. code:: ipython3

    m.fs.M101 = Mixer(default={"property_package": m.fs.thermo_params,
                               "inlet_list": ["reagent_feed", "catalyst_feed"]})
    m.fs.H101 = Heater(default={"property_package": m.fs.thermo_params,
                                "has_pressure_change": False,
                                "has_phase_equilibrium": False})

.. code:: ipython3

    m.fs.R101 = StoichiometricReactor(
                default={"property_package": m.fs.thermo_params,
                         "reaction_package": m.fs.reaction_params,
                         "has_heat_of_reaction": True,
                         "has_heat_transfer": True,
                         "has_pressure_change": False})

Connecting Unit Models using Arcs
---------------------------------

We have now added all the unit models we need to the flowsheet. However,
we have not yet specifed how the units are to be connected. To do this,
we will be using the ``Arc`` which is a pyomo component that takes in
two arguments: ``source`` and ``destination``. Let us connect the outlet
of the mixer(M101) to the inlet of the heater(H101), and the outlet of
the heater(H101) to the inlet of the reactor(R101).

.. code:: ipython3

    m.fs.s03 = Arc(source=m.fs.M101.outlet, destination=m.fs.H101.inlet)
    m.fs.s04 = Arc(source=m.fs.H101.outlet, destination=m.fs.R101.inlet)

We have now connected the unit model block using the arcs. However, each
of these arcs link to ports on the two unit models that are connected.
In this case, the ports consist of the state variables that need to be
linked between the unit models. Pyomo provides a convenient method to
write these equality constraints for us between two ports and this is
done as follows:

.. code:: ipython3

    TransformationFactory("network.expand_arcs").apply_to(m)

Adding expressions to compute operating costs
---------------------------------------------

In this section, we will add a few Expressions that allows us to
evaluate the performance. Expressions provide a convenient way of
calculating certain values that are a function of the variables defined
in the model. For more details on Expressions, please refer to:
https://pyomo.readthedocs.io/en/latest/pyomo_modeling_components/Expressions.html

For this flowsheet, we are interested in computing ethylene glycol
production in millions of pounds per year, as well as the total costs
due to cooling and heating utilities:

Let us first add an Expression to convert the product flow from mol/s to
MM lb/year of ethylene glycol. We see that our molecular weight exists
in the thermo property package, so we may use that value for our
calculations.

.. code:: ipython3

    m.fs.eg_prod = Expression(expr=pyunits.convert(m.fs.R101.outlet.flow_mol_phase_comp[0, "Liq", "ethylene_glycol"]
                                                   *m.fs.thermo_params.ethylene_glycol.mw, # MW defined in properties as kg/mol
                                                   to_units=pyunits.Mlb/pyunits.yr)) # converting kg/s to MM lb/year

Now, let us add expressions to compute the reactor cooling cost
(\\\ :math:`/s) assuming a cost of 0.212E-4 \\`/kW, and the heating
utility cost (\\\ :math:`/s) assuming 2.2E-4 \\`/kW. Note that the heat
duty is in units of watt (J/s). The total operating cost will be the sum
of the two, expressed in \\$/year assuming 8000 operating hours per year
(~10% downtime, which is fairly common for small scale chemical plants):

.. code:: ipython3

    m.fs.cooling_cost = Expression(expr=0.212e-7 * (-m.fs.R101.heat_duty[0]))  # the reaction is exothermic, so R101 duty is negative
    m.fs.heating_cost = Expression(expr=2.2e-7 * m.fs.H101.heat_duty[0])  # the stream must be heated to T_rxn, so H101 duty is positive
    m.fs.operating_cost = Expression(expr=(3600 * 8000 *(m.fs.heating_cost + m.fs.cooling_cost)))

Fixing feed conditions
----------------------

Let us first check how many degrees of freedom exist for this flowsheet
using the ``degrees_of_freedom`` tool we imported earlier. We expect
each stream to have 6 degrees of freedom, the mixer to have 0 (after
both streams are accounted for), the heater to have 1 (just the duty,
since the inlet is also the outlet of M101), and the reactor to have 1
(duty or overall conversion, since the inlet is also the outlet of
H101). In this case, the reactor has an extra degree of freedom since we
have not yet defined the yield of the sole rate-kinetics reaction.
Therefore, we have 15 degrees of freedom to specify: temperature,
pressure and flow of all four components on both streams; outlet heater
temperature; reactor conversion and duty.

.. code:: ipython3

    print(degrees_of_freedom(m))


.. parsed-literal::

    15
    

We will now be fixing the feed stream to the conditions shown in the
flowsheet above. As mentioned in other tutorials, the IDAES framework
expects a time index value for every referenced internal stream or unit
variable, even in steady-state systems with a single time point $ t = 0
$. The non-present components in each stream are assigned a very small
non-zero value to help with convergence and initializing. Based on
stoichiometric ratios for the reaction, 80% conversion and 200 MM
lb/year (46.4 mol/s) of ethylene glycol, we will initialize our
simulation with the following calculated values:

.. code:: ipython3

    m.fs.M101.reagent_feed.flow_mol_phase_comp[0, "Liq", "ethylene_oxide"].fix(58.0*pyunits.mol/pyunits.s)
    m.fs.M101.reagent_feed.flow_mol_phase_comp[0, "Liq", "water"].fix(39.6*pyunits.mol/pyunits.s)  # calculated from 16.1 mol EO / cudm in stream
    m.fs.M101.reagent_feed.flow_mol_phase_comp[0, "Liq", "sulfuric_acid"].fix(1e-5*pyunits.mol/pyunits.s)
    m.fs.M101.reagent_feed.flow_mol_phase_comp[0, "Liq", "ethylene_glycol"].fix(1e-5*pyunits.mol/pyunits.s)
    m.fs.M101.reagent_feed.temperature.fix(298.15*pyunits.K)
    m.fs.M101.reagent_feed.pressure.fix(1e5*pyunits.Pa)
    
    m.fs.M101.catalyst_feed.flow_mol_phase_comp[0, "Liq", "ethylene_oxide"].fix(1e-5*pyunits.mol/pyunits.s)
    m.fs.M101.catalyst_feed.flow_mol_phase_comp[0, "Liq", "water"].fix(200*pyunits.mol/pyunits.s)
    m.fs.M101.catalyst_feed.flow_mol_phase_comp[0, "Liq", "sulfuric_acid"].fix(0.334*pyunits.mol/pyunits.s)  # calculated from 0.9 wt% SA in stream
    m.fs.M101.catalyst_feed.flow_mol_phase_comp[0, "Liq", "ethylene_glycol"].fix(1e-5*pyunits.mol/pyunits.s)
    m.fs.M101.catalyst_feed.temperature.fix(298.15*pyunits.K)
    m.fs.M101.catalyst_feed.pressure.fix(1e5*pyunits.Pa)

Fixing unit model specifications
--------------------------------

Now that we have fixed our inlet feed conditions, we will now be fixing
the operating conditions for the unit models in the flowsheet. Let us
fix the outlet temperature of H101 to 328.15 K.

.. code:: ipython3

    m.fs.H101.outlet.temperature.fix(328.15*pyunits.K)

Unlike the previous two examples, we will need to specify both initial
conversion and heat duty values (these are the only two free variables
to choose from). Since heat duty and the outlet reactor temperature are
interdependent, we can choose to specify this quantity instead. While
the reaction kinetic parameters exist in the property package, we also
do not need to add a rate constant expression since generation is
explicitly defined through the conversion/yield. Note that our initial
problem will solve with zero *temperature change* but will be infeasible
with zero *heat duty*; this is due to the heat of reaction enforced by
allowing heat transfer and mandating a non-zero conversion.

.. code:: ipython3

    m.fs.R101.conversion = Var(initialize=0.80, bounds=(0, 1), units=pyunits.dimensionless)  # fraction
    
    m.fs.R101.conv_constraint = Constraint(
        expr=m.fs.R101.conversion*m.fs.R101.inlet.
        flow_mol_phase_comp[0, "Liq", "ethylene_oxide"] ==
        (m.fs.R101.inlet.flow_mol_phase_comp[0, "Liq", "ethylene_oxide"] -
         m.fs.R101.outlet.flow_mol_phase_comp[0, "Liq", "ethylene_oxide"]))
    
    m.fs.R101.conversion.fix(0.80)
    
    m.fs.R101.outlet.temperature.fix(328.15*pyunits.K)  # equal inlet reactor temperature

.. code:: ipython3

    print(degrees_of_freedom(m))


.. parsed-literal::

    0
    

Finally, we need to initialize the each unit operation in sequence to
solve the flowsheet. In best practice, unit operations are initialized
or solved, and outlet properties are propagated to connected inlet
streams via arc definitions as follows:

.. code:: ipython3

    # Initialize and solve each unit operation
    m.fs.M101.initialize()
    propagate_state(arc=m.fs.s03)
    
    m.fs.H101.initialize()
    propagate_state(arc=m.fs.s04)
    
    m.fs.R101.initialize()
    
    # set solver
    solver = get_solver()


.. parsed-literal::

    2021-11-30 12:29:36 [INFO] idaes.init.fs.M101.reagent_feed_state: Starting initialization
    2021-11-30 12:29:36 [INFO] idaes.init.fs.M101.reagent_feed_state: Property initialization: optimal - Optimal Solution Found.
    2021-11-30 12:29:36 [INFO] idaes.init.fs.M101.catalyst_feed_state: Starting initialization
    2021-11-30 12:29:36 [INFO] idaes.init.fs.M101.catalyst_feed_state: Property initialization: optimal - Optimal Solution Found.
    2021-11-30 12:29:36 [INFO] idaes.init.fs.M101.mixed_state: Starting initialization
    2021-11-30 12:29:36 [INFO] idaes.init.fs.M101.mixed_state: Property initialization: optimal - Optimal Solution Found.
    2021-11-30 12:29:36 [INFO] idaes.init.fs.M101.mixed_state: Property package initialization: optimal - Optimal Solution Found.
    2021-11-30 12:29:36 [INFO] idaes.init.fs.M101: Initialization Complete: optimal - Optimal Solution Found
    2021-11-30 12:29:36 [INFO] idaes.init.fs.H101.control_volume.properties_in: Starting initialization
    2021-11-30 12:29:36 [INFO] idaes.init.fs.H101.control_volume.properties_in: Property initialization: optimal - Optimal Solution Found.
    2021-11-30 12:29:36 [INFO] idaes.init.fs.H101.control_volume.properties_out: Starting initialization
    2021-11-30 12:29:37 [INFO] idaes.init.fs.H101.control_volume.properties_out: Property initialization: optimal - Optimal Solution Found.
    2021-11-30 12:29:37 [INFO] idaes.init.fs.H101.control_volume: Initialization Complete
    2021-11-30 12:29:37 [INFO] idaes.init.fs.H101: Initialization Complete: optimal - Optimal Solution Found
    2021-11-30 12:29:37 [INFO] idaes.init.fs.R101.control_volume.properties_in: Starting initialization
    2021-11-30 12:29:37 [INFO] idaes.init.fs.R101.control_volume.properties_in: Property initialization: optimal - Optimal Solution Found.
    2021-11-30 12:29:37 [INFO] idaes.init.fs.R101.control_volume.properties_out: Starting initialization
    2021-11-30 12:29:37 [INFO] idaes.init.fs.R101.control_volume.properties_out: Property initialization: optimal - Optimal Solution Found.
    2021-11-30 12:29:37 [INFO] idaes.init.fs.R101.control_volume.reactions: Initialization Complete.
    2021-11-30 12:29:37 [INFO] idaes.init.fs.R101.control_volume: Initialization Complete
    2021-11-30 12:29:37 [INFO] idaes.init.fs.R101: Initialization Complete: optimal - Optimal Solution Found
    

.. code:: ipython3

    # Solve the model
    results = solver.solve(m, tee=True)


.. parsed-literal::

    Ipopt 3.13.2: nlp_scaling_method=gradient-based
    tol=1e-06
    
    
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
    
    Number of nonzeros in equality constraint Jacobian...:      232
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:      221
    
    Total number of variables............................:       65
                         variables with only lower bounds:        0
                    variables with lower and upper bounds:       56
                         variables with only upper bounds:        0
    Total number of equality constraints.................:       65
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  0.0000000e+00 1.76e+06 0.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  0.0000000e+00 1.77e+04 6.74e-04  -1.0 1.00e-02    -  9.90e-01 9.90e-01h  1
       2  0.0000000e+00 1.67e+02 5.83e-01  -1.0 1.00e-04    -  9.90e-01 9.91e-01h  1
       3  0.0000000e+00 1.10e-06 8.89e+02  -1.0 9.43e-07    -  9.91e-01 1.00e+00h  1
    Cannot recompute multipliers for feasibility problem.  Error in eq_mult_calculator
    
    Number of Iterations....: 3
    
                                       (scaled)                 (unscaled)
    Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    Dual infeasibility......:   1.0094083379891624e+05    1.0094083379891624e+05
    Constraint violation....:   7.2759576141834259e-12    1.1026859283447266e-06
    Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    Overall NLP error.......:   7.2759576141834259e-12    1.0094083379891624e+05
    
    
    Number of objective function evaluations             = 4
    Number of objective gradient evaluations             = 4
    Number of equality constraint evaluations            = 4
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 4
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 3
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.004
    Total CPU secs in NLP function evaluations           =      0.000
    
    EXIT: Optimal Solution Found.
    

Analyze the results of the square problem
-----------------------------------------

What is the total operating cost?

.. code:: ipython3

    print('operating cost = $', value(m.fs.operating_cost), ' per year')


.. parsed-literal::

    operating cost = $ 3458139.885712622  per year
    

For this operating cost, what conversion did we achieve of ethylene
oxide to ethylene glycol?

.. code:: ipython3

    m.fs.R101.report()
    
    print()
    print('Conversion achieved = ', value(m.fs.R101.conversion)*100, '%')


.. parsed-literal::

    
    ====================================================================================
    Unit : fs.R101                                                             Time: 0.0
    ------------------------------------------------------------------------------------
        Unit Performance
    
        Variables: 
    
        Key                  : Value       : Fixed : Bounds
                   Heat Duty : -5.6566e+06 : False : (None, None)
        Reaction Extent [R1] :      46.400 : False : (None, None)
    
    ------------------------------------------------------------------------------------
        Stream Table
                                                     Inlet     Outlet  
        Molar Flowrate ('Liq', 'ethylene_oxide')      58.000     11.600
        Molar Flowrate ('Liq', 'water')               239.60     193.20
        Molar Flowrate ('Liq', 'sulfuric_acid')      0.33401    0.33401
        Molar Flowrate ('Liq', 'ethylene_glycol') 2.0000e-05     46.400
        Temperature                                   328.15     328.15
        Pressure                                  1.0000e+05 1.0000e+05
    ====================================================================================
    
    Conversion achieved =  80.0 %
    

Optimizing Ethylene Glycol Production
-------------------------------------

Now that the flowsheet has been squared and solved, we can run a small
optimization problem to minimize our production costs. Suppose we
require at least 200 million pounds/year of ethylene glycol produced and
at least 90% conversion of ethylene oxide, allowing for variable reactor
heat duty/outlet temperature and reactor temperature (heater outlet).

Let us declare our objective function for this problem.

.. code:: ipython3

    m.fs.objective = Objective(expr=m.fs.operating_cost)

Now, we need to add the design constraints and unfix the decision
variables as we had solved a square problem (degrees of freedom = 0)
until now, as well as set bounds for the design variables (reactor
outlet temperature is set by state variable bounds in property package):

.. code:: ipython3

    m.fs.eg_prod_con = Constraint(expr=m.fs.eg_prod >= 200*pyunits.Mlb/pyunits.yr)  # MM lb/year
    
    m.fs.conversion_con = Constraint(expr=m.fs.R101.conversion >= 0.90*pyunits.dimensionless)  # fraction
    
    m.fs.R101.conversion.unfix()
    
    m.fs.H101.outlet.temperature.unfix()
    m.fs.H101.outlet.temperature[0].setlb(328.15*pyunits.K)
    m.fs.H101.outlet.temperature[0].setub(470.45*pyunits.K)  # highest component boiling point (ethylene glycol)
    
    m.fs.R101.outlet.temperature.unfix()

We have now defined the optimization problem and we are now ready to
solve this problem.

.. code:: ipython3

    results = solver.solve(m, tee=True)


.. parsed-literal::

    Ipopt 3.13.2: nlp_scaling_method=gradient-based
    tol=1e-06
    
    
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
    
    Number of nonzeros in equality constraint Jacobian...:      236
    Number of nonzeros in inequality constraint Jacobian.:        2
    Number of nonzeros in Lagrangian Hessian.............:      242
    
    Total number of variables............................:       68
                         variables with only lower bounds:        0
                    variables with lower and upper bounds:       59
                         variables with only upper bounds:        0
    Total number of equality constraints.................:       65
    Total number of inequality constraints...............:        2
            inequality constraints with only lower bounds:        2
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  3.4581399e+06 1.76e+06 6.34e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  3.4240101e+06 1.76e+06 2.48e+01  -1.0 9.58e+08    -  1.46e-03 5.83e-05f  1
       2  3.4236530e+06 1.76e+06 1.56e+03  -1.0 3.49e+05    -  9.07e-02 1.48e-03f  1
       3  3.4237931e+06 1.73e+06 1.90e+03  -1.0 1.56e+04    -  4.15e-02 1.78e-02h  1
       4  3.4237959e+06 1.73e+06 2.80e+06  -1.0 2.65e+04    -  9.44e-01 1.81e-04h  1
       5  3.4238022e+06 1.73e+06 1.32e+06  -1.0 6.76e+05    -  7.17e-06 1.54e-05h  1
       6  3.4269127e+06 1.72e+06 3.06e+06  -1.0 7.60e+05    -  1.98e-05 6.70e-03h  1
       7  3.8613024e+06 1.12e+05 3.91e+05  -1.0 7.55e+05    -  6.70e-03 9.42e-01h  1
       8  3.8850641e+06 1.26e+04 2.25e+04  -1.0 4.38e+04    -  9.90e-01 8.88e-01h  1
       9  3.8880240e+06 1.15e+02 2.06e+02  -1.0 4.89e+03    -  9.90e-01 9.91e-01h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  3.8880508e+06 7.17e-05 1.49e-01  -3.8 4.39e+01    -  1.00e+00 1.00e+00h  1
      11  3.8880508e+06 1.46e-08 4.13e-07  -7.0 5.17e-04    -  1.00e+00 1.00e+00f  1
    
    Number of Iterations....: 11
    
                                       (scaled)                 (unscaled)
    Objective...............:   3.8880507982492442e+06    3.8880507982492442e+06
    Dual infeasibility......:   4.1300224027586061e-07    4.1300224027586061e-07
    Constraint violation....:   7.2759576141834259e-12    1.4627630662289448e-08
    Complementarity.........:   9.1066277388888538e-08    9.1066277388888538e-08
    Overall NLP error.......:   2.5311985076228633e-10    4.1300224027586061e-07
    
    
    Number of objective function evaluations             = 12
    Number of objective gradient evaluations             = 12
    Number of equality constraint evaluations            = 12
    Number of inequality constraint evaluations          = 12
    Number of equality constraint Jacobian evaluations   = 12
    Number of inequality constraint Jacobian evaluations = 12
    Number of Lagrangian Hessian evaluations             = 11
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.016
    Total CPU secs in NLP function evaluations           =      0.000
    
    EXIT: Optimal Solution Found.
    

.. code:: ipython3

    print('operating cost = $', value(m.fs.operating_cost), 'per year')
    
    print()
    print('Heater results')
    
    m.fs.H101.report()
    
    print()
    print('Stoichiometric reactor results')
    
    m.fs.R101.report()


.. parsed-literal::

    operating cost = $ 3888050.798249244 per year
    
    Heater results
    
    ====================================================================================
    Unit : fs.H101                                                             Time: 0.0
    ------------------------------------------------------------------------------------
        Unit Performance
    
        Variables: 
    
        Key       : Value  : Fixed : Bounds
        Heat Duty : 699.26 : False : (None, None)
    
    ------------------------------------------------------------------------------------
        Stream Table
                                                     Inlet     Outlet  
        Molar Flowrate ('Liq', 'ethylene_oxide')      58.000     58.000
        Molar Flowrate ('Liq', 'water')               239.60     239.60
        Molar Flowrate ('Liq', 'sulfuric_acid')      0.33401    0.33401
        Molar Flowrate ('Liq', 'ethylene_glycol') 2.0000e-05 2.0000e-05
        Temperature                                   298.15     328.15
        Pressure                                  1.0000e+05 1.0000e+05
    ====================================================================================
    
    Stoichiometric reactor results
    
    ====================================================================================
    Unit : fs.R101                                                             Time: 0.0
    ------------------------------------------------------------------------------------
        Unit Performance
    
        Variables: 
    
        Key                  : Value       : Fixed : Bounds
                   Heat Duty : -6.3608e+06 : False : (None, None)
        Reaction Extent [R1] :      52.200 : False : (None, None)
    
    ------------------------------------------------------------------------------------
        Stream Table
                                                     Inlet     Outlet  
        Molar Flowrate ('Liq', 'ethylene_oxide')      58.000     5.8000
        Molar Flowrate ('Liq', 'water')               239.60     187.40
        Molar Flowrate ('Liq', 'sulfuric_acid')      0.33401    0.33401
        Molar Flowrate ('Liq', 'ethylene_glycol') 2.0000e-05     52.200
        Temperature                                   328.15     450.00
        Pressure                                  1.0000e+05 1.0000e+05
    ====================================================================================
    

Display optimal values for the decision variables and design variables:

.. code:: ipython3

    print('Optimal Values')
    print()
    
    print('H101 outlet temperature = ', value(m.fs.H101.outlet.temperature[0]), 'K')
    
    print()
    print('R101 outlet temperature = ', value(m.fs.R101.outlet.temperature[0]), 'K')
    
    print()
    print('Ethylene glycol produced = ', value(m.fs.eg_prod), 'MM lb/year')
    
    print()
    print('Conversion achieved = ', value(m.fs.R101.conversion)*100, ' %')


.. parsed-literal::

    Optimal Values
    
    H101 outlet temperature =  328.15 K
    
    R101 outlet temperature =  450.0 K
    
    Ethylene glycol produced =  225.41546823488122 MM lb/year
    
    Conversion achieved =  89.9999990000021  %
    

