Tutorial: Mixer Unit Model with Ideal Property Package
======================================================

.. image:: mixer.svg

Learning Outcomes
-----------------

-  Demonstrate use of the Mixer unit model in IDAES
-  Demonstrate different options available

Problem Statement
-----------------

In this example, we will be mixing liquid benzene and liquid toluene
streams to form a mixture. The inlet conditions are as follows:

**Stream 1:**

Benzene Flow Rate = 100 mol/s

Pressure = 101325 Pa

Temperature = 353 K

**Stream 2**

Toluene Flow Rate = 100 mol/s

Pressure = 202650 Pa

Temperature = 356 K

We will look at two cases in this tutorial:

-  Case 1: Specify the number of inlets to the mixer, and set the
   ``momentum_mixing`` type set to “minimize”

-  Case 2: Specify the inlet names, and set ``momentum_mixing`` type set
   to “equality” (in this case, pressure will be specified for only one
   inlet stream)

**Note: When the momentum mixing type is set to ‘minimize’, the mixed
stream pressure takes the minimum value among all inlet stream
pressures. When the momentum mixing type is set to ‘equality’, the mixed
stream, along with all inlet streams have the same value of pressure.**

For more details, please refer to the IDAES documentation:
https://idaes-pse.readthedocs.io/en/stable

Setting up the problem in IDAES
-------------------------------

In the following cell, we will be importing the necessary components
from Pyomo and IDAES.

.. code:: ipython3

    # Import objects from pyomo package 
    from pyomo.environ import ConcreteModel, SolverFactory, value
    
    # Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
    from idaes.core import FlowsheetBlock
    
    # Import the mixer unit model
    from idaes.generic_models.unit_models import Mixer, MomentumMixingType
    
    # Import idaes logger to set output levels
    import idaes.logger as idaeslog
    
    # Import the BTX_ideal property package to create a properties block for the flowsheet
    from idaes.generic_models.properties.activity_coeff_models import BTX_activity_coeff_VLE
    
    # Import the degrees_of_freedom function from the idaes.core.util.model_statistics package
    # DOF = Number of Model Variables - Number of Model Constraints
    from idaes.core.util.model_statistics import degrees_of_freedom
    
    # Create the ConcreteModel and the FlowsheetBlock objects, and attach the flowsheet block to it.
    m = ConcreteModel()
    
    m.fs = FlowsheetBlock(default={"dynamic": False}) # dynamic or ss flowsheet needs to be specified here
    
    # Add properties parameter block to the flowsheet with specifications
    m.fs.properties = BTX_activity_coeff_VLE.BTXParameterBlock(default={"valid_phase": 'Liq',
                                                                        "activity_coeff_model":
                                                                        "Ideal"})
    

Case 1:
-------

Specify the number of inlets to the mixer, and set the
``momentum_mixing`` type set to “minimize”.

.. code:: ipython3

    # Create an instance of the mixer unit, attaching it to the flowsheet
    # Specify that the property package to be used with the mixer is the
    # one we created earlier, the number of mixer inlets is 2, and momentum
    # mixing type is minimize
    
    m.fs.mixer_1 = Mixer(default={"property_package": m.fs.properties,
                                       "num_inlets":2,
                                       "momentum_mixing_type":MomentumMixingType.minimize})
    
    # Call the degrees_of_freedom function, get initial DOF
    DOF_initial = degrees_of_freedom(m)
    print('The initial degrees of freedom is: {0}'.format(DOF_initial))


.. parsed-literal::

    The initial degrees of freedom is: 10
    

For case 1, we chose to specify only the number of inlets and names were
not specified. When this option is selected, the inlets are named as
“inlet_1”, “inlet_2” and so on depending on the number of inlets
specified. In the following cell, we will use this naming convention to
specify the inlet conditions.

.. code:: ipython3

    # Fix the inlet conditions
    
    # Benzene stream
    m.fs.mixer_1.inlet_1.flow_mol.fix(100) # converting to mol/s as unit basis is mol/s
    m.fs.mixer_1.inlet_1.mole_frac_comp[0, "benzene"].fix(1)
    m.fs.mixer_1.inlet_1.mole_frac_comp[0, "toluene"].fix(0)
    m.fs.mixer_1.inlet_1.pressure.fix(101325) # Pa
    m.fs.mixer_1.inlet_1.temperature.fix(353) # K
    
    # Toluene stream
    m.fs.mixer_1.inlet_2.flow_mol.fix(100) # converting to mol/s as unit basis is mol/s
    m.fs.mixer_1.inlet_2.mole_frac_comp[0, "benzene"].fix(0)
    m.fs.mixer_1.inlet_2.mole_frac_comp[0, "toluene"].fix(1)
    m.fs.mixer_1.inlet_2.pressure.fix(202650) # Pa
    m.fs.mixer_1.inlet_2.temperature.fix(356) # K
    
    # Call the degrees_of_freedom function, get final DOF
    DOF_final = degrees_of_freedom(m)
    print('The final degrees of freedom is: {0}'.format(DOF_final))


.. parsed-literal::

    The final degrees of freedom is: 0
    

Flowsheet Initialization
~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    # Initialize the flowsheet, and set the output at WARNING
    m.fs.mixer_1.initialize(outlvl=idaeslog.WARNING)

Obtaining Simulation Results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    # Solve the simulation using ipopt
    # Note: If the degrees of freedom = 0, we have a square problem
    opt = SolverFactory('ipopt')
    result = opt.solve(m, tee=True)


.. parsed-literal::

    Ipopt 3.13.2: 
    
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
    
    Number of nonzeros in equality constraint Jacobian...:       68
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:       19
    
    Total number of variables............................:       25
                         variables with only lower bounds:        3
                    variables with lower and upper bounds:        8
                         variables with only upper bounds:        0
    Total number of equality constraints.................:       25
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  0.0000000e+00 3.58e+02 1.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  0.0000000e+00 3.58e+00 1.00e-02  -1.0 1.00e-02    -  9.90e-01 9.90e-01h  1
       2  0.0000000e+00 3.54e-02 1.98e-03  -1.0 1.00e-04    -  9.90e-01 9.90e-01h  1
       3  0.0000000e+00 7.28e-12 1.32e+01  -1.0 9.90e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 3
    
                                       (scaled)                 (unscaled)
    Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    Constraint violation....:   2.0968859831870735e-12    7.2759576141834259e-12
    Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    Overall NLP error.......:   2.0968859831870735e-12    7.2759576141834259e-12
    
    
    Number of objective function evaluations             = 4
    Number of objective gradient evaluations             = 4
    Number of equality constraint evaluations            = 4
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 4
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 3
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.003
    Total CPU secs in NLP function evaluations           =      0.000
    
    EXIT: Optimal Solution Found.
    

View Results
~~~~~~~~~~~~

.. code:: ipython3

    # Display output report
    m.fs.mixer_1.report()


.. parsed-literal::

    
    ====================================================================================
    Unit : fs.mixer_1                                                          Time: 0.0
    ------------------------------------------------------------------------------------
        Stream Table
                                inlet_1  inlet_2   Outlet  
        flow_mol                   100      100      200.00
        mole_frac_comp benzene       1        0     0.50000
        mole_frac_comp toluene       0        1     0.50000
        temperature                353      356      354.61
        pressure                101325   202650  1.0133e+05
    ====================================================================================
    

Case 2
------

For case 2, we will specify the inlet names for the two inlets, and set
``momentum_mixing`` type set to “equality” (in this case, pressure will
be specified for only one inlet stream). We will name the 2 inlets as
“benzene_inlet” and “toluene_inlet”.

.. code:: ipython3

    # Create an instance of another mixer unit, attaching it to the same flowsheet. 
    # Specify that the property package to be used with the mixer is the one we created earlier,
    # inlet list is specified but names are specified, and momentum mixing type is equality
    
    m.fs.mixer_2 = Mixer(default={"property_package": m.fs.properties,
                                  "inlet_list":["benzene_inlet","toluene_inlet"],
                                  "momentum_mixing_type":MomentumMixingType.equality})

.. code:: ipython3

    # Check the required degrees of freedom
    DOF_init = degrees_of_freedom(m.fs.mixer_2)
    print('The initial degrees of freedom is: {0}'.format(DOF_init))


.. parsed-literal::

    The initial degrees of freedom is: 9
    

We see that the degrees of freedom has dropped by 1 to 9 when compared
with case 1. This is because we selected the ``momentum_mixing_type`` as
``MomentumMixingType.equality`` which basically adds a constraint that
equates the pressure between all inlets and the outlet. Therefore, when
we specify the inlet confitions in the next cell, we will define the
pressure for only the ``benzene_inlet`` stream.

.. code:: ipython3

    # Fix the stream inlet conditions
    
    # Benzene stream
    m.fs.mixer_2.benzene_inlet.flow_mol.fix(100) # converting to mol/s as unit basis is mol/s
    m.fs.mixer_2.benzene_inlet.mole_frac_comp[0, "benzene"].fix(1)
    m.fs.mixer_2.benzene_inlet.mole_frac_comp[0, "toluene"].fix(0)
    m.fs.mixer_2.benzene_inlet.pressure.fix(101325) # Pa , Another option is m1.fs.mixer2.inlet2.pressure.fix(202650)
    m.fs.mixer_2.benzene_inlet.temperature.fix(353) # K
    
    # Toluene stream
    m.fs.mixer_2.toluene_inlet.flow_mol.fix(100) # converting to mol/s as unit basis is mol/s
    m.fs.mixer_2.toluene_inlet.mole_frac_comp[0, "benzene"].fix(0)
    m.fs.mixer_2.toluene_inlet.mole_frac_comp[0, "toluene"].fix(1)
    m.fs.mixer_2.toluene_inlet.temperature.fix(356) # K
    
    DOF_final = degrees_of_freedom(m.fs.mixer_2)
    print('The final degrees of freedom is: {0}'.format(DOF_final))


.. parsed-literal::

    The final degrees of freedom is: 0
    

Flowsheet Initialization
~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    #Initialize the flowsheet, and set the output at WARNING
    
    m.fs.mixer_2.initialize(outlvl=idaeslog.WARNING)

Obtaining Simulation Results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    # Solve the simulation using ipopt
    # Note: If the degrees of freedom = 0, we have a square problem
    opt = SolverFactory('ipopt')
    result = opt.solve(m.fs.mixer_2, tee=True)


.. parsed-literal::

    Ipopt 3.13.2: 
    
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
    
    Number of nonzeros in equality constraint Jacobian...:       66
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:       18
    
    Total number of variables............................:       24
                         variables with only lower bounds:        4
                    variables with lower and upper bounds:        8
                         variables with only upper bounds:        0
    Total number of equality constraints.................:       24
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  0.0000000e+00 3.58e+02 1.00e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  0.0000000e+00 3.58e+00 1.00e-02  -1.0 1.00e-02    -  9.90e-01 9.90e-01h  1
       2  0.0000000e+00 3.54e-02 1.98e-03  -1.0 1.00e-04    -  9.90e-01 9.90e-01h  1
       3  0.0000000e+00 7.28e-12 1.32e+01  -1.0 9.90e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 3
    
                                       (scaled)                 (unscaled)
    Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
    Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
    Constraint violation....:   2.0968859831870735e-12    7.2759576141834259e-12
    Complementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00
    Overall NLP error.......:   2.0968859831870735e-12    7.2759576141834259e-12
    
    
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
    

View Results
~~~~~~~~~~~~

.. code:: ipython3

    # Display a readable report
    m.fs.mixer_2.report()


.. parsed-literal::

    
    ====================================================================================
    Unit : fs.mixer_2                                                          Time: 0.0
    ------------------------------------------------------------------------------------
        Stream Table
                                benzene_inlet  toluene_inlet   Outlet  
        flow_mol                      100           100.00       200.00
        mole_frac_comp benzene          1           0.0000      0.50000
        mole_frac_comp toluene          0           1.0000      0.50000
        temperature                   353           356.00       354.61
        pressure                   101325       1.0132e+05   1.0132e+05
    ====================================================================================
    

