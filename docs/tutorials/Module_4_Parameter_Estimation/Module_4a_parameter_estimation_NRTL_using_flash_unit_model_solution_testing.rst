Module 4a: Parameter Estimation Using Flash Unit Model
------------------------------------------------------

In this module, we will be using Pyomo’s ``parmest`` tool in conjuction
with IDAES models for parameter estimation. We demonstrate these tools
by estimating the parameters associated with the NRTL property model for
a benzene-toluene mixture. The NRTL model has 2 sets of parameters: the
non-randomness parameter (``alpha_ij``) and the binary interaction
parameter (``tau_ij``), where ``i`` and ``j`` is the pure component
species. In this example, we will be only estimate the binary
interaction parameter (``tau_ij``) for a given dataset. When estimating
parameters associated with the property package, IDAES provides the
flexibility of doing the parameter estimation by just using the state
block or by using a unit model with a specified property package. This
module will demonstrate parameter estimation by using the flash unit
model with the NRTL property package.

We will complete the following tasks: \* Set up a method to return an
initialized model \* Set up the parameter estimation problem using
``parmest`` \* Analyze the results \* Demonstrate advanced features from
``parmest``

Key links to documentation:
---------------------------

-  NRTL Model -
   https://idaes-pse.readthedocs.io/en/latest/model_libraries/core_library/property_models/activity_coefficient.html
-  parmest -
   https://pyomo.readthedocs.io/en/stable/contributed_packages/parmest/index.html

.. raw:: html

   <div class="alert alert-block alert-info">

Inline Exercise: import ``ConcreteModel`` from Pyomo, ``FlowsheetBlock``
and ``Flash`` from IDAES.

.. raw:: html

   </div>

.. code:: ipython3

    # Todo: import ConcreteModel from pyomo.environ
    from pyomo.environ import ConcreteModel, value
    
    # Todo: import FlowsheetBlock from idaes.core
    from idaes.core import FlowsheetBlock
    
    # Todo: import Flash unit model from idaes.generic_models.unit_models
    from idaes.generic_models.unit_models import Flash
    

In the next cell, we will be importing the parameter block that we will
be using in this module and the idaes logger.

.. code:: ipython3

    from idaes.generic_models.properties.activity_coeff_models.\
        BTX_activity_coeff_VLE import BTXParameterBlock
    import idaes.logger as idaeslog

In the next cell, we import ``parmest`` from Pyomo and the ``pandas``
package. We need ``pandas`` as ``parmest`` uses ``pandas.dataframe`` for
handling the input data and the results.

.. code:: ipython3

    import pyomo.contrib.parmest.parmest as parmest
    import pandas as pd

Setting up an initialized model
-------------------------------

We need to provide a method that returns an initialized model to the
``parmest`` tool in Pyomo.

.. raw:: html

   <div class="alert alert-block alert-info">

Inline Exercise: Using what you have learned from previous modules, fill
in the missing code below to return an initialized IDAES model.

.. raw:: html

   </div>

.. code:: ipython3

    def NRTL_model(data):
        
        #Todo: Create a ConcreteModel object
        m = ConcreteModel()
        
        #Todo: Create FlowsheetBlock object
        m.fs = FlowsheetBlock(default={"dynamic": False})
        
    
        #Todo: Create a properties parameter object with the following options:
        # "valid_phase": ('Liq', 'Vap')
        # "activity_coeff_model": 'NRTL'
        m.fs.properties = BTXParameterBlock(default={"valid_phase":
                                                     ('Liq', 'Vap'),
                                                     "activity_coeff_model":
                                                     'NRTL'})
        m.fs.flash = Flash(default={"property_package": m.fs.properties})
    
        # Initialize at a certain inlet condition
        m.fs.flash.inlet.flow_mol.fix(1)
        m.fs.flash.inlet.temperature.fix(368)
        m.fs.flash.inlet.pressure.fix(101325)
        m.fs.flash.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.flash.inlet.mole_frac_comp[0, "toluene"].fix(0.5)
    
        # Set Flash unit specifications
        m.fs.flash.heat_duty.fix(0)
        m.fs.flash.deltaP.fix(0)
    
        # Fix NRTL specific variables
        # alpha values (set at 0.3)
        m.fs.properties.\
            alpha["benzene", "benzene"].fix(0)
        m.fs.properties.\
            alpha["benzene", "toluene"].fix(0.3)
        m.fs.properties.\
            alpha["toluene", "toluene"].fix(0)
        m.fs.properties.\
            alpha["toluene", "benzene"].fix(0.3)
    
        # initial tau values
        m.fs.properties.\
            tau["benzene", "benzene"].fix(0)
        m.fs.properties.\
            tau["benzene", "toluene"].fix(0.1690)
        m.fs.properties.\
            tau["toluene", "toluene"].fix(0)
        m.fs.properties.\
            tau["toluene", "benzene"].fix(-0.1559)
    
        # Initialize the flash unit
        m.fs.flash.initialize(outlvl=idaeslog.INFO_LOW)
    
        # Fix at actual temperature
        m.fs.flash.inlet.temperature.fix(float(data["temperature"]))
    
        # Set bounds on variables to be estimated
        m.fs.properties.\
            tau["benzene", "toluene"].setlb(-5)
        m.fs.properties.\
            tau["benzene", "toluene"].setub(5)
    
        m.fs.properties.\
            tau["toluene", "benzene"].setlb(-5)
        m.fs.properties.\
            tau["toluene", "benzene"].setub(5)
    
        # Return initialized flash model
        return m
    

.. code:: ipython3

    from idaes.core.util.model_statistics import degrees_of_freedom
    import pytest
    
    # Testing the initialized model
    test_data = {"temperature": 368}
    
    m = NRTL_model(test_data)
    
    # Check that degrees of freedom is 0
    assert degrees_of_freedom(m) == 0
    
    # Check for output values
    assert value(m.fs.flash.liq_outlet.mole_frac_comp[0, 'benzene']) == pytest.approx(0.4105, abs=1e-3)
    assert value(m.fs.flash.vap_outlet.mole_frac_comp[0, 'benzene']) == pytest.approx(0.6326, abs=1e-3)
    
    assert value(m.fs.flash.liq_outlet.mole_frac_comp[0, 'toluene']) == pytest.approx(0.5895, abs=1e-3)
    assert value(m.fs.flash.vap_outlet.mole_frac_comp[0, 'toluene']) == pytest.approx(0.3673, abs=1e-3)

Parameter estimation using parmest
----------------------------------

In addition to providing a method to return an initialized model, the
``parmest`` tool needs the following:

-  List of variable names to be estimated
-  Dataset with multiple scenarios
-  Expression to compute the sum of squared errors

In this example, we only estimate the binary interaction parameter
(``tau_ij``). Given that this variable is usually indexed as
``tau_ij = Var(component_list, component_list)``, there are 2*2=4
degrees of freedom. However, when i=j, the binary interaction parameter
is 0. Therefore, in this problem, we estimate the binary interaction
parameter for the following variables only:

-  fs.properties.tau[‘benzene’, ‘toluene’]
-  fs.properties.tau[‘toluene’, ‘benzene’]

.. raw:: html

   <div class="alert alert-block alert-info">

Inline Exercise: Create a list called ``variable_name`` with the
above-mentioned variables declared as strings.

.. raw:: html

   </div>

.. code:: ipython3

    # Todo: Create a list of vars to estimate
    variable_name = ["fs.properties.tau['benzene', 'toluene']",
                     "fs.properties.tau['toluene', 'benzene']"]
    

Pyomo’s ``parmest`` tool supports the following data formats: - pandas
dataframe - list of dictionaries - list of json file names.

Please see the documentation for more details.

For this example, we load data from the csv file
``BT_NRTL_dataset.csv``. The dataset consists of fifty data points which
provide the mole fraction of benzene in the vapor and liquid phase as a
function of temperature.

.. code:: ipython3

    # Load data from csv
    data = pd.read_csv('BT_NRTL_dataset.csv')
    
    # Display the dataset
    display(data)



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
          <th>temperature</th>
          <th>liq_benzene</th>
          <th>vap_benzene</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>365.500000</td>
          <td>0.490769</td>
          <td>0.706235</td>
        </tr>
        <tr>
          <th>1</th>
          <td>365.617647</td>
          <td>0.486783</td>
          <td>0.702841</td>
        </tr>
        <tr>
          <th>2</th>
          <td>365.735294</td>
          <td>0.482812</td>
          <td>0.699436</td>
        </tr>
        <tr>
          <th>3</th>
          <td>365.852941</td>
          <td>0.478855</td>
          <td>0.696018</td>
        </tr>
        <tr>
          <th>4</th>
          <td>365.970588</td>
          <td>0.474912</td>
          <td>0.692587</td>
        </tr>
        <tr>
          <th>5</th>
          <td>366.088235</td>
          <td>0.470984</td>
          <td>0.689144</td>
        </tr>
        <tr>
          <th>6</th>
          <td>366.205882</td>
          <td>0.467069</td>
          <td>0.685689</td>
        </tr>
        <tr>
          <th>7</th>
          <td>366.323529</td>
          <td>0.463169</td>
          <td>0.682221</td>
        </tr>
        <tr>
          <th>8</th>
          <td>366.441177</td>
          <td>0.459282</td>
          <td>0.678741</td>
        </tr>
        <tr>
          <th>9</th>
          <td>366.558823</td>
          <td>0.455409</td>
          <td>0.675248</td>
        </tr>
        <tr>
          <th>10</th>
          <td>366.676471</td>
          <td>0.451550</td>
          <td>0.671743</td>
        </tr>
        <tr>
          <th>11</th>
          <td>366.794118</td>
          <td>0.447705</td>
          <td>0.668225</td>
        </tr>
        <tr>
          <th>12</th>
          <td>366.911765</td>
          <td>0.443873</td>
          <td>0.664694</td>
        </tr>
        <tr>
          <th>13</th>
          <td>367.029412</td>
          <td>0.440055</td>
          <td>0.661151</td>
        </tr>
        <tr>
          <th>14</th>
          <td>367.147059</td>
          <td>0.436250</td>
          <td>0.657595</td>
        </tr>
        <tr>
          <th>15</th>
          <td>367.264706</td>
          <td>0.432459</td>
          <td>0.654025</td>
        </tr>
        <tr>
          <th>16</th>
          <td>367.382353</td>
          <td>0.428681</td>
          <td>0.650444</td>
        </tr>
        <tr>
          <th>17</th>
          <td>367.500000</td>
          <td>0.424916</td>
          <td>0.646849</td>
        </tr>
        <tr>
          <th>18</th>
          <td>367.617647</td>
          <td>0.421164</td>
          <td>0.643241</td>
        </tr>
        <tr>
          <th>19</th>
          <td>367.735294</td>
          <td>0.417426</td>
          <td>0.639620</td>
        </tr>
        <tr>
          <th>20</th>
          <td>367.852941</td>
          <td>0.413700</td>
          <td>0.635986</td>
        </tr>
        <tr>
          <th>21</th>
          <td>367.970588</td>
          <td>0.409987</td>
          <td>0.632339</td>
        </tr>
        <tr>
          <th>22</th>
          <td>368.000000</td>
          <td>0.409061</td>
          <td>0.631426</td>
        </tr>
        <tr>
          <th>23</th>
          <td>368.088235</td>
          <td>0.406287</td>
          <td>0.628679</td>
        </tr>
        <tr>
          <th>24</th>
          <td>368.205882</td>
          <td>0.402600</td>
          <td>0.625006</td>
        </tr>
        <tr>
          <th>25</th>
          <td>368.323529</td>
          <td>0.398926</td>
          <td>0.621320</td>
        </tr>
        <tr>
          <th>26</th>
          <td>368.441177</td>
          <td>0.395264</td>
          <td>0.617620</td>
        </tr>
        <tr>
          <th>27</th>
          <td>368.558823</td>
          <td>0.391615</td>
          <td>0.613907</td>
        </tr>
        <tr>
          <th>28</th>
          <td>368.676471</td>
          <td>0.387978</td>
          <td>0.610180</td>
        </tr>
        <tr>
          <th>29</th>
          <td>368.794118</td>
          <td>0.384353</td>
          <td>0.606440</td>
        </tr>
        <tr>
          <th>30</th>
          <td>368.911765</td>
          <td>0.380741</td>
          <td>0.602687</td>
        </tr>
        <tr>
          <th>31</th>
          <td>369.029412</td>
          <td>0.377141</td>
          <td>0.598920</td>
        </tr>
        <tr>
          <th>32</th>
          <td>369.147059</td>
          <td>0.373553</td>
          <td>0.595140</td>
        </tr>
        <tr>
          <th>33</th>
          <td>369.264706</td>
          <td>0.369978</td>
          <td>0.591346</td>
        </tr>
        <tr>
          <th>34</th>
          <td>369.382353</td>
          <td>0.366414</td>
          <td>0.587538</td>
        </tr>
        <tr>
          <th>35</th>
          <td>369.500000</td>
          <td>0.362862</td>
          <td>0.583717</td>
        </tr>
        <tr>
          <th>36</th>
          <td>369.617647</td>
          <td>0.359323</td>
          <td>0.579882</td>
        </tr>
        <tr>
          <th>37</th>
          <td>369.735294</td>
          <td>0.355795</td>
          <td>0.576033</td>
        </tr>
        <tr>
          <th>38</th>
          <td>369.852941</td>
          <td>0.352278</td>
          <td>0.572171</td>
        </tr>
        <tr>
          <th>39</th>
          <td>369.970588</td>
          <td>0.348774</td>
          <td>0.568294</td>
        </tr>
        <tr>
          <th>40</th>
          <td>370.088235</td>
          <td>0.345281</td>
          <td>0.564404</td>
        </tr>
        <tr>
          <th>41</th>
          <td>370.205882</td>
          <td>0.341799</td>
          <td>0.560500</td>
        </tr>
        <tr>
          <th>42</th>
          <td>370.323529</td>
          <td>0.338329</td>
          <td>0.556581</td>
        </tr>
        <tr>
          <th>43</th>
          <td>370.441177</td>
          <td>0.334871</td>
          <td>0.552649</td>
        </tr>
        <tr>
          <th>44</th>
          <td>370.558823</td>
          <td>0.331423</td>
          <td>0.548702</td>
        </tr>
        <tr>
          <th>45</th>
          <td>370.676471</td>
          <td>0.327987</td>
          <td>0.544741</td>
        </tr>
        <tr>
          <th>46</th>
          <td>370.794118</td>
          <td>0.324563</td>
          <td>0.540766</td>
        </tr>
        <tr>
          <th>47</th>
          <td>370.911765</td>
          <td>0.321149</td>
          <td>0.536777</td>
        </tr>
        <tr>
          <th>48</th>
          <td>371.029412</td>
          <td>0.317746</td>
          <td>0.532774</td>
        </tr>
        <tr>
          <th>49</th>
          <td>371.147059</td>
          <td>0.314354</td>
          <td>0.528756</td>
        </tr>
      </tbody>
    </table>
    </div>


We need to provide a method to return an expression to compute the sum
of squared errors that will be used as the objective in solving the
parameter estimation problem. For this problem, the error will be
computed for the mole fraction of benzene in the vapor and liquid phase
between the model prediction and data.

.. raw:: html

   <div class="alert alert-block alert-info">

Inline Exercise: Complete the following cell by adding an expression to
compute the sum of square errors.

.. raw:: html

   </div>

.. code:: ipython3

    # Create method to return an expression that computes the sum of squared error
    def SSE(m, data):
        # Todo: Add expression for computing the sum of squared errors in mole fraction of benzene in the liquid
        # and vapor phase. For example, the squared error for the vapor phase is:
        # (float(data["vap_benzene"]) - m.fs.flash.vap_outlet.mole_frac_comp[0, "benzene"])**2
        expr = ((float(data["vap_benzene"]) -
                 m.fs.flash.vap_outlet.mole_frac_comp[0, "benzene"])**2 +
                (float(data["liq_benzene"]) -
                 m.fs.flash.liq_outlet.mole_frac_comp[0, "benzene"])**2)
        return expr*1E4

.. raw:: html

   <div class="alert alert-block alert-warning">

Note: Notice that we have scaled the expression up by a factor of 10000
as the SSE computed here will be an extremely small number given that we
are using the difference in mole fraction in our expression. A
well-scaled objective will help improve solve robustness when using
IPOPT.

.. raw:: html

   </div>

We are now ready to set up the parameter estimation problem. We will
create a parameter estimation object called ``pest``. As shown below, we
pass the method that returns an initialized model, dataset, list of
variable names to estimate, and the SSE expression to the Estimator
object. ``tee=True`` will print the solver output after solving the
parameter estimation problem.

.. code:: ipython3

    # Initialize a parameter estimation object
    pest = parmest.Estimator(NRTL_model, data, variable_name, SSE, tee=True)
    
    # Run parameter estimation using all data
    obj_value, parameters = pest.theta_est()


.. parsed-literal::

    Ipopt 3.12.13: max_iter=6000
    
    
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
    
    Number of nonzeros in equality constraint Jacobian...:    10950
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:     6000
    
    Total number of variables............................:     2952
                         variables with only lower bounds:      150
                    variables with lower and upper bounds:      600
                         variables with only upper bounds:        0
    Total number of equality constraints.................:     2950
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  5.5857491e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  8.0006441e-02 1.79e+03 1.38e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  7.2326943e-03 5.70e+02 1.39e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  6.6862864e-03 1.29e+02 2.40e-03  -1.0 2.48e+00    -  9.95e-01 1.00e+00h  1
       4  6.4759005e-03 6.81e-02 5.12e-04  -1.7 9.86e-01    -  1.00e+00 1.00e+00h  1
       5  6.4501818e-03 1.27e+00 2.93e-04  -3.8 1.58e-01    -  1.00e+00 1.00e+00h  1
       6  6.4443106e-03 5.22e-02 1.01e-04  -3.8 3.15e-02  -4.0 1.00e+00 1.00e+00h  1
       7  4.9202577e-03 9.26e+03 5.32e-02  -5.7 1.34e+01    -  8.28e-01 1.00e+00h  1
       8  5.6550947e-03 5.14e+03 2.96e-02  -5.7 7.74e+00    -  9.04e-01 5.00e-01h  2
       9  5.2465884e-03 1.15e+03 1.08e-02  -5.7 4.67e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.6784682e-03 5.72e-01 7.56e-04  -5.7 5.53e+00    -  1.00e+00 1.00e+00h  1
      11  4.6633539e-03 5.56e-02 3.71e-06  -5.7 5.02e-01    -  1.00e+00 1.00e+00h  1
      12  4.6633491e-03 3.70e-05 1.06e-09  -5.7 6.32e-03    -  1.00e+00 1.00e+00h  1
      13  4.6633488e-03 3.10e-05 1.34e-09  -8.6 7.46e-03    -  1.00e+00 1.00e+00h  1
      14  4.6633488e-03 7.45e-09 2.51e-14  -8.6 6.81e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 14
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.6633488370399628e-03    4.6633488370399628e-03
    Dual infeasibility......:   2.5085432767634323e-14    2.5085432767634323e-14
    Constraint violation....:   1.4104644499482428e-11    7.4505805969238281e-09
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    7.4505805969238281e-09
    
    
    Number of objective function evaluations             = 16
    Number of objective gradient evaluations             = 15
    Number of equality constraint evaluations            = 16
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 15
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 14
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.043
    Total CPU secs in NLP function evaluations           =      0.051
    
    EXIT: Optimal Solution Found.
    

.. code:: ipython3

    # Check for values of the parameter estimation problem
    assert obj_value == pytest.approx(0.004663, 1e-3)
    assert parameters["fs.properties.tau[('benzene', 'toluene')]"] == pytest.approx(0.47811, 1e-3) 
    assert parameters["fs.properties.tau[('toluene', 'benzene')]"] == pytest.approx(-0.40924, 1e-3)

You will notice that the resulting parameter estimation problem, when
using the flash unit model, will have 2952 variables and 2950
constraints. This is because the unit models in IDAES use control volume
blocks which have two state blocks attached; one at the inlet and one at
the outlet. Even though there are two state blocks, they still use the
same parameter block i.e. ``m.fs.properties`` in our example which is
where our parameters that need to be estimated exist.

Let us display the results by running the next cell.

.. code:: ipython3

    print("The SSE at the optimal solution is %0.6f" % obj_value)
    print()
    print("The values for the parameters are as follows:")
    for k,v in parameters.items():
        print(k, "=", v)


.. parsed-literal::

    The SSE at the optimal solution is 0.004663
    
    The values for the parameters are as follows:
    fs.properties.tau[('benzene', 'toluene')] = 0.47810867841010257
    fs.properties.tau[('toluene', 'benzene')] = -0.4092446537759336
    

Using the data that was provided, we have estimated the binary
interaction parameters in the NRTL model for a benzene-toluene mixture.
Although the dataset that was provided was temperature dependent, in
this example we have estimated a single value that fits best for all
temperatures.

Advanced options for parmest: bootstrapping
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pyomo’s ``parmest`` tool allows for bootstrapping where the parameter
estimation is repeated over ``n`` samples with resampling from the
original data set. Parameter estimation with bootstrap resampling can be
used to identify confidence regions around each parameter estimate. This
analysis can be slow given the increased number of model instances that
need to be solved. In the following cell, we run the parameter
estimation with 10 bootstrap samples from the given dataset. We then
plot the parameter estimates along with an confidence regions using
rectangular and multivariate normal distributions.

.. code:: ipython3

    # Run parameter estimation using bootstrap resample of the data (10 samples),
    # plot results along with confidence regions
    bootstrap_theta = pest.theta_est_bootstrap(10)


.. parsed-literal::

    Ipopt 3.12.13: max_iter=6000
    
    
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
    
    Number of nonzeros in equality constraint Jacobian...:    10950
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:     6000
    
    Total number of variables............................:     2952
                         variables with only lower bounds:      150
                    variables with lower and upper bounds:      600
                         variables with only upper bounds:        0
    Total number of equality constraints.................:     2950
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  6.4113312e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  8.8129365e-02 1.79e+03 1.58e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  8.7956558e-03 3.85e+02 1.25e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  8.1371059e-03 4.37e+01 1.27e-03  -1.0 2.40e+00    -  9.95e-01 1.00e+00h  1
       4  7.8355611e-03 3.83e-02 4.59e-04  -2.5 3.42e+00    -  1.00e+00 1.00e+00h  1
       5  7.4944400e-03 1.27e+02 1.23e-03  -3.8 1.50e+00    -  1.00e+00 1.00e+00h  1
       6  7.4839166e-03 4.14e+00 1.26e-03  -3.8 4.44e-02  -4.0 1.00e+00 1.00e+00h  1
       7  5.9486770e-03 1.78e+03 1.99e-02  -3.8 8.87e+00    -  1.00e+00 1.00e+00h  1
       8  5.2416399e-03 1.64e+02 4.46e-03  -3.8 1.04e+01    -  1.00e+00 1.00e+00h  1
       9  5.1602419e-03 6.38e-01 5.00e-05  -3.8 1.47e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  5.1600918e-03 7.79e-03 1.35e-07  -3.8 5.30e-02    -  1.00e+00 1.00e+00h  1
      11  5.1587813e-03 1.33e-01 7.15e-06  -5.7 4.97e-01    -  1.00e+00 1.00e+00h  1
      12  5.1587627e-03 2.46e-05 4.29e-10  -5.7 2.90e-03    -  1.00e+00 1.00e+00h  1
      13  5.1587625e-03 1.94e-05 1.06e-09  -8.6 6.09e-03    -  1.00e+00 1.00e+00h  1
      14  5.1587625e-03 7.45e-09 2.51e-14  -8.6 4.28e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 14
    
                                       (scaled)                 (unscaled)
    Objective...............:   5.1587624626448459e-03    5.1587624626448459e-03
    Dual infeasibility......:   2.5085432767634323e-14    2.5085432767634323e-14
    Constraint violation....:   1.4104644499482428e-11    7.4505805969238281e-09
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    7.4505805969238281e-09
    
    
    Number of objective function evaluations             = 15
    Number of objective gradient evaluations             = 15
    Number of equality constraint evaluations            = 15
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 15
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 14
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.048
    Total CPU secs in NLP function evaluations           =      0.082
    
    EXIT: Optimal Solution Found.
    Ipopt 3.12.13: max_iter=6000
    
    
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
    
    Number of nonzeros in equality constraint Jacobian...:    10950
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:     6000
    
    Total number of variables............................:     2952
                         variables with only lower bounds:      150
                    variables with lower and upper bounds:      600
                         variables with only upper bounds:        0
    Total number of equality constraints.................:     2950
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  5.7021605e+01 5.20e+05 3.84e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  7.0218344e-02 1.53e+03 1.49e+00  -1.0 1.27e+04    -  9.82e-01 1.00e+00h  1
       2  7.9578832e-03 2.72e+02 1.00e-01  -1.0 4.05e+02    -  9.90e-01 1.00e+00h  1
       3  7.4460193e-03 7.75e+00 4.44e-04  -1.0 1.97e+00    -  9.96e-01 1.00e+00h  1
       4  7.0005419e-03 7.53e-02 4.13e-04  -2.5 5.26e+00    -  1.00e+00 1.00e+00h  1
       5  6.7277182e-03 9.41e+01 8.93e-04  -3.8 1.17e+00    -  1.00e+00 1.00e+00h  1
       6  6.7188663e-03 1.47e+00 6.85e-04  -3.8 4.00e-02  -4.0 1.00e+00 1.00e+00h  1
       7  5.4239656e-03 1.79e+03 2.01e-02  -3.8 7.85e+00    -  1.00e+00 1.00e+00h  1
       8  4.8111880e-03 1.65e+02 4.28e-03  -3.8 9.56e+00    -  1.00e+00 1.00e+00h  1
       9  4.7459397e-03 5.78e-01 3.81e-05  -3.8 1.14e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.7458001e-03 5.83e-03 6.80e-08  -3.8 3.33e-02    -  1.00e+00 1.00e+00h  1
      11  4.7443696e-03 1.68e-01 8.25e-06  -5.7 5.02e-01    -  1.00e+00 1.00e+00h  1
      12  4.7443510e-03 4.03e-05 6.33e-10  -5.7 3.32e-03    -  1.00e+00 1.00e+00h  1
      13  4.7443508e-03 2.42e-05 1.22e-09  -8.6 6.15e-03    -  1.00e+00 1.00e+00h  1
      14  4.7443508e-03 1.49e-08 2.51e-14  -8.6 4.88e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 14
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.7443508083237001e-03    4.7443508083237001e-03
    Dual infeasibility......:   2.5085432767634323e-14    2.5085432767634323e-14
    Constraint violation....:   1.4104644499482428e-11    1.4901161193847656e-08
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    1.4901161193847656e-08
    
    
    Number of objective function evaluations             = 15
    Number of objective gradient evaluations             = 15
    Number of equality constraint evaluations            = 15
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 15
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 14
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.043
    Total CPU secs in NLP function evaluations           =      0.052
    
    EXIT: Optimal Solution Found.
    Ipopt 3.12.13: max_iter=6000
    
    
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
    
    Number of nonzeros in equality constraint Jacobian...:    10950
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:     6000
    
    Total number of variables............................:     2952
                         variables with only lower bounds:      150
                    variables with lower and upper bounds:      600
                         variables with only upper bounds:        0
    Total number of equality constraints.................:     2950
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  5.1307463e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  7.3366710e-02 1.79e+03 1.41e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  6.6506372e-03 5.56e+02 1.39e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  6.1467594e-03 1.18e+02 2.29e-03  -1.0 2.43e+00    -  9.95e-01 1.00e+00h  1
       4  5.9581355e-03 5.18e-02 4.85e-04  -1.7 1.14e+00    -  1.00e+00 1.00e+00h  1
       5  5.9360594e-03 1.08e+00 2.39e-04  -3.8 1.49e-01    -  1.00e+00 1.00e+00h  1
       6  5.9310010e-03 4.75e-02 8.66e-05  -3.8 2.78e-02  -4.0 1.00e+00 1.00e+00h  1
       7  4.6470908e-03 7.47e+03 3.88e-02  -5.7 1.14e+01    -  8.44e-01 1.00e+00h  1
       8  5.4564702e-03 1.71e+03 1.28e-02  -5.7 5.93e+00    -  9.26e-01 1.00e+00h  1
       9  4.4559748e-03 1.04e+02 1.97e-03  -5.7 8.47e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.4007638e-03 2.85e+00 4.93e-05  -5.7 7.39e-01    -  1.00e+00 1.00e+00h  1
      11  4.4004673e-03 5.83e-04 5.78e-08  -5.7 9.74e-02    -  1.00e+00 1.00e+00h  1
      12  4.4004667e-03 3.74e-05 1.54e-09  -8.6 7.36e-03    -  1.00e+00 1.00e+00h  1
      13  4.4004667e-03 7.45e-09 2.51e-14  -8.6 7.50e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.4004667292155686e-03    4.4004667292155686e-03
    Dual infeasibility......:   2.5085432767634323e-14    2.5085432767634323e-14
    Constraint violation....:   1.4104644499482428e-11    7.4505805969238281e-09
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    7.4505805969238281e-09
    
    
    Number of objective function evaluations             = 14
    Number of objective gradient evaluations             = 14
    Number of equality constraint evaluations            = 14
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 14
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 13
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.056
    Total CPU secs in NLP function evaluations           =      0.054
    
    EXIT: Optimal Solution Found.
    Ipopt 3.12.13: max_iter=6000
    
    
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
    
    Number of nonzeros in equality constraint Jacobian...:    10950
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:     6000
    
    Total number of variables............................:     2952
                         variables with only lower bounds:      150
                    variables with lower and upper bounds:      600
                         variables with only upper bounds:        0
    Total number of equality constraints.................:     2950
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  5.4917004e+01 5.41e+05 4.00e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  7.5308186e-02 1.66e+03 1.37e+00  -1.0 1.32e+04    -  9.82e-01 1.00e+00h  1
       2  7.3667672e-03 5.21e+02 1.31e-01  -1.0 4.39e+02    -  9.90e-01 1.00e+00h  1
       3  6.7766948e-03 1.02e+02 2.02e-03  -1.0 2.54e+00    -  9.96e-01 1.00e+00h  1
       4  6.5541391e-03 2.90e-02 4.68e-04  -1.7 1.69e+00    -  1.00e+00 1.00e+00h  1
       5  6.5238503e-03 1.40e+00 2.77e-04  -3.8 3.29e-01    -  1.00e+00 1.00e+00h  1
       6  6.5176806e-03 5.43e-02 1.01e-04  -3.8 3.04e-02  -4.0 1.00e+00 1.00e+00h  1
       7  4.9143017e-03 9.45e+03 4.91e-02  -5.7 1.31e+01    -  8.26e-01 1.00e+00h  1
       8  5.6376525e-03 5.24e+03 2.74e-02  -5.7 7.89e+00    -  9.02e-01 5.00e-01h  2
       9  5.2370045e-03 1.18e+03 1.01e-02  -5.7 4.39e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.6731761e-03 1.18e+00 6.69e-04  -5.7 4.94e+00    -  1.00e+00 1.00e+00h  1
      11  4.6578807e-03 2.81e-02 2.51e-06  -5.7 4.58e-01    -  1.00e+00 1.00e+00h  1
      12  4.6578802e-03 1.48e-05 4.67e-10  -5.7 4.06e-03    -  1.00e+00 1.00e+00h  1
      13  4.6578799e-03 2.95e-05 1.31e-09  -8.6 6.93e-03    -  1.00e+00 1.00e+00h  1
      14  4.6578799e-03 1.49e-08 2.51e-14  -8.6 6.15e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 14
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.6578799257086089e-03    4.6578799257086089e-03
    Dual infeasibility......:   2.5140943918865581e-14    2.5140943918865581e-14
    Constraint violation....:   1.4104644499482428e-11    1.4901161193847656e-08
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    1.4901161193847656e-08
    
    
    Number of objective function evaluations             = 16
    Number of objective gradient evaluations             = 15
    Number of equality constraint evaluations            = 16
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 15
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 14
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.045
    Total CPU secs in NLP function evaluations           =      0.054
    
    EXIT: Optimal Solution Found.
    Ipopt 3.12.13: max_iter=6000
    
    
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
    
    Number of nonzeros in equality constraint Jacobian...:    10950
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:     6000
    
    Total number of variables............................:     2952
                         variables with only lower bounds:      150
                    variables with lower and upper bounds:      600
                         variables with only upper bounds:        0
    Total number of equality constraints.................:     2950
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  6.4847156e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  9.3365252e-02 1.79e+03 1.39e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  8.9854203e-03 8.99e+02 1.75e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  7.9665258e-03 3.87e+02 5.45e-03  -1.0 3.77e+00    -  9.95e-01 1.00e+00h  1
       4  7.5848257e-03 1.54e+01 2.05e-04  -1.0 1.70e+00    -  1.00e+00 1.00e+00h  1
       5  7.6556038e-03 6.12e-02 2.87e-04  -2.5 8.02e-01    -  1.00e+00 1.00e+00h  1
       6  7.4007874e-03 8.34e+01 9.45e-04  -3.8 1.30e+00    -  1.00e+00 1.00e+00h  1
       7  7.3909784e-03 9.38e-01 6.65e-04  -3.8 4.57e-02  -4.0 1.00e+00 1.00e+00h  1
       8  6.0059589e-03 2.08e+03 2.13e-02  -3.8 8.55e+00    -  1.00e+00 1.00e+00h  1
       9  5.2598390e-03 2.10e+02 5.01e-03  -3.8 1.16e+01    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  5.1588405e-03 7.22e-01 5.07e-05  -3.8 1.61e+00    -  1.00e+00 1.00e+00h  1
      11  5.1585375e-03 1.16e-02 1.46e-07  -3.8 5.63e-02    -  1.00e+00 1.00e+00h  1
      12  5.1571783e-03 1.47e-01 7.12e-06  -5.7 5.48e-01    -  1.00e+00 1.00e+00h  1
      13  5.1571573e-03 3.07e-05 4.89e-10  -5.7 3.37e-03    -  1.00e+00 1.00e+00h  1
      14  5.1571571e-03 2.14e-05 1.06e-09  -8.6 6.72e-03    -  1.00e+00 1.00e+00h  1
      15  5.1571571e-03 7.45e-09 2.51e-14  -8.6 4.96e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 15
    
                                       (scaled)                 (unscaled)
    Objective...............:   5.1571571053064024e-03    5.1571571053064024e-03
    Dual infeasibility......:   2.5140943918865581e-14    2.5140943918865581e-14
    Constraint violation....:   1.4104644499482428e-11    7.4505805969238281e-09
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    7.4505805969238281e-09
    
    
    Number of objective function evaluations             = 16
    Number of objective gradient evaluations             = 16
    Number of equality constraint evaluations            = 16
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 16
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 15
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.045
    Total CPU secs in NLP function evaluations           =      0.078
    
    EXIT: Optimal Solution Found.
    Ipopt 3.12.13: max_iter=6000
    
    
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
    
    Number of nonzeros in equality constraint Jacobian...:    10950
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:     6000
    
    Total number of variables............................:     2952
                         variables with only lower bounds:      150
                    variables with lower and upper bounds:      600
                         variables with only upper bounds:        0
    Total number of equality constraints.................:     2950
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  6.6613428e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  9.7953195e-02 1.79e+03 1.36e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  1.0123950e-02 1.83e+03 2.47e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  7.4368602e-03 1.46e+03 2.26e-02  -1.0 9.52e+00    -  9.95e-01 1.00e+00h  1
       4  7.0754661e-03 1.05e+02 1.06e-03  -1.0 1.51e+00    -  1.00e+00 1.00e+00h  1
       5  7.7918612e-03 4.21e-03 3.77e-04  -1.7 3.59e+00    -  1.00e+00 1.00e+00h  1
       6  7.7517255e-03 1.86e+00 3.46e-04  -3.8 2.13e-01    -  1.00e+00 1.00e+00h  1
       7  5.6647899e-03 1.17e+04 7.17e-02  -3.8 5.01e+02    -  1.36e-01 3.08e-02h  3
       8  6.9641788e-03 8.98e+03 5.58e-02  -3.8 3.53e+01    -  1.00e+00 2.50e-01h  3
       9  7.0208099e-03 8.84e+03 5.49e-02  -3.8 2.44e+01    -  1.00e+00 1.56e-02h  7
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  7.0487290e-03 8.77e+03 5.45e-02  -3.8 2.39e+01    -  1.00e+00 7.81e-03h  8
      11  7.0504634e-03 8.76e+03 5.44e-02  -3.8 2.38e+01    -  1.00e+00 4.88e-04h 12
      12  7.0508968e-03 8.76e+03 5.44e-02  -3.8 2.38e+01    -  1.00e+00 1.22e-04h 14
      13  1.2279886e-02 2.82e+03 2.77e-02  -3.8 2.37e+01    -  1.00e+00 1.00e+00h  1
      14  1.1777947e-02 5.78e+03 2.13e+00  -3.8 5.78e+00  -4.0 1.00e+00 5.00e-01h  2
      15  1.2254597e-02 5.69e+03 2.10e+00  -3.8 1.40e+02    -  1.00e+00 1.56e-02h  7
      16  1.2292058e-02 5.69e+03 2.10e+00  -3.8 1.55e+02    -  1.00e+00 9.77e-04h 11
      17  1.2311151e-02 5.68e+03 2.10e+00  -3.8 1.56e+02    -  1.00e+00 4.88e-04h 12
      18  1.2312355e-02 5.68e+03 2.10e+00  -3.8 1.57e+02    -  1.00e+00 3.05e-05h 16
      19  1.2312957e-02 5.68e+03 2.10e+00  -3.8 1.57e+02    -  1.00e+00 1.53e-05h 17
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      20  3.2523863e-02 5.05e+03 1.93e+00  -3.8 1.57e+02    -  1.00e+00 2.50e-01h  3
      21  3.6599411e-02 8.00e+02 8.60e-01  -3.8 6.24e+00  -4.5 1.00e+00 1.00e+00h  1
      22  3.7088953e-02 7.82e+02 7.70e-01  -3.8 1.21e+01  -4.1 1.00e+00 6.25e-02h  5
      23  3.7299464e-02 7.73e+02 7.13e-01  -3.8 1.01e+01  -3.6 1.00e+00 3.12e-02h  6
      24  3.7353476e-02 7.73e+02 7.13e-01  -3.8 1.00e+03  -4.1 1.33e-01 6.93e-05h 12
      25  3.7540946e-02 7.70e+02 6.94e-01  -3.8 3.20e+01  -3.7 1.00e+00 7.81e-03h  8
      26  3.7635287e-02 7.68e+02 6.70e-01  -3.8 1.66e+01  -3.3 1.00e+00 7.81e-03h  8
      27  3.7691373e-02 7.66e+02 6.59e-01  -3.8 1.02e+01  -2.8 1.00e+00 7.81e-03h  8
      28  3.7712405e-02 7.65e+02 6.56e-01  -3.8 7.93e+00  -2.4 1.00e+00 3.91e-03h  9
      29  3.7727422e-02 7.65e+02 6.56e-01  -3.8 1.99e+02  -2.9 2.91e-01 9.67e-05h 12
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      30  3.7750003e-02 7.65e+02 6.55e-01  -3.8 2.40e+01  -2.4 1.00e+00 1.26e-03h 10
      31  3.7759625e-02 7.65e+02 6.54e-01  -3.8 1.10e+01  -2.0 6.16e-01 1.25e-03h 10
      32  3.9047737e-02 5.77e+04 4.83e+07  -3.8 9.10e+00  -1.6 1.00e+00 5.66e-01w  1
      33  2.2123190e+02 5.40e+04 1.21e+09  -3.8 3.69e+05  -2.1 5.39e-05 4.42e-02w  1
      34  3.1400350e+02 5.30e+04 1.15e+09  -3.8 3.42e+05  -2.6 4.66e-02 1.67e-02w  1
      35  3.7762112e-02 7.65e+02 6.53e-01  -3.8 5.20e+05  -3.0 1.00e+00 1.11e-03h  9
      36  3.7765169e-02 7.65e+02 6.53e-01  -3.8 1.43e+01  -1.7 3.39e-01 3.50e-04h 11
      37  3.7764337e-02 7.64e+02 6.53e-01  -3.8 1.33e+01  -1.3 1.00e+00 3.76e-04h 11
      38  3.7766965e-02 7.64e+02 6.53e-01  -3.8 4.91e+01  -1.8 1.41e-01 7.24e-05h 12
      39  3.7765931e-02 7.64e+02 6.52e-01  -3.8 1.84e+01  -1.3 1.00e+00 2.72e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      40  3.7754276e-02 7.65e+02 7.94e-01  -3.8 2.09e+01  -0.9 3.18e-01 6.53e-04h 10
      41  7.6303864e-03 1.06e+02 1.05e+00  -3.8 5.76e+01    -  1.00e+00 1.00e+00h  1
      42  7.5438847e-03 2.72e+02 3.00e+01  -3.8 5.00e-01  -1.4 9.99e-01 1.00e+00h  1
      43  7.5439291e-03 4.26e+01 1.15e+01  -3.8 7.45e-02  -1.9 1.00e+00 1.00e+00h  1
      44  7.3820017e-03 7.35e-03 2.57e-03  -3.8 3.53e+00    -  1.00e+00 1.00e+00h  1
      45  7.3821808e-03 2.26e-04 4.85e-06  -3.8 1.04e-03  -2.3 1.00e+00 1.00e+00h  1
      46  5.3334940e-03 6.09e+03 6.40e-02  -5.7 1.81e+01    -  8.56e-01 1.00e+00h  1
      47  6.3493242e-03 1.38e+03 6.96e-03  -5.7 3.81e+00    -  9.36e-01 1.00e+00h  1
      48  5.2784720e-03 1.30e+02 1.60e-03  -5.7 5.62e+00    -  1.00e+00 1.00e+00h  1
      49  5.2175165e-03 6.53e-01 1.09e-05  -5.7 3.60e-01    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      50  5.2173550e-03 1.61e-03 7.67e-09  -5.7 2.35e-03    -  1.00e+00 1.00e+00h  1
      51  5.2173545e-03 2.04e-05 1.02e-09  -8.6 6.55e-03    -  1.00e+00 1.00e+00h  1
      52  5.2173545e-03 7.45e-09 2.51e-14  -8.6 4.74e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 52
    
                                       (scaled)                 (unscaled)
    Objective...............:   5.2173545122733110e-03    5.2173545122733110e-03
    Dual infeasibility......:   2.5085432767634323e-14    2.5085432767634323e-14
    Constraint violation....:   1.4104644499482428e-11    7.4505805969238281e-09
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    7.4505805969238281e-09
    
    
    Number of objective function evaluations             = 323
    Number of objective gradient evaluations             = 53
    Number of equality constraint evaluations            = 323
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 53
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 52
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.223
    Total CPU secs in NLP function evaluations           =      0.541
    
    EXIT: Optimal Solution Found.
    Ipopt 3.12.13: max_iter=6000
    
    
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
    
    Number of nonzeros in equality constraint Jacobian...:    10950
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:     6000
    
    Total number of variables............................:     2952
                         variables with only lower bounds:      150
                    variables with lower and upper bounds:      600
                         variables with only upper bounds:        0
    Total number of equality constraints.................:     2950
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  4.3600768e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  7.0635325e-02 1.79e+03 1.37e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  5.2648042e-03 4.47e+02 1.23e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  4.9626977e-03 5.98e+01 1.41e-03  -1.0 1.99e+00    -  9.95e-01 1.00e+00h  1
       4  4.8622210e-03 8.67e-03 2.30e-04  -1.7 8.97e-01    -  1.00e+00 1.00e+00h  1
       5  4.8483397e-03 6.56e-01 2.11e-04  -3.8 1.70e-01    -  1.00e+00 1.00e+00h  1
       6  3.9036317e-03 1.00e+04 6.33e-02  -3.8 2.83e+01    -  1.00e+00 5.00e-01h  2
       7  4.0163765e-03 9.42e+03 5.94e-02  -3.8 3.29e+01    -  1.00e+00 6.25e-02h  5
       8  4.1217865e-03 8.84e+03 5.58e-02  -3.8 2.93e+01    -  1.00e+00 6.25e-02h  5
       9  4.1232600e-03 8.83e+03 5.58e-02  -3.8 2.60e+01    -  1.00e+00 9.77e-04h 11
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.1239957e-03 8.83e+03 5.57e-02  -3.8 2.60e+01    -  1.00e+00 4.88e-04h 12
      11  6.6351510e-03 2.59e+03 2.84e-02  -3.8 2.60e+01    -  1.00e+00 1.00e+00s 22
      12r 6.6351510e-03 2.59e+03 9.99e+02   1.0 0.00e+00  -4.0 0.00e+00 0.00e+00R  1
      13r 6.4224721e-03 2.29e+03 9.97e+02   1.0 4.16e+04    -  4.67e-03 9.91e-04f  1
      14r 6.1694827e-03 1.83e+03 9.93e+02   1.0 3.83e+03    -  3.41e-03 2.58e-03f  1
      15r 6.1322543e-03 9.21e+02 9.45e+02   1.0 2.74e+03    -  5.57e-02 6.91e-03f  1
      16r 1.0604780e-02 1.83e+02 8.63e+02   1.0 4.25e+02    -  4.48e-01 5.22e-02f  1
      17r 3.3435672e-02 1.45e+01 1.74e+02   1.0 7.30e+01    -  1.00e+00 5.10e-01f  1
      18r 8.7410423e-03 3.10e+01 3.91e+01   0.3 4.09e+01    -  1.00e+00 7.97e-01f  1
      19r 5.0256480e-03 6.86e+00 5.96e-01  -0.4 3.50e+01    -  1.00e+00 9.91e-01f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      20r 4.8168005e-03 6.00e-01 7.83e+00  -2.7 3.62e+01    -  1.00e+00 9.14e-01f  1
      21r 4.2778578e-03 1.81e-03 1.07e-02  -2.7 1.77e+01    -  1.00e+00 1.00e+00f  1
      22  3.8343705e-03 3.09e+01 6.90e-03  -3.8 9.45e+00    -  1.00e+00 1.00e+00h  1
      23  3.8325110e-03 1.06e+00 1.05e-05  -3.8 2.22e-01    -  1.00e+00 1.00e+00h  1
      24  3.8324434e-03 6.14e-04 4.45e-09  -3.8 5.98e-03    -  1.00e+00 1.00e+00h  1
      25  3.8304547e-03 4.37e-01 1.42e-05  -5.7 7.26e-01    -  1.00e+00 1.00e+00h  1
      26  3.8304226e-03 3.32e-04 3.71e-09  -5.7 8.28e-03    -  1.00e+00 1.00e+00h  1
      27  3.8304222e-03 6.01e-05 2.03e-09  -8.6 8.80e-03    -  1.00e+00 1.00e+00h  1
      28  3.8304222e-03 1.49e-08 2.51e-14  -8.6 1.17e-06    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 28
    
                                       (scaled)                 (unscaled)
    Objective...............:   3.8304222329600532e-03    3.8304222329600532e-03
    Dual infeasibility......:   2.5085432767634323e-14    2.5085432767634323e-14
    Constraint violation....:   1.4104644499482428e-11    1.4901161193847656e-08
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    1.4901161193847656e-08
    
    
    Number of objective function evaluations             = 82
    Number of objective gradient evaluations             = 21
    Number of equality constraint evaluations            = 82
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 30
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 28
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.083
    Total CPU secs in NLP function evaluations           =      0.157
    
    EXIT: Optimal Solution Found.
    Ipopt 3.12.13: max_iter=6000
    
    
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
    
    Number of nonzeros in equality constraint Jacobian...:    10950
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:     6000
    
    Total number of variables............................:     2952
                         variables with only lower bounds:      150
                    variables with lower and upper bounds:      600
                         variables with only upper bounds:        0
    Total number of equality constraints.................:     2950
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  6.0208549e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  1.0838816e-01 1.79e+03 1.48e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  7.3390871e-03 2.00e+03 2.57e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  5.7286917e-03 1.56e+03 2.82e-02  -1.0 1.01e+01    -  9.95e-01 1.00e+00h  1
       4  5.5549160e-03 9.65e+01 1.15e-03  -1.0 1.32e+00    -  1.00e+00 1.00e+00h  1
       5  5.7605292e-03 1.76e-02 3.72e-04  -1.7 5.30e+00    -  1.00e+00 1.00e+00h  1
       6  5.7307692e-03 7.51e-01 2.53e-04  -3.8 1.10e+00    -  1.00e+00 1.00e+00h  1
       7  4.7173790e-03 6.21e+03 4.36e-02  -3.8 1.11e+01    -  1.00e+00 1.00e+00h  1
       8  5.1039514e-03 4.77e+03 3.38e-02  -3.8 1.76e+01    -  1.00e+00 2.50e-01h  3
       9  5.4171659e-03 3.67e+03 2.62e-02  -3.8 1.33e+01    -  1.00e+00 2.50e-01h  3
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  5.4679876e-03 3.56e+03 2.54e-02  -3.8 1.63e+01    -  1.00e+00 3.12e-02h  6
      11  5.4717180e-03 3.55e+03 2.53e-02  -3.8 1.88e+01    -  1.00e+00 1.95e-03h 10
      12  5.4726630e-03 3.55e+03 2.53e-02  -3.8 1.90e+01    -  1.00e+00 4.88e-04h 12
      13  8.6414205e-03 1.72e+03 2.15e-02  -3.8 1.91e+01    -  1.00e+00 1.00e+00h  1
      14  8.5762591e-03 1.64e+03 1.82e-01  -3.8 2.46e+00  -4.0 1.00e+00 2.50e-01h  3
      15  8.7223329e-03 1.48e+03 2.16e-01  -3.8 1.45e+00  -4.5 1.00e+00 5.00e-01h  2
      16  8.8541226e-03 1.39e+03 1.15e-01  -3.8 2.68e+00  -4.1 1.00e+00 2.50e-01h  3
      17  8.9644485e-03 1.36e+03 1.92e-01  -3.8 6.53e+00  -3.6 1.00e+00 6.25e-02h  5
      18  8.9187921e-03 1.37e+03 6.08e+00  -3.8 1.24e+01  -1.4 8.79e-01 1.52e-02h  7
      19  1.6528504e-02 1.17e+02 1.16e+01  -3.8 3.34e+01    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      20  1.6386744e-02 1.48e+02 1.57e+01  -3.8 5.92e-01  -1.9 1.00e+00 1.00e+00h  1
      21  5.6886661e-03 6.76e+00 4.25e+00  -3.8 3.81e+01    -  1.00e+00 1.00e+00h  1
      22  5.6985390e-03 1.03e-03 1.12e-03  -3.8 3.53e-02  -2.3 1.00e+00 1.00e+00h  1
      23  5.6982040e-03 8.06e-06 7.26e-05  -5.7 1.71e-03  -2.8 1.00e+00 1.00e+00h  1
      24  4.4298681e-03 4.65e+02 2.75e-03  -5.7 1.32e+01    -  1.00e+00 1.00e+00f  1
      25  4.4219273e-03 5.08e+01 2.87e-04  -5.7 1.69e-01    -  1.00e+00 1.00e+00h  1
      26  4.4111787e-03 3.37e-01 2.56e-06  -5.7 9.10e-02    -  1.00e+00 1.00e+00h  1
      27  4.4111115e-03 1.05e-07 4.36e-11  -5.7 1.22e-03    -  1.00e+00 1.00e+00h  1
      28  4.4111112e-03 4.56e-05 1.52e-09  -8.6 8.97e-03    -  1.00e+00 1.00e+00h  1
      29  4.4111112e-03 7.45e-09 2.51e-14  -8.6 1.03e-06    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 29
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.4111112080836680e-03    4.4111112080836680e-03
    Dual infeasibility......:   2.5085432767634323e-14    2.5085432767634323e-14
    Constraint violation....:   1.4104644499482428e-11    7.4505805969238281e-09
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    7.4505805969238281e-09
    
    
    Number of objective function evaluations             = 79
    Number of objective gradient evaluations             = 30
    Number of equality constraint evaluations            = 79
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 30
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 29
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.111
    Total CPU secs in NLP function evaluations           =      0.171
    
    EXIT: Optimal Solution Found.
    Ipopt 3.12.13: max_iter=6000
    
    
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
    
    Number of nonzeros in equality constraint Jacobian...:    10950
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:     6000
    
    Total number of variables............................:     2952
                         variables with only lower bounds:      150
                    variables with lower and upper bounds:      600
                         variables with only upper bounds:        0
    Total number of equality constraints.................:     2950
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  5.1309447e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  7.4280191e-02 1.79e+03 1.42e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  6.5366378e-03 4.31e+02 1.23e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  6.1213267e-03 5.66e+01 1.39e-03  -1.0 2.07e+00    -  9.95e-01 1.00e+00h  1
       4  5.9948140e-03 1.83e-03 2.55e-04  -1.7 9.23e-01    -  1.00e+00 1.00e+00h  1
       5  5.9715797e-03 1.09e+00 2.66e-04  -3.8 1.86e-01    -  1.00e+00 1.00e+00h  1
       6  5.9664567e-03 5.00e-02 9.24e-05  -3.8 2.94e-02  -4.0 1.00e+00 1.00e+00h  1
       7  4.6751281e-03 7.49e+03 4.29e-02  -5.7 1.18e+01    -  8.44e-01 1.00e+00h  1
       8  5.5169707e-03 1.70e+03 1.31e-02  -5.7 6.11e+00    -  9.27e-01 1.00e+00h  1
       9  4.4746521e-03 9.59e+01 2.20e-03  -5.7 9.31e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.4130116e-03 2.47e+00 4.93e-05  -5.7 9.78e-01    -  1.00e+00 1.00e+00h  1
      11  4.4127379e-03 1.37e-04 2.83e-08  -5.7 1.00e-01    -  1.00e+00 1.00e+00h  1
      12  4.4127375e-03 3.66e-05 1.51e-09  -8.6 7.68e-03    -  1.00e+00 1.00e+00h  1
      13  4.4127375e-03 1.49e-08 2.51e-14  -8.6 7.67e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.4127374700608252e-03    4.4127374700608252e-03
    Dual infeasibility......:   2.5113188343249952e-14    2.5113188343249952e-14
    Constraint violation....:   1.4104644499482428e-11    1.4901161193847656e-08
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    1.4901161193847656e-08
    
    
    Number of objective function evaluations             = 14
    Number of objective gradient evaluations             = 14
    Number of equality constraint evaluations            = 14
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 14
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 13
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.037
    Total CPU secs in NLP function evaluations           =      0.047
    
    EXIT: Optimal Solution Found.
    Ipopt 3.12.13: max_iter=6000
    
    
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
    
    Number of nonzeros in equality constraint Jacobian...:    10950
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:     6000
    
    Total number of variables............................:     2952
                         variables with only lower bounds:      150
                    variables with lower and upper bounds:      600
                         variables with only upper bounds:        0
    Total number of equality constraints.................:     2950
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  5.7202259e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  8.3796778e-02 1.79e+03 1.40e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  7.4096591e-03 5.10e+02 1.33e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  6.8779909e-03 9.77e+01 1.99e-03  -1.0 2.37e+00    -  9.95e-01 1.00e+00h  1
       4  6.6898021e-03 3.05e-02 4.06e-04  -1.7 1.03e+00    -  1.00e+00 1.00e+00h  1
       5  6.6611005e-03 1.38e+00 3.02e-04  -3.8 1.62e-01    -  1.00e+00 1.00e+00h  1
       6  6.6548288e-03 5.56e-02 1.06e-04  -3.8 3.22e-02  -4.0 1.00e+00 1.00e+00h  1
       7  5.0192491e-03 9.98e+03 5.66e-02  -5.7 1.39e+01    -  8.22e-01 1.00e+00h  1
       8  5.9374745e-03 5.54e+03 3.14e-02  -5.7 8.76e+00    -  8.96e-01 5.00e-01h  2
       9  5.4714275e-03 1.27e+03 1.18e-02  -5.7 4.65e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.7731041e-03 1.77e+00 9.31e-04  -5.7 6.02e+00    -  1.00e+00 1.00e+00h  1
      11  4.7527121e-03 2.97e-02 3.03e-06  -5.7 5.42e-01    -  1.00e+00 1.00e+00h  1
      12  4.7527121e-03 1.96e-05 6.39e-10  -5.7 5.04e-03    -  1.00e+00 1.00e+00h  1
      13  4.7527118e-03 2.87e-05 1.29e-09  -8.6 7.27e-03    -  1.00e+00 1.00e+00h  1
      14  4.7527118e-03 1.49e-08 2.51e-14  -8.6 6.36e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 14
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.7527118447067912e-03    4.7527118447067912e-03
    Dual infeasibility......:   2.5085432767634323e-14    2.5085432767634323e-14
    Constraint violation....:   1.4104644499482428e-11    1.4901161193847656e-08
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    1.4901161193847656e-08
    
    
    Number of objective function evaluations             = 16
    Number of objective gradient evaluations             = 15
    Number of equality constraint evaluations            = 16
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 15
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 14
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.058
    Total CPU secs in NLP function evaluations           =      0.074
    
    EXIT: Optimal Solution Found.
    

.. code:: ipython3

    display(bootstrap_theta)



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
          <th>fs.properties.tau[('benzene', 'toluene')]</th>
          <th>fs.properties.tau[('toluene', 'benzene')]</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>0.508343</td>
          <td>-0.431948</td>
        </tr>
        <tr>
          <th>1</th>
          <td>0.490949</td>
          <td>-0.418673</td>
        </tr>
        <tr>
          <th>2</th>
          <td>0.464479</td>
          <td>-0.398766</td>
        </tr>
        <tr>
          <th>3</th>
          <td>0.480199</td>
          <td>-0.410741</td>
        </tr>
        <tr>
          <th>4</th>
          <td>0.504068</td>
          <td>-0.428914</td>
        </tr>
        <tr>
          <th>5</th>
          <td>0.506408</td>
          <td>-0.430699</td>
        </tr>
        <tr>
          <th>6</th>
          <td>0.430020</td>
          <td>-0.372046</td>
        </tr>
        <tr>
          <th>7</th>
          <td>0.455372</td>
          <td>-0.392116</td>
        </tr>
        <tr>
          <th>8</th>
          <td>0.465412</td>
          <td>-0.399478</td>
        </tr>
        <tr>
          <th>9</th>
          <td>0.483388</td>
          <td>-0.413251</td>
        </tr>
      </tbody>
    </table>
    </div>

