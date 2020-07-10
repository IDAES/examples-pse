Parameter Estimation Using Flash Unit Model
-------------------------------------------

In this module, we will be using Pyomo's ``parmest`` tool in conjuction
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

-  NRTL Model - https://idaes-pse.readthedocs.io/en/stable/
-  parmest -
   https://pyomo.readthedocs.io/en/stable/contributed\_packages/parmest/index.html

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


Parameter estimation using parmest
----------------------------------

In addition to providing a method to return an initialized model, the
``parmest`` tool needs the following:

-  List of variable names to be estimated
-  Dataset with multiple scenarios
-  Expression to compute the sum of squared errors

In this example, we only estimate the binary interaction parameter
(``tau_ij``). Given that this variable is usually indexed as
``tau_ij = Var(component_list, component_list)``, there are 2\*2=4
degrees of freedom. However, when i=j, the binary interaction parameter
is 0. Therefore, in this problem, we estimate the binary interaction
parameter for the following variables only:

-  fs.properties.tau['benzene', 'toluene']
-  fs.properties.tau['toluene', 'benzene']

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


Pyomo's ``parmest`` tool supports the following data formats: - pandas
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
      14  4.6633488e-03 1.49e-08 2.51e-14  -8.6 6.81e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 14
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.6633488370395396e-03    4.6633488370395396e-03
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.050
    Total CPU secs in NLP function evaluations           =      0.052
    
    EXIT: Optimal Solution Found.
    

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
    fs.properties.tau[('benzene', 'toluene')] = 0.47810867841011423
    fs.properties.tau[('toluene', 'benzene')] = -0.40924465377594205


Using the data that was provided, we have estimated the binary
interaction parameters in the NRTL model for a benzene-toluene mixture.
Although the dataset that was provided was temperature dependent, in
this example we have estimated a single value that fits best for all
temperatures.

Advanced options for parmest: bootstrapping
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pyomo's ``parmest`` tool allows for bootstrapping where the parameter
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
       0  5.4707694e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  7.3519732e-02 1.79e+03 1.33e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  6.8692631e-03 4.98e+02 1.27e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  6.4294329e-03 8.78e+01 1.78e-03  -1.0 2.15e+00    -  9.95e-01 1.00e+00h  1
       4  6.2636513e-03 2.21e-02 3.22e-04  -1.7 1.33e+00    -  1.00e+00 1.00e+00h  1
       5  6.2398781e-03 1.13e+00 2.68e-04  -3.8 2.83e-01    -  1.00e+00 1.00e+00h  1
       6  6.2344454e-03 5.17e-02 9.59e-05  -3.8 3.02e-02  -4.0 1.00e+00 1.00e+00h  1
       7  4.8479214e-03 8.40e+03 4.75e-02  -5.7 1.24e+01    -  8.35e-01 1.00e+00h  1
       8  6.1204569e-03 1.93e+03 1.48e-02  -5.7 6.86e+00    -  9.15e-01 1.00e+00h  1
       9  4.6688343e-03 1.55e+02 2.94e-03  -5.7 1.05e+01    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.5791378e-03 5.95e+00 8.58e-05  -5.7 9.12e-01    -  1.00e+00 1.00e+00h  1
      11  4.5784214e-03 7.07e-03 3.78e-07  -5.7 1.73e-01    -  1.00e+00 1.00e+00h  1
      12  4.5784192e-03 3.82e-08 1.84e-11  -5.7 3.30e-04    -  1.00e+00 1.00e+00h  1
      13  4.5784190e-03 3.40e-05 1.40e-09  -8.6 7.59e-03    -  1.00e+00 1.00e+00h  1
      14  4.5784190e-03 1.49e-08 2.51e-14  -8.6 7.34e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 14
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.5784189532737853e-03    4.5784189532737853e-03
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.059
    Total CPU secs in NLP function evaluations           =      0.062
    
    EXIT: Optimal Solution Found.
    Ipopt 3.13.2: 
    
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
       0  5.7325316e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  8.8966801e-02 1.79e+03 1.45e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  7.8164291e-03 7.45e+02 1.64e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  7.0330176e-03 2.51e+02 4.03e-03  -1.0 3.24e+00    -  9.95e-01 1.00e+00h  1
       4  6.7077678e-03 3.20e-01 1.03e-03  -1.7 1.90e+00    -  1.00e+00 1.00e+00h  1
       5  6.6817265e-03 1.50e+00 2.90e-04  -3.8 2.16e-01    -  1.00e+00 1.00e+00h  1
       6  6.6753354e-03 4.89e-02 1.00e-04  -3.8 3.11e-02  -4.0 1.00e+00 1.00e+00h  1
       7  4.9745524e-03 1.03e+04 5.35e-02  -5.7 1.41e+01    -  8.19e-01 1.00e+00h  1
       8  5.9256561e-03 5.70e+03 2.96e-02  -5.7 8.76e+00    -  8.91e-01 5.00e-01h  2
       9  5.4841133e-03 1.35e+03 1.17e-02  -5.7 4.16e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.7741679e-03 7.55e+00 9.12e-04  -5.7 5.71e+00    -  1.00e+00 1.00e+00h  1
      11  4.7520303e-03 2.81e-02 2.98e-06  -5.7 4.62e-01    -  1.00e+00 1.00e+00h  1
      12  4.7520316e-03 7.57e-06 4.00e-10  -5.7 4.65e-03    -  1.00e+00 1.00e+00h  1
      13  4.7520313e-03 2.79e-05 1.29e-09  -8.6 6.90e-03    -  1.00e+00 1.00e+00h  1
      14  4.7520313e-03 1.12e-08 2.51e-14  -8.6 5.93e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 14
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.7520313407266564e-03    4.7520313407266564e-03
    Dual infeasibility......:   2.5085432767634323e-14    2.5085432767634323e-14
    Constraint violation....:   1.4104644499482428e-11    1.1175870895385742e-08
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    1.1175870895385742e-08
    
    
    Number of objective function evaluations             = 16
    Number of objective gradient evaluations             = 15
    Number of equality constraint evaluations            = 16
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 15
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 14
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.055
    Total CPU secs in NLP function evaluations           =      0.061
    
    EXIT: Optimal Solution Found.
    Ipopt 3.13.2: 
    
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
       0  5.9995772e+01 5.41e+05 4.00e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  8.2802289e-02 1.66e+03 1.45e+00  -1.0 1.32e+04    -  9.82e-01 1.00e+00h  1
       2  8.2317854e-03 4.45e+02 1.26e-01  -1.0 4.39e+02    -  9.90e-01 1.00e+00h  1
       3  7.5746608e-03 6.83e+01 1.59e-03  -1.0 2.53e+00    -  9.96e-01 1.00e+00h  1
       4  7.2952130e-03 1.02e-03 4.06e-04  -1.7 2.51e+00    -  1.00e+00 1.00e+00h  1
       5  7.2510625e-03 1.87e+00 3.24e-04  -3.8 6.12e-01    -  1.00e+00 1.00e+00h  1
       6  7.2434048e-03 6.32e-02 1.24e-04  -3.8 3.33e-02  -4.0 1.00e+00 1.00e+00h  1
       7  5.2123829e-03 1.20e+04 6.18e-02  -5.7 1.51e+01    -  8.06e-01 1.00e+00h  1
       8  6.6050838e-03 6.65e+03 3.41e-02  -5.7 1.17e+01    -  8.73e-01 5.00e-01h  2
       9  6.8609204e-03 5.82e+03 2.94e-02  -5.7 3.06e+02    -  1.00e+00 1.25e-01h  4
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  5.9549584e-03 1.46e+03 2.80e-02  -5.7 3.02e+01    -  1.00e+00 1.00e+00h  1
      11  4.9807954e-03 5.67e+00 1.97e-03  -5.7 1.09e+01    -  1.00e+00 1.00e+00h  1
      12  4.9466179e-03 2.83e-02 4.68e-06  -5.7 8.15e-01    -  1.00e+00 1.00e+00h  1
      13  4.9466268e-03 3.09e-05 1.23e-09  -5.7 6.83e-03    -  1.00e+00 1.00e+00h  1
      14  4.9466266e-03 2.28e-05 1.14e-09  -8.6 6.44e-03    -  1.00e+00 1.00e+00h  1
      15  4.9466265e-03 1.12e-08 2.51e-14  -8.6 4.94e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 15
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.9466265479011706e-03    4.9466265479011706e-03
    Dual infeasibility......:   2.5140943918865581e-14    2.5140943918865581e-14
    Constraint violation....:   1.4104644499482428e-11    1.1175870895385742e-08
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    1.1175870895385742e-08
    
    
    Number of objective function evaluations             = 20
    Number of objective gradient evaluations             = 16
    Number of equality constraint evaluations            = 20
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 16
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 15
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.054
    Total CPU secs in NLP function evaluations           =      0.065
    
    EXIT: Optimal Solution Found.
    Ipopt 3.13.2: 
    
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
       0  5.5655824e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  8.5677817e-02 1.79e+03 1.39e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  7.3017904e-03 6.72e+02 1.52e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  6.6769461e-03 1.94e+02 3.24e-03  -1.0 2.82e+00    -  9.95e-01 1.00e+00h  1
       4  6.4171889e-03 1.88e-01 7.67e-04  -1.7 1.26e+00    -  1.00e+00 1.00e+00h  1
       5  6.3934774e-03 1.27e+00 2.91e-04  -3.8 1.38e-01    -  1.00e+00 1.00e+00h  1
       6  6.3876910e-03 4.90e-02 9.65e-05  -3.8 3.11e-02  -4.0 1.00e+00 1.00e+00h  1
       7  4.8740303e-03 9.19e+03 5.24e-02  -5.7 1.35e+01    -  8.28e-01 1.00e+00h  1
       8  5.5902708e-03 5.11e+03 2.91e-02  -5.7 7.50e+00    -  9.04e-01 5.00e-01h  2
       9  5.2016521e-03 1.15e+03 1.06e-02  -5.7 4.51e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.6510048e-03 1.73e+00 7.50e-04  -5.7 5.49e+00    -  1.00e+00 1.00e+00h  1
      11  4.6361885e-03 1.95e-02 2.05e-06  -5.7 4.62e-01    -  1.00e+00 1.00e+00h  1
      12  4.6361891e-03 9.60e-06 3.08e-10  -5.7 3.51e-03    -  1.00e+00 1.00e+00h  1
      13  4.6361889e-03 3.14e-05 1.36e-09  -8.6 7.46e-03    -  1.00e+00 1.00e+00h  1
      14  4.6361889e-03 1.49e-08 2.51e-14  -8.6 6.85e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 14
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.6361888910869207e-03    4.6361888910869207e-03
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.052
    Total CPU secs in NLP function evaluations           =      0.071
    
    EXIT: Optimal Solution Found.
    Ipopt 3.13.2: 
    
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
       0  5.8105601e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  9.4072625e-02 1.79e+03 1.30e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  8.0897952e-03 1.54e+03 2.21e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  6.5495440e-03 1.01e+03 1.35e-02  -1.0 6.47e+00    -  9.95e-01 1.00e+00h  1
       4  6.2100676e-03 5.12e+01 5.52e-04  -1.0 1.79e+00    -  1.00e+00 1.00e+00h  1
       5  6.4946297e-03 1.44e-02 2.37e-04  -2.5 2.02e+00    -  1.00e+00 1.00e+00h  1
       6  6.3305112e-03 5.16e+01 5.87e-04  -3.8 1.02e+00    -  1.00e+00 1.00e+00h  1
       7  6.3237071e-03 1.33e-01 2.06e-04  -3.8 3.61e-02  -4.0 1.00e+00 1.00e+00h  1
       8  5.1808541e-03 2.88e+03 2.06e-02  -5.7 8.52e+00    -  9.00e-01 1.00e+00h  1
       9  4.8657317e-03 4.58e+02 6.17e-03  -5.7 7.59e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.6875700e-03 3.61e-01 1.34e-04  -5.7 2.90e+00    -  1.00e+00 1.00e+00h  1
      11  4.6852088e-03 4.86e-02 1.48e-06  -5.7 2.39e-01    -  1.00e+00 1.00e+00h  1
      12  4.6852027e-03 8.92e-06 1.89e-10  -5.7 2.31e-03    -  1.00e+00 1.00e+00h  1
      13  4.6852025e-03 3.17e-05 1.33e-09  -8.6 7.13e-03    -  1.00e+00 1.00e+00h  1
      14  4.6852025e-03 7.45e-09 2.51e-14  -8.6 6.63e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 14
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.6852024567361631e-03    4.6852024567361631e-03
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.053
    Total CPU secs in NLP function evaluations           =      0.066
    
    EXIT: Optimal Solution Found.
    Ipopt 3.13.2: 
    
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
       0  6.3810367e+01 5.41e+05 4.00e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  1.1333874e-01 1.66e+03 1.52e+00  -1.0 1.32e+04    -  9.82e-01 1.00e+00h  1
       2  8.4347241e-03 3.58e+03 3.61e-01  -1.0 4.39e+02    -  9.90e-01 1.00e+00h  1
       3  6.1759862e-03 3.86e+03 1.12e-01  -1.0 3.63e+01    -  9.96e-01 1.00e+00h  1
       4  5.0634150e-03 4.30e+02 5.11e-03  -1.0 6.32e+00    -  1.00e+00 1.00e+00h  1
       5  6.3891147e-03 1.03e+01 1.68e-03  -1.0 2.13e+01    -  1.00e+00 1.00e+00h  1
       6  5.9591355e-03 1.69e-02 3.48e-04  -2.5 5.48e+00    -  1.00e+00 1.00e+00h  1
       7  5.8495788e-03 3.42e+01 5.12e-04  -3.8 1.03e+00    -  1.00e+00 1.00e+00h  1
       8  5.8444434e-03 2.32e-02 1.03e-04  -3.8 3.60e-02  -4.0 1.00e+00 1.00e+00h  1
       9  4.8991017e-03 3.43e+03 2.66e-02  -5.7 9.50e+00    -  8.91e-01 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.8132355e-03 6.42e+02 7.63e-03  -5.7 5.99e+00    -  9.91e-01 1.00e+00h  1
      11  4.5444624e-03 7.94e-01 3.29e-04  -5.7 4.58e+00    -  1.00e+00 1.00e+00h  1
      12  4.5393218e-03 4.66e-02 1.91e-06  -5.7 4.08e-01    -  1.00e+00 1.00e+00h  1
      13  4.5393172e-03 1.92e-05 4.14e-10  -5.7 4.43e-03    -  1.00e+00 1.00e+00h  1
      14  4.5393169e-03 4.22e-05 1.36e-09  -8.6 8.94e-03    -  1.00e+00 1.00e+00h  1
      15  4.5393169e-03 1.12e-08 2.51e-14  -8.6 9.82e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 15
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.5393169202680153e-03    4.5393169202680153e-03
    Dual infeasibility......:   2.5085432767634323e-14    2.5085432767634323e-14
    Constraint violation....:   1.4104644499482428e-11    1.1175870895385742e-08
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    1.1175870895385742e-08
    
    
    Number of objective function evaluations             = 16
    Number of objective gradient evaluations             = 16
    Number of equality constraint evaluations            = 16
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 16
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 15
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.057
    Total CPU secs in NLP function evaluations           =      0.065
    
    EXIT: Optimal Solution Found.
    Ipopt 3.13.2: 
    
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
       0  5.4257179e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  7.1794849e-02 1.79e+03 1.49e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  7.1002290e-03 3.37e+02 1.13e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  6.6779684e-03 2.35e+01 8.27e-04  -1.0 1.90e+00    -  9.95e-01 1.00e+00h  1
       4  6.5489646e-03 3.27e-02 2.53e-04  -2.5 1.86e+00    -  1.00e+00 1.00e+00h  1
       5  6.3471407e-03 6.90e+01 7.36e-04  -3.8 1.10e+00    -  1.00e+00 1.00e+00h  1
       6  6.3397608e-03 7.04e-01 4.64e-04  -3.8 3.80e-02  -4.0 1.00e+00 1.00e+00h  1
       7  5.2156454e-03 2.00e+03 1.88e-02  -3.8 7.96e+00    -  1.00e+00 1.00e+00h  1
       8  4.7402891e-03 2.30e+02 4.54e-03  -3.8 8.73e+00    -  1.00e+00 1.00e+00h  1
       9  4.6477092e-03 5.01e-01 5.47e-05  -3.8 1.70e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.6475861e-03 1.82e-02 3.14e-07  -3.8 8.49e-02    -  1.00e+00 1.00e+00h  1
      11  4.6460490e-03 2.06e-01 9.19e-06  -5.7 5.65e-01    -  1.00e+00 1.00e+00h  1
      12  4.6460244e-03 6.38e-05 9.24e-10  -5.7 4.22e-03    -  1.00e+00 1.00e+00h  1
      13  4.6460241e-03 2.95e-05 1.36e-09  -8.6 6.91e-03    -  1.00e+00 1.00e+00h  1
      14  4.6460241e-03 1.49e-08 2.51e-14  -8.6 6.17e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 14
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.6460241092730701e-03    4.6460241092730701e-03
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
    Total CPU secs in NLP function evaluations           =      0.050
    
    EXIT: Optimal Solution Found.
    Ipopt 3.13.2: 
    
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
       0  6.3446419e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  9.6896457e-02 1.79e+03 1.35e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  8.8416342e-03 1.29e+03 2.03e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  7.4631062e-03 7.70e+02 1.11e-02  -1.0 5.28e+00    -  9.95e-01 1.00e+00h  1
       4  7.0430673e-03 3.86e+01 5.06e-04  -1.0 2.07e+00    -  1.00e+00 1.00e+00h  1
       5  7.2619465e-03 3.37e-02 2.61e-04  -2.5 1.32e+00    -  1.00e+00 1.00e+00h  1
       6  7.0484847e-03 6.82e+01 8.30e-04  -3.8 1.23e+00    -  1.00e+00 1.00e+00h  1
       7  7.0398190e-03 4.50e-01 4.60e-04  -3.8 4.40e-02  -4.0 1.00e+00 1.00e+00h  1
       8  5.7572393e-03 2.32e+03 2.09e-02  -3.8 8.69e+00    -  1.00e+00 1.00e+00h  1
       9  5.1531534e-03 2.83e+02 5.38e-03  -3.8 1.06e+01    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  5.0236976e-03 4.59e-01 7.56e-05  -3.8 2.20e+00    -  1.00e+00 1.00e+00h  1
      11  5.0233365e-03 2.67e-02 4.75e-07  -3.8 1.18e-01    -  1.00e+00 1.00e+00h  1
      12  5.0219051e-03 1.68e-01 7.59e-06  -5.7 5.84e-01    -  1.00e+00 1.00e+00h  1
      13  5.0218823e-03 4.16e-05 6.14e-10  -5.7 3.87e-03    -  1.00e+00 1.00e+00h  1
      14  5.0218821e-03 2.44e-05 1.13e-09  -8.6 7.15e-03    -  1.00e+00 1.00e+00h  1
      15  5.0218821e-03 7.45e-09 2.51e-14  -8.6 5.70e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 15
    
                                       (scaled)                 (unscaled)
    Objective...............:   5.0218820705220570e-03    5.0218820705220570e-03
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.056
    Total CPU secs in NLP function evaluations           =      0.058
    
    EXIT: Optimal Solution Found.
    Ipopt 3.13.2: 
    
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
       0  4.3353769e+01 4.99e+05 3.68e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  4.6880789e-02 1.41e+03 1.20e+00  -1.0 1.22e+04    -  9.83e-01 1.00e+00h  1
       2  5.5691074e-03 1.69e+02 6.86e-02  -1.0 3.73e+02    -  9.90e-01 1.00e+00h  1
       3  5.4031906e-03 1.34e+00 7.80e-05  -1.0 9.16e-01    -  9.96e-01 1.00e+00h  1
       4  5.3431431e-03 3.95e-02 2.33e-04  -2.5 2.13e+00    -  1.00e+00 1.00e+00h  1
       5  5.2190274e-03 4.03e+01 4.93e-04  -3.8 8.08e-01    -  1.00e+00 1.00e+00h  1
       6  5.2141809e-03 4.52e-02 1.14e-04  -3.8 3.14e-02  -4.0 1.00e+00 1.00e+00h  1
       7  4.3271825e-03 2.84e+03 1.82e-02  -5.7 8.51e+00    -  9.01e-01 1.00e+00h  1
       8  4.2159497e-03 4.95e+02 5.12e-03  -5.7 5.19e+00    -  1.00e+00 1.00e+00h  1
       9  4.0409050e-03 3.51e-02 1.63e-04  -5.7 3.10e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.0382819e-03 3.96e-02 1.40e-06  -5.7 2.59e-01    -  1.00e+00 1.00e+00h  1
      11  4.0382778e-03 1.01e-05 2.09e-10  -5.7 2.54e-03    -  1.00e+00 1.00e+00h  1
      12  4.0382775e-03 4.66e-05 1.66e-09  -8.6 7.90e-03    -  1.00e+00 1.00e+00h  1
      13  4.0382775e-03 1.49e-08 2.51e-14  -8.6 9.21e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.0382775416164711e-03    4.0382775416164711e-03
    Dual infeasibility......:   2.5140943918865581e-14    2.5140943918865581e-14
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.050
    Total CPU secs in NLP function evaluations           =      0.053
    
    EXIT: Optimal Solution Found.
    Ipopt 3.13.2: 
    
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
       0  5.7110358e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  8.6903910e-02 1.79e+03 1.39e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  7.5805999e-03 7.31e+02 1.58e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  6.8869724e-03 2.37e+02 3.75e-03  -1.0 3.02e+00    -  9.95e-01 1.00e+00h  1
       4  6.5938952e-03 2.84e-01 9.28e-04  -1.7 1.36e+00    -  1.00e+00 1.00e+00h  1
       5  6.5701114e-03 1.37e+00 2.71e-04  -3.8 1.41e-01    -  1.00e+00 1.00e+00h  1
       6  6.5640211e-03 4.89e-02 9.55e-05  -3.8 3.03e-02  -4.0 1.00e+00 1.00e+00h  1
       7  4.9528977e-03 9.90e+03 5.08e-02  -5.7 1.36e+01    -  8.22e-01 1.00e+00h  1
       8  5.8448294e-03 5.50e+03 2.83e-02  -5.7 8.24e+00    -  8.95e-01 5.00e-01h  2
       9  5.4104744e-03 1.28e+03 1.09e-02  -5.7 4.05e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.7420081e-03 4.74e+00 7.79e-04  -5.7 5.31e+00    -  1.00e+00 1.00e+00h  1
      11  4.7222975e-03 3.50e-04 6.17e-07  -5.7 4.68e-01    -  1.00e+00 1.00e+00h  1
      12  4.7223020e-03 3.57e-05 1.56e-09  -8.6 7.68e-03    -  1.00e+00 1.00e+00h  1
      13  4.7223020e-03 1.49e-08 2.51e-14  -8.6 8.46e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.7223020139861649e-03    4.7223020139861649e-03
    Dual infeasibility......:   2.5140943918865581e-14    2.5140943918865581e-14
    Constraint violation....:   1.4104644499482428e-11    1.4901161193847656e-08
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    1.4901161193847656e-08
    
    
    Number of objective function evaluations             = 15
    Number of objective gradient evaluations             = 14
    Number of equality constraint evaluations            = 15
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 14
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 13
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.060
    Total CPU secs in NLP function evaluations           =      0.068
    
    EXIT: Optimal Solution Found.
    

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
          <td>0.472341</td>
          <td>-0.404884</td>
        </tr>
        <tr>
          <th>1</th>
          <td>0.484398</td>
          <td>-0.413968</td>
        </tr>
        <tr>
          <th>2</th>
          <td>0.497148</td>
          <td>-0.423533</td>
        </tr>
        <tr>
          <th>3</th>
          <td>0.476884</td>
          <td>-0.408297</td>
        </tr>
        <tr>
          <th>4</th>
          <td>0.477237</td>
          <td>-0.408654</td>
        </tr>
        <tr>
          <th>5</th>
          <td>0.461687</td>
          <td>-0.397034</td>
        </tr>
        <tr>
          <th>6</th>
          <td>0.480036</td>
          <td>-0.410605</td>
        </tr>
        <tr>
          <th>7</th>
          <td>0.495566</td>
          <td>-0.422568</td>
        </tr>
        <tr>
          <th>8</th>
          <td>0.447394</td>
          <td>-0.385450</td>
        </tr>
        <tr>
          <th>9</th>
          <td>0.481446</td>
          <td>-0.411780</td>
        </tr>
      </tbody>
    </table>
    </div>

