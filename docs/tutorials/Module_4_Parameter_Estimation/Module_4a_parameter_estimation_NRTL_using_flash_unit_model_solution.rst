Module 4a: Parameter Estimation Using Flash Unit Model
------------------------------------------------------

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

-  NRTL Model -
   https://idaes-pse.readthedocs.io/en/latest/model\_libraries/core\_library/property\_models/activity\_coefficient.html
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

    Ipopt 3.13.2: max_iter=6000
    
    
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.043
    Total CPU secs in NLP function evaluations           =      0.046
    
    EXIT: Optimal Solution Found.


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

    Ipopt 3.13.2: max_iter=6000
    
    
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
       0  5.5826379e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  8.1712714e-02 1.79e+03 1.36e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  7.1460920e-03 5.65e+02 1.37e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  6.6259153e-03 1.26e+02 2.33e-03  -1.0 2.42e+00    -  9.95e-01 1.00e+00h  1
       4  6.4210938e-03 6.68e-02 4.87e-04  -1.7 1.20e+00    -  1.00e+00 1.00e+00h  1
       5  6.3961907e-03 1.23e+00 2.90e-04  -3.8 2.13e-01    -  1.00e+00 1.00e+00h  1
       6  6.3904474e-03 5.17e-02 1.01e-04  -3.8 3.13e-02  -4.0 1.00e+00 1.00e+00h  1
       7  4.9034739e-03 9.07e+03 5.27e-02  -5.7 1.33e+01    -  8.29e-01 1.00e+00h  1
       8  5.6088113e-03 5.04e+03 2.94e-02  -5.7 7.51e+00    -  9.06e-01 5.00e-01h  2
       9  5.2069389e-03 1.12e+03 1.06e-02  -5.7 4.70e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.6594284e-03 8.48e-02 7.20e-04  -5.7 5.45e+00    -  1.00e+00 1.00e+00h  1
      11  4.6451964e-03 7.79e-02 4.34e-06  -5.7 5.10e-01    -  1.00e+00 1.00e+00h  1
      12  4.6451879e-03 5.63e-05 1.52e-09  -5.7 7.47e-03    -  1.00e+00 1.00e+00h  1
      13  4.6451877e-03 3.17e-05 1.35e-09  -8.6 7.56e-03    -  1.00e+00 1.00e+00h  1
      14  4.6451877e-03 1.12e-08 2.51e-14  -8.6 7.01e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 14
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.6451876538814735e-03    4.6451876538814735e-03
    Dual infeasibility......:   2.5140943918865581e-14    2.5140943918865581e-14
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.047
    Total CPU secs in NLP function evaluations           =      0.046
    
    EXIT: Optimal Solution Found.
    Ipopt 3.13.2: max_iter=6000
    
    
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
       0  5.6849119e+01 5.41e+05 4.00e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  7.9048528e-02 1.66e+03 1.33e+00  -1.0 1.32e+04    -  9.82e-01 1.00e+00h  1
       2  7.5634600e-03 5.55e+02 1.33e-01  -1.0 4.39e+02    -  9.90e-01 1.00e+00h  1
       3  6.9528802e-03 1.23e+02 2.24e-03  -1.0 2.59e+00    -  9.96e-01 1.00e+00h  1
       4  6.7255563e-03 5.45e-02 5.17e-04  -1.7 1.43e+00    -  1.00e+00 1.00e+00h  1
       5  6.6956422e-03 1.45e+00 3.03e-04  -3.8 2.37e-01    -  1.00e+00 1.00e+00h  1
       6  6.6892408e-03 5.50e-02 1.05e-04  -3.8 3.21e-02  -4.0 1.00e+00 1.00e+00h  1
       7  5.0119333e-03 1.01e+04 5.62e-02  -5.7 1.40e+01    -  8.21e-01 1.00e+00h  1
       8  5.9323319e-03 5.62e+03 3.12e-02  -5.7 8.63e+00    -  8.94e-01 5.00e-01h  2
       9  5.4774765e-03 1.30e+03 1.19e-02  -5.7 4.56e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.7761000e-03 3.17e+00 9.58e-04  -5.7 6.01e+00    -  1.00e+00 1.00e+00h  1
      11  4.7550531e-03 4.56e-03 1.05e-06  -5.7 5.02e-01    -  1.00e+00 1.00e+00h  1
      12  4.7550571e-03 1.48e-05 6.99e-10  -8.6 5.30e-03    -  1.00e+00 1.00e+00h  1
      13  4.7550571e-03 1.49e-08 2.51e-14  -8.6 1.98e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.7550570758593024e-03    4.7550570758593024e-03
    Dual infeasibility......:   2.5085432767634323e-14    2.5085432767634323e-14
    Constraint violation....:   1.4104644499482428e-11    1.4901161193847656e-08
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    1.4901161193847656e-08
    
    
    Number of objective function evaluations             = 15
    Number of objective gradient evaluations             = 14
    Number of equality constraint evaluations            = 15
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 14
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 13
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.038
    Total CPU secs in NLP function evaluations           =      0.039
    
    EXIT: Optimal Solution Found.
    Ipopt 3.13.2: max_iter=6000
    
    
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
       0  5.2366267e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  7.5573044e-02 1.79e+03 1.49e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  6.8173475e-03 3.69e+02 1.18e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  6.3934364e-03 3.37e+01 1.03e-03  -1.0 2.00e+00    -  9.95e-01 1.00e+00h  1
       4  6.2549488e-03 2.03e-02 2.82e-04  -2.5 1.79e+00    -  1.00e+00 1.00e+00h  1
       5  6.0714088e-03 6.28e+01 7.14e-04  -3.8 1.07e+00    -  1.00e+00 1.00e+00h  1
       6  6.0647110e-03 5.65e-01 4.22e-04  -3.8 3.77e-02  -4.0 1.00e+00 1.00e+00h  1
       7  5.0139033e-03 2.02e+03 1.81e-02  -3.8 8.08e+00    -  1.00e+00 1.00e+00h  1
       8  4.6037737e-03 2.44e+02 4.43e-03  -3.8 8.26e+00    -  1.00e+00 1.00e+00h  1
       9  4.5124439e-03 3.81e-01 5.44e-05  -3.8 1.84e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.5107538e-03 1.08e-01 5.88e-06  -5.7 4.96e-01    -  1.00e+00 1.00e+00h  1
      11  4.5107427e-03 1.72e-05 2.23e-10  -5.7 2.08e-03    -  1.00e+00 1.00e+00h  1
      12  4.5107425e-03 3.27e-05 1.45e-09  -8.6 7.29e-03    -  1.00e+00 1.00e+00h  1
      13  4.5107425e-03 1.12e-08 2.51e-14  -8.6 6.87e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.5107424554100004e-03    4.5107424554100004e-03
    Dual infeasibility......:   2.5085432767634323e-14    2.5085432767634323e-14
    Constraint violation....:   1.4104644499482428e-11    1.1175870895385742e-08
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    1.1175870895385742e-08
    
    
    Number of objective function evaluations             = 14
    Number of objective gradient evaluations             = 14
    Number of equality constraint evaluations            = 14
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 14
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 13
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.042
    Total CPU secs in NLP function evaluations           =      0.042
    
    EXIT: Optimal Solution Found.
    Ipopt 3.13.2: max_iter=6000
    
    
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
       0  5.0175845e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  6.9312151e-02 1.79e+03 1.54e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  6.5356379e-03 2.71e+02 1.03e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  6.2022672e-03 6.86e+00 4.24e-04  -1.0 1.62e+00    -  9.95e-01 1.00e+00h  1
       4  6.0805463e-03 4.38e-02 2.37e-04  -2.5 2.38e+00    -  1.00e+00 1.00e+00h  1
       5  5.9076947e-03 5.77e+01 6.52e-04  -3.8 9.81e-01    -  1.00e+00 1.00e+00h  1
       6  5.9012856e-03 3.18e-01 3.11e-04  -3.8 3.63e-02  -4.0 1.00e+00 1.00e+00h  1
       7  4.8414406e-03 2.38e+03 1.94e-02  -5.7 8.40e+00    -  9.09e-01 1.00e+00h  1
       8  4.5261115e-03 3.39e+02 5.34e-03  -5.7 7.50e+00    -  1.00e+00 1.00e+00h  1
       9  4.4057985e-03 7.30e-02 8.84e-05  -5.7 2.41e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.4046004e-03 3.62e-02 9.23e-07  -5.7 1.66e-01    -  1.00e+00 1.00e+00h  1
      11  4.4045959e-03 3.64e-06 7.15e-11  -5.7 1.30e-03    -  1.00e+00 1.00e+00h  1
      12  4.4045956e-03 3.46e-05 1.52e-09  -8.6 7.27e-03    -  1.00e+00 1.00e+00h  1
      13  4.4045956e-03 1.49e-08 2.51e-14  -8.6 7.08e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.4045956234040895e-03    4.4045956234040895e-03
    Dual infeasibility......:   2.5085432767634323e-14    2.5085432767634323e-14
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.051
    Total CPU secs in NLP function evaluations           =      0.057
    
    EXIT: Optimal Solution Found.
    Ipopt 3.13.2: max_iter=6000
    
    
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
       0  5.1271707e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  6.8173142e-02 1.79e+03 1.67e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  6.9814283e-03 2.59e+02 1.07e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  6.5852295e-03 5.08e+00 3.90e-04  -1.0 1.75e+00    -  9.95e-01 1.00e+00h  1
       4  6.2696025e-03 5.98e-02 3.56e-04  -2.5 4.42e+00    -  1.00e+00 1.00e+00h  1
       5  6.0638498e-03 6.93e+01 7.09e-04  -3.8 1.02e+00    -  1.00e+00 1.00e+00h  1
       6  6.0567790e-03 5.53e-01 4.01e-04  -3.8 3.67e-02  -4.0 1.00e+00 1.00e+00h  1
       7  4.9546992e-03 1.98e+03 1.98e-02  -3.8 7.77e+00    -  1.00e+00 1.00e+00h  1
       8  4.5173097e-03 2.30e+02 4.68e-03  -3.8 8.07e+00    -  1.00e+00 1.00e+00h  1
       9  4.4372668e-03 3.79e-01 5.15e-05  -3.8 1.57e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.4356102e-03 1.13e-01 6.55e-06  -5.7 4.70e-01    -  1.00e+00 1.00e+00h  1
      11  4.4355997e-03 1.92e-05 2.75e-10  -5.7 2.13e-03    -  1.00e+00 1.00e+00h  1
      12  4.4355995e-03 3.11e-05 1.49e-09  -8.6 6.71e-03    -  1.00e+00 1.00e+00h  1
      13  4.4355995e-03 1.12e-08 2.51e-14  -8.6 6.13e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.4355994922417388e-03    4.4355994922417388e-03
    Dual infeasibility......:   2.5140943918865581e-14    2.5140943918865581e-14
    Constraint violation....:   1.4104644499482428e-11    1.1175870895385742e-08
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    1.1175870895385742e-08
    
    
    Number of objective function evaluations             = 14
    Number of objective gradient evaluations             = 14
    Number of equality constraint evaluations            = 14
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 14
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 13
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.041
    Total CPU secs in NLP function evaluations           =      0.042
    
    EXIT: Optimal Solution Found.
    Ipopt 3.13.2: max_iter=6000
    
    
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
       0  4.4736806e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  6.1203628e-02 1.79e+03 1.55e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  5.7965141e-03 2.45e+02 9.92e-02  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  5.5315657e-03 2.70e+00 2.72e-04  -1.0 1.44e+00    -  9.95e-01 1.00e+00h  1
       4  5.4152770e-03 3.92e-02 2.15e-04  -2.5 2.54e+00    -  1.00e+00 1.00e+00h  1
       5  5.2802059e-03 4.44e+01 5.00e-04  -3.8 8.14e-01    -  1.00e+00 1.00e+00h  1
       6  5.2751723e-03 1.00e-01 1.54e-04  -3.8 3.11e-02  -4.0 1.00e+00 1.00e+00h  1
       7  4.3691461e-03 2.59e+03 1.83e-02  -5.7 8.08e+00    -  9.06e-01 1.00e+00h  1
       8  4.1915005e-03 4.20e+02 5.27e-03  -5.7 5.64e+00    -  1.00e+00 1.00e+00h  1
       9  4.0631520e-03 2.87e-01 1.10e-04  -5.7 2.58e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.0614978e-03 3.75e-02 1.21e-06  -5.7 2.02e-01    -  1.00e+00 1.00e+00h  1
      11  4.0614942e-03 6.54e-06 1.39e-10  -5.7 1.83e-03    -  1.00e+00 1.00e+00h  1
      12  4.0614939e-03 4.45e-05 1.80e-09  -8.6 7.49e-03    -  1.00e+00 1.00e+00h  1
      13  4.0614939e-03 7.45e-09 2.51e-14  -8.6 8.44e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.0614939265038825e-03    4.0614939265038825e-03
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.034
    Total CPU secs in NLP function evaluations           =      0.041
    
    EXIT: Optimal Solution Found.
    Ipopt 3.13.2: max_iter=6000
    
    
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
       0  4.8808160e+01 5.41e+05 4.00e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  6.3562345e-02 1.66e+03 1.46e+00  -1.0 1.32e+04    -  9.82e-01 1.00e+00h  1
       2  6.4927910e-03 2.97e+02 1.04e-01  -1.0 4.39e+02    -  9.90e-01 1.00e+00h  1
       3  6.1189351e-03 1.21e+01 5.56e-04  -1.0 1.78e+00    -  9.96e-01 1.00e+00h  1
       4  5.9450442e-03 3.92e-02 2.56e-04  -2.5 2.91e+00    -  1.00e+00 1.00e+00h  1
       5  5.7740014e-03 5.74e+01 6.08e-04  -3.8 9.33e-01    -  1.00e+00 1.00e+00h  1
       6  5.7678293e-03 3.27e-01 2.90e-04  -3.8 3.42e-02  -4.0 1.00e+00 1.00e+00h  1
       7  4.7348312e-03 2.31e+03 1.88e-02  -5.7 8.06e+00    -  9.11e-01 1.00e+00h  1
       8  4.4306928e-03 3.27e+02 5.07e-03  -5.7 7.08e+00    -  1.00e+00 1.00e+00h  1
       9  4.3238797e-03 5.98e-02 7.40e-05  -5.7 2.19e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.3228362e-03 3.24e-02 7.92e-07  -5.7 1.49e-01    -  1.00e+00 1.00e+00h  1
      11  4.3228325e-03 2.90e-06 5.50e-11  -5.7 1.10e-03    -  1.00e+00 1.00e+00h  1
      12  4.3228323e-03 3.59e-05 1.54e-09  -8.6 7.04e-03    -  1.00e+00 1.00e+00h  1
      13  4.3228323e-03 1.49e-08 2.51e-14  -8.6 7.00e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.3228322723800209e-03    4.3228322723800209e-03
    Dual infeasibility......:   2.5085432767634323e-14    2.5085432767634323e-14
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.045
    Total CPU secs in NLP function evaluations           =      0.041
    
    EXIT: Optimal Solution Found.
    Ipopt 3.13.2: max_iter=6000
    
    
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
       0  6.1350570e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  8.6637305e-02 1.79e+03 1.54e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  8.6384257e-03 6.37e+02 1.58e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  7.7740780e-03 1.81e+02 3.31e-03  -1.0 3.20e+00    -  9.95e-01 1.00e+00h  1
       4  7.3877602e-03 1.42e-01 8.46e-04  -1.7 2.61e+00    -  1.00e+00 1.00e+00h  1
       5  7.3459512e-03 1.97e+00 3.34e-04  -3.8 5.31e-01    -  1.00e+00 1.00e+00h  1
       6  7.3381344e-03 5.64e-02 1.22e-04  -3.8 3.37e-02  -4.0 1.00e+00 1.00e+00h  1
       7  5.2368049e-03 1.25e+04 6.40e-02  -5.7 1.58e+01    -  8.03e-01 1.00e+00h  1
       8  6.7607039e-03 6.89e+03 3.49e-02  -5.7 1.24e+01    -  8.67e-01 5.00e-01h  2
       9  7.0355639e-03 6.47e+03 3.23e-02  -5.7 4.10e+02    -  1.00e+00 6.25e-02h  5
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  7.2127426e-03 6.27e+03 3.10e-02  -5.7 4.93e+02    -  1.00e+00 3.12e-02h  6
      11  7.3084955e-03 6.17e+03 3.05e-02  -5.7 5.82e+02    -  1.00e+00 1.56e-02h  7
      12  7.3203565e-03 6.16e+03 3.04e-02  -5.7 7.00e+02    -  1.00e+00 1.95e-03h 10
      13  3.7972674e-02 6.83e+02 3.44e-01  -5.7 7.27e+02    -  1.00e+00 1.00e+00h  1
      14  5.2107669e-03 1.02e+02 1.48e-01  -5.7 7.77e+02    -  1.00e+00 1.00e+00h  1
      15  5.0043345e-03 1.38e+01 1.06e-03  -5.7 4.43e+00    -  1.00e+00 1.00e+00h  1
      16  5.0019078e-03 1.69e-01 2.70e-06  -5.7 2.18e-01    -  1.00e+00 1.00e+00h  1
      17  5.0018809e-03 1.15e-05 9.63e-11  -5.7 1.03e-03    -  1.00e+00 1.00e+00h  1
      18  5.0018807e-03 2.21e-05 1.14e-09  -8.6 6.41e-03    -  1.00e+00 1.00e+00h  1
      19  5.0018807e-03 1.49e-08 2.51e-14  -8.6 4.82e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 19
    
                                       (scaled)                 (unscaled)
    Objective...............:   5.0018806566479454e-03    5.0018806566479454e-03
    Dual infeasibility......:   2.5140943918865581e-14    2.5140943918865581e-14
    Constraint violation....:   1.4104644499482428e-11    1.4901161193847656e-08
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    1.4901161193847656e-08
    
    
    Number of objective function evaluations             = 45
    Number of objective gradient evaluations             = 20
    Number of equality constraint evaluations            = 45
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 20
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 19
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.063
    Total CPU secs in NLP function evaluations           =      0.080
    
    EXIT: Optimal Solution Found.
    Ipopt 3.13.2: max_iter=6000
    
    
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
       0  6.1794939e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  8.7233770e-02 1.79e+03 1.36e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  8.1206319e-03 6.21e+02 1.44e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  7.4626777e-03 1.66e+02 2.84e-03  -1.0 2.70e+00    -  9.95e-01 1.00e+00h  1
       4  7.1964953e-03 1.27e-01 6.47e-04  -1.7 1.21e+00    -  1.00e+00 1.00e+00h  1
       5  7.1648985e-03 1.62e+00 3.32e-04  -3.8 2.08e-01    -  1.00e+00 1.00e+00h  1
       6  7.1577644e-03 5.67e-02 1.15e-04  -3.8 3.43e-02  -4.0 1.00e+00 1.00e+00h  1
       7  5.2510992e-03 1.21e+04 6.87e-02  -5.7 1.57e+01    -  8.06e-01 1.00e+00h  1
       8  6.8576704e-03 6.70e+03 3.77e-02  -5.7 1.16e+01    -  8.71e-01 5.00e-01h  2
       9  7.0073487e-03 6.60e+03 3.71e-02  -5.7 8.92e+02    -  1.00e+00 1.56e-02h  7
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  7.1130854e-03 6.39e+03 8.35e-02  -5.7 9.38e+01  -4.5 1.00e+00 3.12e-02h  6
      11  7.1438030e-03 6.34e+03 8.99e-02  -5.7 1.87e+02  -5.0 1.00e+00 7.81e-03h  8
      12  7.1474475e-03 6.34e+03 9.00e-02  -5.7 1.39e+02  -4.5 1.00e+00 9.77e-04h 11
      13  1.4456043e-02 1.56e+03 3.53e+01  -5.7 2.66e+02  -5.0 1.00e+00 1.00e+00h  1
      14  5.4394950e-03 4.07e+01 7.03e+00  -5.7 2.84e+02    -  1.00e+00 1.00e+00h  1
      15  5.0111912e-03 3.33e+01 5.74e-02  -5.7 6.85e+00    -  1.00e+00 1.00e+00h  1
      16  4.9961298e-03 2.80e+00 6.89e-05  -5.7 1.53e+00    -  1.00e+00 1.00e+00h  1
      17  4.9956073e-03 1.45e-02 3.13e-07  -5.7 9.67e-02    -  1.00e+00 1.00e+00h  1
      18  4.9956049e-03 2.04e-07 1.84e-11  -5.7 2.76e-04    -  1.00e+00 1.00e+00h  1
      19  4.9956047e-03 2.44e-05 1.14e-09  -8.6 7.04e-03    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      20  4.9956047e-03 1.49e-08 2.51e-14  -8.6 5.62e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 20
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.9956046942981613e-03    4.9956046942981613e-03
    Dual infeasibility......:   2.5085432767634323e-14    2.5085432767634323e-14
    Constraint violation....:   1.4104644499482428e-11    1.4901161193847656e-08
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    1.4901161193847656e-08
    
    
    Number of objective function evaluations             = 50
    Number of objective gradient evaluations             = 21
    Number of equality constraint evaluations            = 50
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 21
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 20
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.072
    Total CPU secs in NLP function evaluations           =      0.082
    
    EXIT: Optimal Solution Found.
    Ipopt 3.13.2: max_iter=6000
    
    
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
       0  6.9524347e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  1.0476857e-01 1.79e+03 1.48e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  1.0144861e-02 1.01e+03 1.92e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  8.7236266e-03 5.14e+02 7.32e-03  -1.0 4.66e+00    -  9.95e-01 1.00e+00h  1
       4  8.2629544e-03 2.67e+01 3.34e-04  -1.0 1.78e+00    -  1.00e+00 1.00e+00h  1
       5  8.3595735e-03 7.43e-02 3.36e-04  -2.5 2.39e+00    -  1.00e+00 1.00e+00h  1
       6  8.0322662e-03 1.08e+02 1.14e-03  -3.8 1.44e+00    -  1.00e+00 1.00e+00h  1
       7  8.0205315e-03 1.98e+00 9.82e-04  -3.8 4.80e-02  -4.0 1.00e+00 1.00e+00h  1
       8  6.4208242e-03 1.93e+03 2.25e-02  -3.8 8.62e+00    -  1.00e+00 1.00e+00h  1
       9  5.4773293e-03 1.60e+02 4.87e-03  -3.8 1.28e+01    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  5.3943280e-03 6.77e-01 3.50e-05  -3.8 1.14e+00    -  1.00e+00 1.00e+00h  1
      11  5.3940228e-03 3.85e-03 2.20e-08  -3.8 1.81e-02    -  1.00e+00 1.00e+00h  1
      12  5.3927724e-03 1.18e-01 6.36e-06  -5.7 4.98e-01    -  1.00e+00 1.00e+00h  1
      13  5.3927548e-03 1.90e-05 3.35e-10  -5.7 2.71e-03    -  1.00e+00 1.00e+00h  1
      14  5.3927546e-03 1.72e-05 9.48e-10  -8.6 6.11e-03    -  1.00e+00 1.00e+00h  1
      15  5.3927546e-03 1.49e-08 2.51e-14  -8.6 4.01e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 15
    
                                       (scaled)                 (unscaled)
    Objective...............:   5.3927546365642610e-03    5.3927546365642610e-03
    Dual infeasibility......:   2.5140943918865581e-14    2.5140943918865581e-14
    Constraint violation....:   1.4104644499482428e-11    1.4901161193847656e-08
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    1.4901161193847656e-08
    
    
    Number of objective function evaluations             = 16
    Number of objective gradient evaluations             = 16
    Number of equality constraint evaluations            = 16
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 16
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 15
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.048
    Total CPU secs in NLP function evaluations           =      0.048
    
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
          <td>0.476675</td>
          <td>-0.408169</td>
        </tr>
        <tr>
          <th>1</th>
          <td>0.484354</td>
          <td>-0.413949</td>
        </tr>
        <tr>
          <th>2</th>
          <td>0.472694</td>
          <td>-0.404994</td>
        </tr>
        <tr>
          <th>3</th>
          <td>0.468290</td>
          <td>-0.401568</td>
        </tr>
        <tr>
          <th>4</th>
          <td>0.473732</td>
          <td>-0.405595</td>
        </tr>
        <tr>
          <th>5</th>
          <td>0.449693</td>
          <td>-0.387192</td>
        </tr>
        <tr>
          <th>6</th>
          <td>0.464854</td>
          <td>-0.398884</td>
        </tr>
        <tr>
          <th>7</th>
          <td>0.499482</td>
          <td>-0.425310</td>
        </tr>
        <tr>
          <th>8</th>
          <td>0.495345</td>
          <td>-0.422358</td>
        </tr>
        <tr>
          <th>9</th>
          <td>0.517508</td>
          <td>-0.438912</td>
        </tr>
      </tbody>
    </table>
    </div>

