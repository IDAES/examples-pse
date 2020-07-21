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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.173
    Total CPU secs in NLP function evaluations           =      0.135
    
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
       0  5.7090216e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  8.3311134e-02 1.79e+03 1.47e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  7.7435193e-03 6.11e+02 1.50e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  7.0528187e-03 1.58e+02 2.90e-03  -1.0 2.86e+00    -  9.95e-01 1.00e+00h  1
       4  6.7726273e-03 1.07e-01 6.93e-04  -1.7 1.89e+00    -  1.00e+00 1.00e+00h  1
       5  6.7414115e-03 1.53e+00 3.05e-04  -3.8 3.13e-01    -  1.00e+00 1.00e+00h  1
       6  6.7348587e-03 5.36e-02 1.07e-04  -3.8 3.20e-02  -4.0 1.00e+00 1.00e+00h  1
       7  5.0065562e-03 1.03e+04 5.60e-02  -5.7 1.42e+01    -  8.19e-01 1.00e+00h  1
       8  5.9452768e-03 5.73e+03 3.09e-02  -5.7 9.11e+00    -  8.91e-01 5.00e-01h  2
       9  5.4971690e-03 1.34e+03 1.19e-02  -5.7 4.38e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.7851029e-03 5.37e+00 1.01e-03  -5.7 6.10e+00    -  1.00e+00 1.00e+00h  1
      11  4.7630331e-03 4.38e-03 1.31e-06  -5.7 4.81e-01    -  1.00e+00 1.00e+00h  1
      12  4.7630374e-03 4.31e-05 2.00e-09  -8.6 8.94e-03    -  1.00e+00 1.00e+00h  1
      13  4.7630374e-03 7.45e-09 2.51e-14  -8.6 1.22e-06    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.7630374068821117e-03    4.7630374068821117e-03
    Dual infeasibility......:   2.5085432767634323e-14    2.5085432767634323e-14
    Constraint violation....:   1.4104644499482428e-11    7.4505805969238281e-09
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    7.4505805969238281e-09
    
    
    Number of objective function evaluations             = 15
    Number of objective gradient evaluations             = 14
    Number of equality constraint evaluations            = 15
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 14
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 13
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.158
    Total CPU secs in NLP function evaluations           =      0.221
    
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
       0  5.5931959e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  8.2068816e-02 1.79e+03 1.38e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  7.5614433e-03 8.76e+02 1.73e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  6.7535311e-03 3.44e+02 5.00e-03  -1.0 3.46e+00    -  9.95e-01 1.00e+00h  1
       4  6.4404576e-03 1.06e+01 1.50e-04  -1.0 1.64e+00    -  1.00e+00 1.00e+00h  1
       5  6.4856227e-03 4.35e-02 2.34e-04  -2.5 4.45e-01    -  1.00e+00 1.00e+00h  1
       6  6.3090394e-03 5.71e+01 6.64e-04  -3.8 1.05e+00    -  1.00e+00 1.00e+00h  1
       7  6.3020694e-03 2.47e-01 2.95e-04  -3.8 3.82e-02  -4.0 1.00e+00 1.00e+00h  1
       8  5.1601672e-03 2.59e+03 2.03e-02  -5.7 8.61e+00    -  9.06e-01 1.00e+00h  1
       9  4.7938424e-03 3.79e+02 5.78e-03  -5.7 8.20e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.6523482e-03 1.67e-01 1.09e-04  -5.7 2.67e+00    -  1.00e+00 1.00e+00h  1
      11  4.6507994e-03 4.19e-02 1.08e-06  -5.7 1.92e-01    -  1.00e+00 1.00e+00h  1
      12  4.6507941e-03 4.98e-06 9.80e-11  -5.7 1.61e-03    -  1.00e+00 1.00e+00h  1
      13  4.6507939e-03 3.11e-05 1.35e-09  -8.6 7.26e-03    -  1.00e+00 1.00e+00h  1
      14  4.6507939e-03 1.49e-08 2.51e-14  -8.6 6.64e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 14
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.6507938670864403e-03    4.6507938670864403e-03
    Dual infeasibility......:   2.5140943918865581e-14    2.5140943918865581e-14
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.101
    Total CPU secs in NLP function evaluations           =      0.236
    
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
       0  4.3896291e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  5.6747095e-02 1.79e+03 1.34e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  5.1199732e-03 2.13e+02 8.35e-02  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  4.9812505e-03 1.34e-01 1.25e-04  -1.0 9.98e-01    -  9.95e-01 1.00e+00h  1
       4  4.9672666e-03 2.35e-02 1.93e-04  -2.5 6.88e-01    -  1.00e+00 1.00e+00h  1
       5  4.8763863e-03 2.93e+01 4.16e-04  -3.8 7.86e-01    -  1.00e+00 1.00e+00h  1
       6  4.8725692e-03 7.76e-03 5.43e-05  -3.8 2.94e-02  -4.0 1.00e+00 1.00e+00h  1
       7  4.1284115e-03 3.00e+03 1.97e-02  -5.7 8.58e+00    -  8.98e-01 1.00e+00h  1
       8  4.0972663e-03 5.59e+02 5.31e-03  -5.7 4.23e+00    -  1.00e+00 1.00e+00h  1
       9  3.9069185e-03 4.35e-01 2.21e-04  -5.7 3.64e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  3.9036366e-03 4.41e-02 1.77e-06  -5.7 3.18e-01    -  1.00e+00 1.00e+00h  1
      11  3.9036326e-03 1.71e-05 3.59e-10  -5.7 3.50e-03    -  1.00e+00 1.00e+00h  1
      12  3.9036323e-03 5.88e-05 1.98e-09  -8.6 8.92e-03    -  1.00e+00 1.00e+00h  1
      13  3.9036323e-03 7.45e-09 2.51e-14  -8.6 1.19e-06    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   3.9036323435068187e-03    3.9036323435068187e-03
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.113
    Total CPU secs in NLP function evaluations           =      0.193
    
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
       0  4.2470579e+01 4.57e+05 3.36e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  4.7730476e-02 1.18e+03 1.19e+00  -1.0 1.12e+04    -  9.83e-01 1.00e+00h  1
       2  5.8438802e-03 1.33e+02 6.06e-02  -1.0 3.12e+02    -  9.90e-01 1.00e+00h  1
       3  5.6819152e-03 7.35e+00 1.66e-04  -1.0 7.88e-01    -  9.96e-01 1.00e+00h  1
       4  5.2871703e-03 5.11e-02 3.52e-04  -2.5 5.38e+00    -  1.00e+00 1.00e+00h  1
       5  5.1426623e-03 4.70e+01 5.02e-04  -3.8 8.52e-01    -  1.00e+00 1.00e+00h  1
       6  5.1375298e-03 7.49e-02 1.35e-04  -3.8 3.04e-02  -4.0 1.00e+00 1.00e+00h  1
       7  4.2153308e-03 2.61e+03 1.77e-02  -5.7 8.03e+00    -  9.05e-01 1.00e+00h  1
       8  4.0560471e-03 4.26e+02 4.95e-03  -5.7 5.28e+00    -  1.00e+00 1.00e+00h  1
       9  3.9223613e-03 2.32e-01 1.08e-04  -5.7 2.48e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  3.9206100e-03 3.32e-02 1.04e-06  -5.7 1.96e-01    -  1.00e+00 1.00e+00h  1
      11  3.9206067e-03 5.65e-06 1.13e-10  -5.7 1.70e-03    -  1.00e+00 1.00e+00h  1
      12  3.9206065e-03 4.44e-05 1.64e-09  -8.6 7.18e-03    -  1.00e+00 1.00e+00h  1
      13  3.9206065e-03 7.45e-09 2.51e-14  -8.6 8.11e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   3.9206064678035529e-03    3.9206064678035529e-03
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.122
    Total CPU secs in NLP function evaluations           =      0.195
    
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
       0  6.6124320e+01 5.41e+05 4.00e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  1.0445594e-01 1.66e+03 1.54e+00  -1.0 1.32e+04    -  9.82e-01 1.00e+00h  1
       2  9.3378449e-03 2.88e+03 3.20e-01  -1.0 4.39e+02    -  9.90e-01 1.00e+00h  1
       3  5.7961851e-03 2.91e+03 6.94e-02  -1.0 2.35e+01    -  9.96e-01 1.00e+00h  1
       4  5.5473081e-03 2.86e+02 2.62e-03  -1.0 2.85e+00    -  1.00e+00 1.00e+00h  1
       5  7.0594142e-03 2.82e+00 6.41e-04  -1.0 1.27e+01    -  1.00e+00 1.00e+00h  1
       6  6.8544551e-03 3.64e-02 2.48e-04  -2.5 3.72e+00    -  1.00e+00 1.00e+00h  1
       7  6.6900324e-03 5.24e+01 7.16e-04  -3.8 1.19e+00    -  1.00e+00 1.00e+00h  1
       8  6.6827310e-03 1.95e-01 3.09e-04  -3.8 4.27e-02  -4.0 1.00e+00 1.00e+00h  1
       9  5.5031737e-03 2.88e+03 2.52e-02  -5.7 9.47e+00    -  9.00e-01 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  5.1339931e-03 4.49e+02 6.71e-03  -5.7 9.24e+00    -  1.00e+00 1.00e+00h  1
      11  4.9277143e-03 3.21e-01 1.86e-04  -5.7 3.61e+00    -  1.00e+00 1.00e+00h  1
      12  4.9251641e-03 6.21e-02 1.56e-06  -5.7 2.83e-01    -  1.00e+00 1.00e+00h  1
      13  4.9251543e-03 1.21e-05 2.23e-10  -5.7 2.96e-03    -  1.00e+00 1.00e+00h  1
      14  4.9251541e-03 2.93e-05 1.13e-09  -8.6 8.05e-03    -  1.00e+00 1.00e+00h  1
      15  4.9251541e-03 1.49e-08 2.51e-14  -8.6 7.12e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 15
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.9251541053267249e-03    4.9251541053267249e-03
    Dual infeasibility......:   2.5085432767634323e-14    2.5085432767634323e-14
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.251
    Total CPU secs in NLP function evaluations           =      0.216
    
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
       0  6.2185191e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  1.1081138e-01 1.79e+03 1.56e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  8.3714961e-03 2.72e+03 3.11e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  5.4428059e-03 2.59e+03 5.96e-02  -1.0 2.04e+01    -  9.95e-01 1.00e+00h  1
       4  5.3038024e-03 2.24e+02 2.02e-03  -1.0 2.26e+00    -  1.00e+00 1.00e+00h  1
       5  6.2481688e-03 1.30e-01 8.71e-04  -1.7 8.86e+00    -  1.00e+00 1.00e+00h  1
       6  6.2185374e-03 9.63e-01 2.98e-04  -3.8 9.99e-01    -  1.00e+00 1.00e+00h  1
       7  5.3952778e-03 2.34e+03 1.78e-02  -3.8 6.93e+00    -  1.00e+00 1.00e+00h  1
       8  5.7858129e-03 1.35e+03 1.05e-02  -3.8 6.99e+00    -  1.00e+00 5.00e-01h  2
       9  5.7410491e-03 2.95e+02 8.02e-02  -3.8 3.23e-01  -4.0 1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  5.6970234e-03 2.24e-01 3.94e-03  -3.8 3.26e-01  -4.5 1.00e+00 1.00e+00h  1
      11  4.8331955e-03 6.78e+02 3.30e-02  -3.8 9.84e+00    -  1.00e+00 1.00e+00h  1
      12  4.6444269e-03 3.34e+01 8.09e-04  -3.8 6.36e+00    -  1.00e+00 1.00e+00h  1
      13  4.6317329e-03 6.35e-02 2.81e-06  -3.8 4.28e-01    -  1.00e+00 1.00e+00h  1
      14  4.6300562e-03 2.55e-01 9.13e-06  -5.7 7.01e-01    -  1.00e+00 1.00e+00h  1
      15  4.6300263e-03 1.05e-04 1.45e-09  -5.7 5.90e-03    -  1.00e+00 1.00e+00h  1
      16  4.6300261e-03 3.71e-05 1.36e-09  -8.6 8.58e-03    -  1.00e+00 1.00e+00h  1
      17  4.6300261e-03 1.49e-08 2.51e-14  -8.6 8.70e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 17
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.6300260603010070e-03    4.6300260603010070e-03
    Dual infeasibility......:   2.5085432767634323e-14    2.5085432767634323e-14
    Constraint violation....:   1.4104644499482428e-11    1.4901161193847656e-08
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    1.4901161193847656e-08
    
    
    Number of objective function evaluations             = 19
    Number of objective gradient evaluations             = 18
    Number of equality constraint evaluations            = 19
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 18
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 17
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.298
    Total CPU secs in NLP function evaluations           =      0.200
    
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
       0  4.4358477e+01 5.20e+05 3.84e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  5.1535792e-02 1.53e+03 1.37e+00  -1.0 1.27e+04    -  9.82e-01 1.00e+00h  1
       2  5.8630753e-03 2.01e+02 8.21e-02  -1.0 4.05e+02    -  9.90e-01 1.00e+00h  1
       3  5.6293774e-03 9.38e-02 1.10e-04  -1.0 1.22e+00    -  9.96e-01 1.00e+00h  1
       4  5.4780122e-03 4.70e-02 2.34e-04  -2.5 3.19e+00    -  1.00e+00 1.00e+00h  1
       5  5.3366802e-03 4.62e+01 5.06e-04  -3.8 8.05e-01    -  1.00e+00 1.00e+00h  1
       6  5.3314380e-03 9.79e-02 1.54e-04  -3.8 3.12e-02  -4.0 1.00e+00 1.00e+00h  1
       7  4.3970615e-03 2.62e+03 1.82e-02  -5.7 8.09e+00    -  9.05e-01 1.00e+00h  1
       8  4.2155532e-03 4.25e+02 5.22e-03  -5.7 5.68e+00    -  1.00e+00 1.00e+00h  1
       9  4.0792941e-03 2.83e-01 1.10e-04  -5.7 2.56e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.0775337e-03 3.75e-02 1.18e-06  -5.7 2.04e-01    -  1.00e+00 1.00e+00h  1
      11  4.0775299e-03 6.57e-06 1.36e-10  -5.7 1.84e-03    -  1.00e+00 1.00e+00h  1
      12  4.0775296e-03 4.29e-05 1.68e-09  -8.6 7.31e-03    -  1.00e+00 1.00e+00h  1
      13  4.0775296e-03 7.45e-09 2.51e-14  -8.6 8.09e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.0775296347024337e-03    4.0775296347024337e-03
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.198
    Total CPU secs in NLP function evaluations           =      0.172
    
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
       0  5.8523477e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  8.6872497e-02 1.79e+03 1.40e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  7.6410008e-03 5.45e+02 1.37e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  7.0603941e-03 1.18e+02 2.27e-03  -1.0 2.50e+00    -  9.95e-01 1.00e+00h  1
       4  6.8458093e-03 5.43e-02 4.86e-04  -1.7 1.14e+00    -  1.00e+00 1.00e+00h  1
       5  6.8158313e-03 1.47e+00 3.11e-04  -3.8 1.66e-01    -  1.00e+00 1.00e+00h  1
       6  6.8092792e-03 5.60e-02 1.09e-04  -3.8 3.28e-02  -4.0 1.00e+00 1.00e+00h  1
       7  5.0859703e-03 1.06e+04 5.98e-02  -5.7 1.44e+01    -  8.17e-01 1.00e+00h  1
       8  6.1791840e-03 5.88e+03 3.31e-02  -5.7 9.56e+00    -  8.88e-01 5.00e-01h  2
       9  5.6742736e-03 1.39e+03 1.27e-02  -5.7 4.62e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.8501736e-03 4.52e+00 1.16e-03  -5.7 6.60e+00    -  1.00e+00 1.00e+00h  1
      11  4.8241336e-03 1.07e-04 8.41e-07  -5.7 5.45e-01    -  1.00e+00 1.00e+00h  1
      12  4.8241394e-03 2.53e-05 1.15e-09  -8.6 6.90e-03    -  1.00e+00 1.00e+00h  1
      13  4.8241394e-03 1.49e-08 2.51e-14  -8.6 5.40e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.8241394280379895e-03    4.8241394280379895e-03
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.175
    Total CPU secs in NLP function evaluations           =      0.160
    
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
       0  4.6050524e+01 5.41e+05 4.00e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  5.6479690e-02 1.66e+03 1.53e+00  -1.0 1.32e+04    -  9.82e-01 1.00e+00h  1
       2  6.1587989e-03 1.92e+02 8.65e-02  -1.0 4.39e+02    -  9.90e-01 1.00e+00h  1
       3  5.9131736e-03 2.05e-02 1.22e-04  -1.0 1.22e+00    -  9.96e-01 1.00e+00h  1
       4  5.6892114e-03 5.48e-02 2.69e-04  -2.5 3.95e+00    -  1.00e+00 1.00e+00h  1
       5  5.5300646e-03 5.22e+01 5.72e-04  -3.8 8.68e-01    -  1.00e+00 1.00e+00h  1
       6  5.5242965e-03 1.57e-01 2.06e-04  -3.8 3.34e-02  -4.0 1.00e+00 1.00e+00h  1
       7  4.5275302e-03 2.51e+03 1.94e-02  -5.7 8.21e+00    -  9.07e-01 1.00e+00h  1
       8  4.2923630e-03 3.86e+02 5.43e-03  -5.7 6.34e+00    -  1.00e+00 1.00e+00h  1
       9  4.1666900e-03 2.74e-01 9.95e-05  -5.7 2.46e+00    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.1652292e-03 3.60e-02 1.07e-06  -5.7 1.83e-01    -  1.00e+00 1.00e+00h  1
      11  4.1652253e-03 4.88e-06 1.03e-10  -5.7 1.55e-03    -  1.00e+00 1.00e+00h  1
      12  4.1652251e-03 3.89e-05 1.66e-09  -8.6 7.20e-03    -  1.00e+00 1.00e+00h  1
      13  4.1652251e-03 7.45e-09 2.51e-14  -8.6 7.51e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.1652250871936367e-03    4.1652250871936367e-03
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.130
    Total CPU secs in NLP function evaluations           =      0.172
    
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
       0  5.4911681e+01 5.63e+05 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  7.8160707e-02 1.79e+03 1.32e+00  -1.0 1.37e+04    -  9.82e-01 1.00e+00h  1
       2  6.9365117e-03 5.98e+02 1.39e-01  -1.0 4.74e+02    -  9.90e-01 1.00e+00h  1
       3  6.4373281e-03 1.43e+02 2.50e-03  -1.0 2.41e+00    -  9.95e-01 1.00e+00h  1
       4  6.2194344e-03 9.18e-02 5.27e-04  -1.7 1.72e+00    -  1.00e+00 1.00e+00h  1
       5  6.1969348e-03 1.11e+00 2.67e-04  -3.8 3.23e-01    -  1.00e+00 1.00e+00h  1
       6  6.1915911e-03 4.92e-02 9.29e-05  -3.8 3.00e-02  -4.0 1.00e+00 1.00e+00h  1
       7  4.8186156e-03 8.37e+03 4.75e-02  -5.7 1.25e+01    -  8.35e-01 1.00e+00h  1
       8  6.0766611e-03 1.93e+03 1.44e-02  -5.7 6.67e+00    -  9.14e-01 1.00e+00h  1
       9  4.6552471e-03 1.63e+02 3.00e-03  -5.7 1.06e+01    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.5650243e-03 6.36e+00 8.27e-05  -5.7 8.59e-01    -  1.00e+00 1.00e+00h  1
      11  4.5642883e-03 1.09e-02 4.70e-07  -5.7 1.69e-01    -  1.00e+00 1.00e+00h  1
      12  4.5642854e-03 5.48e-08 1.84e-11  -5.7 4.48e-04    -  1.00e+00 1.00e+00h  1
      13  4.5642852e-03 3.45e-05 1.41e-09  -8.6 7.64e-03    -  1.00e+00 1.00e+00h  1
      14  4.5642852e-03 1.49e-08 2.51e-14  -8.6 7.44e-07    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 14
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.5642851869240279e-03    4.5642851869240279e-03
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.135
    Total CPU secs in NLP function evaluations           =      0.242
    
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
          <td>0.485639</td>
          <td>-0.414888</td>
        </tr>
        <tr>
          <th>1</th>
          <td>0.477456</td>
          <td>-0.408739</td>
        </tr>
        <tr>
          <th>2</th>
          <td>0.433285</td>
          <td>-0.374654</td>
        </tr>
        <tr>
          <th>3</th>
          <td>0.447786</td>
          <td>-0.385511</td>
        </tr>
        <tr>
          <th>4</th>
          <td>0.486034</td>
          <td>-0.415542</td>
        </tr>
        <tr>
          <th>5</th>
          <td>0.469280</td>
          <td>-0.402769</td>
        </tr>
        <tr>
          <th>6</th>
          <td>0.451944</td>
          <td>-0.388895</td>
        </tr>
        <tr>
          <th>7</th>
          <td>0.487185</td>
          <td>-0.416137</td>
        </tr>
        <tr>
          <th>8</th>
          <td>0.458341</td>
          <td>-0.393791</td>
        </tr>
        <tr>
          <th>9</th>
          <td>0.471282</td>
          <td>-0.404081</td>
        </tr>
      </tbody>
    </table>
    </div>

