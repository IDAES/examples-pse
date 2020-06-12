Module 4b: Parameter Estimation Using the NRTL State Block
----------------------------------------------------------

In this module, we use Pyomo’s ``parmest`` tool in conjunction with
IDAES models for parameter estimation. We demonstrate these tools by
estimating the parameters associated with the NRTL property model for a
benzene-toluene mixture. The NRTL model has 2 sets of parameters: the
non-randomness parameter (``alpha_ij``) and the binary interaction
parameter (``tau_ij``), where ``i`` and ``j`` are the pure component
species. In this example, we only estimate the binary interaction
parameter (``tau_ij``) for a given dataset. When estimating parameters
associated with the property package, IDAES provides the flexibility of
doing the parameter estimation by just using the state block or by using
a unit model with a specified property package. This module will
demonstrate parameter estimation by using only the state block.

We will complete the following tasks: \* Set up a method to return an
initialized model \* Set up the parameter estimation problem using
``parmest`` \* Analyze the results \* Demonstrate advanced features
using ``parmest``

Key links to documentation:
---------------------------

-  NRTL Model -
   https://idaes-pse.readthedocs.io/en/latest/model_libraries/core_library/property_models/activity_coefficient.html
-  parmest -
   https://pyomo.readthedocs.io/en/stable/contributed_packages/parmest/index.html

.. raw:: html

   <div class="alert alert-block alert-info">

Inline Exercise: import ``ConcreteModel`` from Pyomo and
``FlowsheetBlock`` from IDAES.

.. raw:: html

   </div>

.. code:: ipython3

    # Todo: import ConcreteModel from pyomo.environ
    from pyomo.environ import ConcreteModel, value
    
    # Todo: import FlowsheetBlock from idaes.core
    from idaes.core import FlowsheetBlock

In the next cell, we import the parameter block used in this module and
the idaes logger.

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
        m.fs.state_block = m.fs.properties.state_block_class(
            default={"parameters": m.fs.properties,
                     "defined_state": True})
    
        
        # Todo: Fix the state varaibles on the state block
        # hint: state variables exist on the state block i.e. on m.fs.state_block
        
        m.fs.state_block.flow_mol.fix(1)
        m.fs.state_block.temperature.fix(368)
        m.fs.state_block.pressure.fix(101325)
        m.fs.state_block.mole_frac_comp["benzene"].fix(0.5)
        m.fs.state_block.mole_frac_comp["toluene"].fix(0.5)
    
        # Fix NRTL specific parameters. 
        
        # non-randomness parameter - alpha_ij (set at 0.3, 0 if i=j)
        m.fs.properties.\
            alpha["benzene", "benzene"].fix(0)
        m.fs.properties.\
            alpha["benzene", "toluene"].fix(0.3)
        m.fs.properties.\
            alpha["toluene", "toluene"].fix(0)
        m.fs.properties.\
            alpha["toluene", "benzene"].fix(0.3)
    
        # binary interaction parameter - tau_ij (0 if i=j, else to be estimated later but fixing to initialize)
        m.fs.properties.\
            tau["benzene", "benzene"].fix(0)
        m.fs.properties.\
            tau["benzene", "toluene"].fix(0.1690)
        m.fs.properties.\
            tau["toluene", "toluene"].fix(0)
        m.fs.properties.\
            tau["toluene", "benzene"].fix(-0.1559)
    
        # Initialize the flash unit
        m.fs.state_block.initialize(outlvl=idaeslog.INFO)
    
        # Fix at actual temperature
        m.fs.state_block.temperature.fix(float(data["temperature"]))
    
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
    assert value(m.fs.state_block.mole_frac_phase_comp['Liq', 'benzene']) == pytest.approx(0.4105, abs=1e-3)
    assert value(m.fs.state_block.mole_frac_phase_comp['Vap', 'benzene']) == pytest.approx(0.6326, abs=1e-3)
    
    assert value(m.fs.state_block.mole_frac_phase_comp['Liq', 'toluene']) == pytest.approx(0.5895, abs=1e-3)
    assert value(m.fs.state_block.mole_frac_phase_comp['Vap', 'toluene']) == pytest.approx(0.3673, abs=1e-3)


.. parsed-literal::

    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    

Parameter estimation using parmest
----------------------------------

In addition to providing a method to return an initialized model, the
``parmest`` tool needs the following:

-  List of variable names to be estimated
-  Dataset
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
        # (float(data["vap_benzene"]) - m.fs.state_block.mole_frac_phase_comp["Vap", "benzene"])**2
        expr = ((float(data["vap_benzene"]) -
                 m.fs.state_block.mole_frac_phase_comp["Vap", "benzene"])**2 +
                (float(data["liq_benzene"]) -
                 m.fs.state_block.mole_frac_phase_comp["Liq", "benzene"])**2)
        return expr*1E4

.. raw:: html

   <div class="alert alert-block alert-warning">

Note: Notice that we have scaled the expression up by a factor of 10000
as the SSE computed here will be an extremely small number given that we
are using the difference in mole fraction in our expression. This will
help in using a well-scaled objective to improve solve robustness when
using IPOPT.

.. raw:: html

   </div>

We are now ready to set up the parameter estimation problem. We will
create a parameter estimation object called ``pest``. As shown below, we
pass the method that returns an initialized model, data, variable_name,
and the SSE expression to the Estimator method. ``tee=True`` will print
the solver output after solving the parameter estimation problem.

.. code:: ipython3

    # Initialize a parameter estimation object
    pest = parmest.Estimator(NRTL_model, data, variable_name, SSE, tee=True)
    
    # Run parameter estimation using all data
    obj_value, parameters = pest.theta_est()


.. parsed-literal::

    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
    
    Number of nonzeros in equality constraint Jacobian...:     3750
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:     2200
    
    Total number of variables............................:     1102
                         variables with only lower bounds:        0
                    variables with lower and upper bounds:      300
                         variables with only upper bounds:        0
    Total number of equality constraints.................:     1100
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  5.5857491e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  7.8156477e-02 1.41e+03 4.83e-01  -1.0 1.37e+04    -  9.88e-01 1.00e+00f  1
       2  6.1697942e-03 1.10e+01 1.73e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  6.1984875e-03 8.87e-02 3.38e-05  -2.5 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  4.8761541e-03 2.95e+03 2.50e-02  -3.8 6.42e-01    -  9.33e-01 1.00e+00h  1
       5  5.3296404e-03 7.11e+02 4.21e-03  -3.8 3.17e-01    -  1.00e+00 1.00e+00h  1
       6  4.7282530e-03 1.22e+01 2.04e-03  -3.8 6.72e-02    -  1.00e+00 1.00e+00h  1
       7  4.6651516e-03 1.10e+01 3.42e-04  -3.8 4.35e-02    -  1.00e+00 1.00e+00h  1
       8  4.6648092e-03 3.85e-01 7.85e-06  -3.8 7.85e-03    -  1.00e+00 1.00e+00h  1
       9  4.6633709e-03 1.85e-01 7.65e-06  -5.7 5.64e-03    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.6633491e-03 5.52e-05 7.34e-10  -5.7 9.26e-05    -  1.00e+00 1.00e+00h  1
      11  4.6633488e-03 2.83e-05 1.18e-09  -8.6 6.98e-05    -  1.00e+00 1.00e+00h  1
      12  4.6633488e-03 8.73e-11 2.03e-14  -8.6 1.42e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 12
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.6633488370413792e-03    4.6633488370413792e-03
    Dual infeasibility......:   2.0301671375356738e-14    2.0301671375356738e-14
    Constraint violation....:   3.2290069315681125e-13    8.7311491370201111e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 13
    Number of objective gradient evaluations             = 13
    Number of equality constraint evaluations            = 13
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 13
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 12
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.015
    Total CPU secs in NLP function evaluations           =      0.014
    
    EXIT: Optimal Solution Found.
    

.. code:: ipython3

    # Check for values of the parameter estimation problem
    assert obj_value == pytest.approx(0.004663, 1e-3)
    assert parameters["fs.properties.tau[('benzene', 'toluene')]"] == pytest.approx(0.47811, 1e-3) 
    assert parameters["fs.properties.tau[('toluene', 'benzene')]"] == pytest.approx(-0.40924, 1e-3)

You will notice that the resulting parameter estimation problem will
have 1102 variables and 1100 constraints. Let us display the results by
running the next cell.

.. code:: ipython3

    print("The SSE at the optimal solution is %0.6f" % obj_value)
    print()
    print("The values for the parameters are as follows:")
    for k,v in parameters.items():
        print(k, "=", v)


.. parsed-literal::

    The SSE at the optimal solution is 0.004663
    
    The values for the parameters are as follows:
    fs.properties.tau[('benzene', 'toluene')] = 0.47810868272725465
    fs.properties.tau[('toluene', 'benzene')] = -0.4092446570740113
    

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

    2020-06-11 20:05:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:05:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
    
    Number of nonzeros in equality constraint Jacobian...:     3750
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:     2200
    
    Total number of variables............................:     1102
                         variables with only lower bounds:        0
                    variables with lower and upper bounds:      300
                         variables with only upper bounds:        0
    Total number of equality constraints.................:     1100
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  6.0061000e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  1.0202933e-01 1.41e+03 4.97e-01  -1.0 1.37e+04    -  9.88e-01 1.00e+00f  1
       2  6.2140443e-03 2.81e+01 1.96e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  6.2517079e-03 1.24e-01 4.22e-05  -2.5 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  5.4589717e-03 8.87e+02 8.04e-03  -3.8 3.50e-01    -  9.85e-01 1.00e+00h  1
       5  5.1671497e-03 8.30e+01 7.37e-03  -3.8 1.43e-01    -  1.00e+00 1.00e+00h  1
       6  4.7275762e-03 5.18e+01 1.72e-03  -3.8 9.51e-02    -  1.00e+00 1.00e+00h  1
       7  4.7065145e-03 6.75e+00 1.37e-04  -3.8 3.30e-02    -  1.00e+00 1.00e+00h  1
       8  4.7062467e-03 6.76e-02 1.05e-06  -3.8 3.23e-03    -  1.00e+00 1.00e+00h  1
       9  4.7047843e-03 2.06e-01 7.78e-06  -5.7 5.94e-03    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.7047608e-03 6.79e-05 9.46e-10  -5.7 1.03e-04    -  1.00e+00 1.00e+00h  1
      11  4.7047606e-03 2.97e-05 1.15e-09  -8.6 7.14e-05    -  1.00e+00 1.00e+00h  1
      12  4.7047606e-03 4.37e-11 2.22e-14  -8.6 1.49e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 12
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.7047606108040355e-03    4.7047606108040355e-03
    Dual infeasibility......:   2.2181805979252078e-14    2.2181805979252078e-14
    Constraint violation....:   1.6145034657840562e-13    4.3655745685100555e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 13
    Number of objective gradient evaluations             = 13
    Number of equality constraint evaluations            = 13
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 13
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 12
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.012
    Total CPU secs in NLP function evaluations           =      0.013
    
    EXIT: Optimal Solution Found.
    2020-06-11 20:06:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
    
    Number of nonzeros in equality constraint Jacobian...:     3750
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:     2200
    
    Total number of variables............................:     1102
                         variables with only lower bounds:        0
                    variables with lower and upper bounds:      300
                         variables with only upper bounds:        0
    Total number of equality constraints.................:     1100
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  4.9477362e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  7.2778616e-02 1.42e+03 4.56e-01  -1.0 1.37e+04    -  9.87e-01 1.00e+00f  1
       2  5.6282695e-03 7.16e+00 1.91e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  5.5862977e-03 6.09e+00 3.52e-04  -3.8 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  5.5686440e-03 3.58e-02 2.46e-05  -3.8 2.46e-03  -2.0 1.00e+00 1.00e+00h  1
       5  5.5410457e-03 3.99e-01 2.60e-05  -5.7 7.62e-03  -2.5 1.00e+00 1.00e+00h  1
       6  5.4455933e-03 4.55e+00 5.10e-05  -5.7 2.58e-02  -3.0 1.00e+00 1.00e+00h  1
       7  5.0056872e-03 9.21e+01 9.42e-04  -5.7 1.17e-01  -3.4 1.00e+00 1.00e+00h  1
       8  4.1448608e-03 4.91e+02 1.15e-02  -5.7 4.76e-01  -3.9 1.00e+00 1.00e+00H  1
       9  4.5136970e-03 1.24e+02 7.37e-03  -5.7 2.61e-01    -  1.00e+00 5.00e-01h  2
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.5149649e-03 8.67e+01 6.10e-03  -5.7 6.23e-01    -  1.00e+00 2.50e-01h  3
      11  4.2873775e-03 2.73e+02 4.39e-03  -5.7 2.73e+00    -  1.00e+00 6.25e-02h  5
      12  4.3214531e-03 2.34e+01 3.82e-04  -5.7 6.93e-02    -  1.00e+00 1.00e+00h  1
      13  4.3187301e-03 1.37e+00 1.49e-05  -5.7 1.47e-02    -  1.00e+00 1.00e+00h  1
      14  4.3185595e-03 1.68e-03 1.85e-08  -5.7 5.00e-04    -  1.00e+00 1.00e+00h  1
      15  4.3185591e-03 3.12e-05 1.37e-09  -8.6 7.31e-05    -  1.00e+00 1.00e+00h  1
      16  4.3185591e-03 7.28e-11 2.18e-14  -8.6 1.58e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 16
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.3185590827435208e-03    4.3185590827435208e-03
    Dual infeasibility......:   2.1774589460246032e-14    2.1774589460246032e-14
    Constraint violation....:   2.6908391096400935e-13    7.2759576141834246e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 27
    Number of objective gradient evaluations             = 17
    Number of equality constraint evaluations            = 27
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 17
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 16
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.016
    Total CPU secs in NLP function evaluations           =      0.026
    
    EXIT: Optimal Solution Found.
    2020-06-11 20:06:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
    
    Number of nonzeros in equality constraint Jacobian...:     3750
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:     2200
    
    Total number of variables............................:     1102
                         variables with only lower bounds:        0
                    variables with lower and upper bounds:      300
                         variables with only upper bounds:        0
    Total number of equality constraints.................:     1100
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  6.4861503e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  1.0484932e-01 1.40e+03 5.25e-01  -1.0 1.37e+04    -  9.89e-01 1.00e+00f  1
       2  5.7101353e-03 6.34e+01 2.53e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  5.9002864e-03 7.03e+01 5.80e-04  -1.7 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  5.8432784e-03 9.71e+00 7.76e-05  -1.7 3.74e-02    -  1.00e+00 1.00e+00h  1
       5  5.8178282e-03 5.49e-01 7.27e-05  -2.5 8.84e-03    -  1.00e+00 1.00e+00h  1
       6  5.7010488e-03 4.22e+01 3.99e-04  -3.8 7.47e-02    -  1.00e+00 1.00e+00h  1
       7  5.0827309e-03 1.01e+03 7.52e-03  -3.8 3.73e-01  -4.0 1.00e+00 1.00e+00h  1
       8  4.5070852e-03 4.17e+00 1.37e-03  -3.8 4.02e-02    -  1.00e+00 1.00e+00h  1
       9  4.4956806e-03 5.71e-01 2.17e-05  -3.8 9.76e-03    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.4958978e-03 3.49e-03 7.09e-08  -3.8 7.55e-04    -  1.00e+00 1.00e+00h  1
      11  4.4942463e-03 2.93e-01 8.67e-06  -5.7 7.06e-03    -  1.00e+00 1.00e+00h  1
      12  4.4942139e-03 1.48e-04 1.81e-09  -5.7 1.51e-04    -  1.00e+00 1.00e+00h  1
      13  4.4942137e-03 4.12e-05 1.26e-09  -8.6 8.38e-05    -  1.00e+00 1.00e+00h  1
      14  4.4942137e-03 4.37e-11 2.10e-14  -8.6 2.15e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 14
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.4942136563596267e-03    4.4942136563596267e-03
    Dual infeasibility......:   2.1006479927854700e-14    2.1006479927854700e-14
    Constraint violation....:   1.6145034657840562e-13    4.3655745685100555e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 15
    Number of objective gradient evaluations             = 15
    Number of equality constraint evaluations            = 15
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 15
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 14
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.019
    Total CPU secs in NLP function evaluations           =      0.011
    
    EXIT: Optimal Solution Found.
    2020-06-11 20:06:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
    
    Number of nonzeros in equality constraint Jacobian...:     3750
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:     2200
    
    Total number of variables............................:     1102
                         variables with only lower bounds:        0
                    variables with lower and upper bounds:      300
                         variables with only upper bounds:        0
    Total number of equality constraints.................:     1100
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  4.4929476e+01 3.03e+00 4.00e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  6.0080573e-02 1.32e+03 4.20e-01  -1.0 1.32e+04    -  9.94e-01 1.00e+00f  1
       2  5.1496710e-03 1.58e+01 1.61e-02  -1.7 4.39e+02    -  1.00e+00 1.00e+00h  1
       3  5.1778448e-03 5.91e-03 1.25e-05  -2.5 5.70e-01    -  1.00e+00 1.00e+00h  1
       4  4.6396972e-03 4.66e+02 3.73e-03  -3.8 2.55e-01    -  1.00e+00 1.00e+00h  1
       5  4.0942161e-03 7.24e+02 8.09e-03  -3.8 3.47e+00    -  1.00e+00 6.25e-02h  5
       6  4.0936827e-03 1.24e+02 8.60e-04  -3.8 1.33e-01    -  1.00e+00 1.00e+00h  1
       7  4.0690559e-03 2.28e+00 1.78e-05  -3.8 1.72e-02    -  1.00e+00 1.00e+00h  1
       8  4.0684746e-03 4.41e-04 2.72e-08  -3.8 3.35e-04    -  1.00e+00 1.00e+00h  1
       9  4.0668450e-03 2.91e-01 1.04e-05  -5.7 7.01e-03    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.0668195e-03 1.44e-04 1.60e-09  -5.7 1.48e-04    -  1.00e+00 1.00e+00h  1
      11  4.0668193e-03 4.09e-05 1.51e-09  -8.6 8.33e-05    -  1.00e+00 1.00e+00h  1
      12  4.0668193e-03 2.91e-11 2.12e-14  -8.6 2.12e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 12
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.0668192796870558e-03    4.0668192796870558e-03
    Dual infeasibility......:   2.1203332014342407e-14    2.1203332014342407e-14
    Constraint violation....:   5.6843418860808015e-14    2.9103830456733704e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 18
    Number of objective gradient evaluations             = 13
    Number of equality constraint evaluations            = 18
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 13
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 12
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.012
    Total CPU secs in NLP function evaluations           =      0.018
    
    EXIT: Optimal Solution Found.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
    
    Number of nonzeros in equality constraint Jacobian...:     3750
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:     2200
    
    Total number of variables............................:     1102
                         variables with only lower bounds:        0
                    variables with lower and upper bounds:      300
                         variables with only upper bounds:        0
    Total number of equality constraints.................:     1100
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  5.6092960e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  8.4804180e-02 1.41e+03 4.97e-01  -1.0 1.37e+04    -  9.88e-01 1.00e+00f  1
       2  5.7876676e-03 1.81e+01 1.75e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  5.8267302e-03 2.67e-03 1.80e-05  -2.5 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  5.1330206e-03 7.35e+02 6.89e-03  -3.8 3.20e-01    -  9.91e-01 1.00e+00h  1
       5  5.1078874e-03 1.62e+02 9.50e-03  -3.8 1.85e-01    -  1.00e+00 1.00e+00h  1
       6  4.5438636e-03 7.62e+01 2.31e-03  -3.8 1.14e-01    -  1.00e+00 1.00e+00h  1
       7  4.5078991e-03 1.16e+01 2.18e-04  -3.8 4.31e-02    -  1.00e+00 1.00e+00h  1
       8  4.5072371e-03 1.91e-01 2.86e-06  -3.8 5.41e-03    -  1.00e+00 1.00e+00h  1
       9  4.5056835e-03 2.38e-01 8.54e-06  -5.7 6.38e-03    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.5056573e-03 9.51e-05 1.30e-09  -5.7 1.21e-04    -  1.00e+00 1.00e+00h  1
      11  4.5056570e-03 3.46e-05 1.27e-09  -8.6 7.70e-05    -  1.00e+00 1.00e+00h  1
      12  4.5056570e-03 8.73e-11 1.92e-14  -8.6 1.77e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 12
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.5056570354069687e-03    4.5056570354069687e-03
    Dual infeasibility......:   1.9214097277204777e-14    1.9214097277204777e-14
    Constraint violation....:   3.2290069315681125e-13    8.7311491370201111e-11
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    2.5059035596800626e-09
    
    
    Number of objective function evaluations             = 13
    Number of objective gradient evaluations             = 13
    Number of equality constraint evaluations            = 13
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 13
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 12
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.012
    Total CPU secs in NLP function evaluations           =      0.016
    
    EXIT: Optimal Solution Found.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
    
    Number of nonzeros in equality constraint Jacobian...:     3750
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:     2200
    
    Total number of variables............................:     1102
                         variables with only lower bounds:        0
                    variables with lower and upper bounds:      300
                         variables with only upper bounds:        0
    Total number of equality constraints.................:     1100
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  5.1019924e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  7.4738160e-02 1.42e+03 4.69e-01  -1.0 1.37e+04    -  9.87e-01 1.00e+00f  1
       2  5.7666680e-03 8.04e+00 1.81e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  5.8027032e-03 2.70e-02 2.03e-05  -2.5 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  4.6690080e-03 2.04e+03 1.71e-02  -3.8 5.35e-01    -  9.51e-01 1.00e+00h  1
       5  4.8027577e-03 4.33e+02 2.35e-03  -3.8 2.43e-01    -  1.00e+00 1.00e+00h  1
       6  4.4617780e-03 2.76e+01 1.96e-03  -3.8 8.04e-02    -  1.00e+00 1.00e+00h  1
       7  4.4290485e-03 9.49e+00 2.54e-04  -3.8 3.99e-02    -  1.00e+00 1.00e+00h  1
       8  4.4289350e-03 2.21e-01 4.10e-06  -3.8 5.90e-03    -  1.00e+00 1.00e+00h  1
       9  4.4274484e-03 2.13e-01 8.79e-06  -5.7 6.03e-03    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.4274262e-03 7.38e-05 9.72e-10  -5.7 1.07e-04    -  1.00e+00 1.00e+00h  1
      11  4.4274260e-03 3.14e-05 1.32e-09  -8.6 7.33e-05    -  1.00e+00 1.00e+00h  1
      12  4.4274260e-03 2.91e-11 2.08e-14  -8.6 1.59e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 12
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.4274259707645022e-03    4.4274259707645022e-03
    Dual infeasibility......:   2.0811329193895519e-14    2.0811329193895519e-14
    Constraint violation....:   1.0763356438560375e-13    2.9103830456733704e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 13
    Number of objective gradient evaluations             = 13
    Number of equality constraint evaluations            = 13
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 13
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 12
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.013
    Total CPU secs in NLP function evaluations           =      0.015
    
    EXIT: Optimal Solution Found.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:06:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
    
    Number of nonzeros in equality constraint Jacobian...:     3750
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:     2200
    
    Total number of variables............................:     1102
                         variables with only lower bounds:        0
                    variables with lower and upper bounds:      300
                         variables with only upper bounds:        0
    Total number of equality constraints.................:     1100
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  5.6333574e+01 3.03e+00 4.00e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  7.7897939e-02 1.31e+03 4.31e-01  -1.0 1.32e+04    -  9.94e-01 1.00e+00f  1
       2  6.4081717e-03 2.53e+01 1.76e-02  -1.7 4.39e+02    -  1.00e+00 1.00e+00h  1
       3  6.4059215e-03 2.61e-01 6.09e-05  -2.5 5.70e-01    -  1.00e+00 1.00e+00h  1
       4  5.2816972e-03 1.71e+03 1.41e-02  -3.8 4.86e-01    -  9.60e-01 1.00e+00h  1
       5  5.0797338e-03 2.80e+02 2.91e-03  -3.8 1.89e-01    -  1.00e+00 1.00e+00h  1
       6  4.8101089e-03 1.10e+02 3.72e-03  -3.8 1.43e-01    -  1.00e+00 1.00e+00h  1
       7  4.7257805e-03 2.52e+01 5.56e-04  -3.8 6.43e-02    -  1.00e+00 1.00e+00h  1
       8  4.7225342e-03 9.15e-01 1.48e-05  -3.8 1.19e-02    -  1.00e+00 1.00e+00h  1
       9  4.7224743e-03 7.59e-04 9.78e-09  -3.8 3.37e-04    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.7210793e-03 1.84e-01 7.52e-06  -5.7 5.63e-03    -  1.00e+00 1.00e+00h  1
      11  4.7210570e-03 5.30e-05 7.03e-10  -5.7 9.07e-05    -  1.00e+00 1.00e+00h  1
      12  4.7210568e-03 2.65e-05 1.11e-09  -8.6 6.76e-05    -  1.00e+00 1.00e+00h  1
      13  4.7210568e-03 4.37e-11 2.14e-14  -8.6 1.32e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.7210567672512803e-03    4.7210567672512803e-03
    Dual infeasibility......:   2.1435085580569746e-14    2.1435085580569746e-14
    Constraint violation....:   1.6145034657840562e-13    4.3655745685100555e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 14
    Number of objective gradient evaluations             = 14
    Number of equality constraint evaluations            = 14
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 14
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 13
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.013
    Total CPU secs in NLP function evaluations           =      0.023
    
    EXIT: Optimal Solution Found.
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
    
    Number of nonzeros in equality constraint Jacobian...:     3750
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:     2200
    
    Total number of variables............................:     1102
                         variables with only lower bounds:        0
                    variables with lower and upper bounds:      300
                         variables with only upper bounds:        0
    Total number of equality constraints.................:     1100
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  5.9652492e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  8.4002981e-02 1.41e+03 4.94e-01  -1.0 1.37e+04    -  9.88e-01 1.00e+00f  1
       2  6.3664448e-03 1.57e+01 1.75e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  6.3896178e-03 1.19e-01 3.78e-05  -2.5 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  5.1259730e-03 2.47e+03 2.16e-02  -3.8 5.87e-01    -  9.42e-01 1.00e+00h  1
       5  5.3794858e-03 5.57e+02 2.77e-03  -3.8 2.77e-01    -  1.00e+00 1.00e+00h  1
       6  4.9670642e-03 8.91e+01 4.72e-03  -3.8 1.38e-01    -  1.00e+00 1.00e+00h  1
       7  4.7910404e-03 3.43e+01 9.05e-04  -3.8 7.62e-02    -  1.00e+00 1.00e+00h  1
       8  4.7840864e-03 2.32e+00 4.18e-05  -3.8 1.91e-02    -  1.00e+00 1.00e+00h  1
       9  4.7839573e-03 6.29e-03 9.21e-08  -3.8 9.79e-04    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.7825394e-03 1.92e-01 7.54e-06  -5.7 5.74e-03    -  1.00e+00 1.00e+00h  1
      11  4.7825152e-03 5.79e-05 8.08e-10  -5.7 9.48e-05    -  1.00e+00 1.00e+00h  1
      12  4.7825150e-03 2.75e-05 1.11e-09  -8.6 6.89e-05    -  1.00e+00 1.00e+00h  1
      13  4.7825149e-03 1.16e-10 2.17e-14  -8.6 1.38e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.7825149479292905e-03    4.7825149479292905e-03
    Dual infeasibility......:   2.1665108228351927e-14    2.1665108228351927e-14
    Constraint violation....:   4.3053425754241500e-13    1.1641532182693481e-10
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    2.5059035596800626e-09
    
    
    Number of objective function evaluations             = 14
    Number of objective gradient evaluations             = 14
    Number of equality constraint evaluations            = 14
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 14
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 13
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.018
    Total CPU secs in NLP function evaluations           =      0.012
    
    EXIT: Optimal Solution Found.
    2020-06-11 20:07:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
    
    Number of nonzeros in equality constraint Jacobian...:     3750
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:     2200
    
    Total number of variables............................:     1102
                         variables with only lower bounds:        0
                    variables with lower and upper bounds:      300
                         variables with only upper bounds:        0
    Total number of equality constraints.................:     1100
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  5.3514250e+01 3.03e+00 4.00e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  6.9507794e-02 1.32e+03 4.21e-01  -1.0 1.32e+04    -  9.94e-01 1.00e+00f  1
       2  6.2046627e-03 1.53e+01 1.60e-02  -1.7 4.39e+02    -  1.00e+00 1.00e+00h  1
       3  6.2077841e-03 1.78e-01 5.03e-05  -2.5 5.70e-01    -  1.00e+00 1.00e+00h  1
       4  4.8825224e-03 2.70e+03 2.16e-02  -3.8 6.13e-01    -  9.38e-01 1.00e+00h  1
       5  5.2286948e-03 6.29e+02 3.24e-03  -3.8 2.95e-01    -  1.00e+00 1.00e+00h  1
       6  4.7263694e-03 4.28e+01 3.52e-03  -3.8 1.03e-01    -  1.00e+00 1.00e+00h  1
       7  4.6117380e-03 2.11e+01 6.30e-04  -3.8 6.01e-02    -  1.00e+00 1.00e+00h  1
       8  4.6092642e-03 1.11e+00 2.22e-05  -3.8 1.33e-02    -  1.00e+00 1.00e+00h  1
       9  4.6092331e-03 1.60e-03 2.48e-08  -3.8 4.95e-04    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.6078258e-03 1.90e-01 7.95e-06  -5.7 5.70e-03    -  1.00e+00 1.00e+00h  1
      11  4.6078039e-03 5.64e-05 7.64e-10  -5.7 9.34e-05    -  1.00e+00 1.00e+00h  1
      12  4.6078037e-03 2.72e-05 1.17e-09  -8.6 6.84e-05    -  1.00e+00 1.00e+00h  1
      13  4.6078037e-03 1.60e-10 2.05e-14  -8.6 1.36e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.6078036578471211e-03    4.6078036578471211e-03
    Dual infeasibility......:   2.0522352825509770e-14    2.0522352825509770e-14
    Constraint violation....:   5.9198460412082065e-13    1.6007106751203537e-10
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 14
    Number of objective gradient evaluations             = 14
    Number of equality constraint evaluations            = 14
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 14
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 13
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.014
    Total CPU secs in NLP function evaluations           =      0.022
    
    EXIT: Optimal Solution Found.
    2020-06-11 20:07:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-06-11 20:07:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-06-11 20:07:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-06-11 20:07:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-06-11 20:07:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
    
    Number of nonzeros in equality constraint Jacobian...:     3750
    Number of nonzeros in inequality constraint Jacobian.:        0
    Number of nonzeros in Lagrangian Hessian.............:     2200
    
    Total number of variables............................:     1102
                         variables with only lower bounds:        0
                    variables with lower and upper bounds:      300
                         variables with only upper bounds:        0
    Total number of equality constraints.................:     1100
    Total number of inequality constraints...............:        0
            inequality constraints with only lower bounds:        0
       inequality constraints with lower and upper bounds:        0
            inequality constraints with only upper bounds:        0
    
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
       0  3.1960790e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  3.6328141e-02 1.42e+03 4.58e-01  -1.0 1.37e+04    -  9.87e-01 1.00e+00f  1
       2  3.7177243e-03 4.36e+00 1.91e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  3.7105275e-03 1.33e+00 1.54e-04  -3.8 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  3.3360154e-03 4.37e+02 3.73e-03  -3.8 2.50e-01  -4.0 1.00e+00 1.00e+00h  1
       5  3.2135952e-03 6.67e+00 3.12e-04  -3.8 3.35e-02    -  1.00e+00 1.00e+00h  1
       6  3.2097108e-03 1.90e-01 3.50e-06  -3.8 5.76e-03    -  1.00e+00 1.00e+00h  1
       7  3.2076689e-03 5.84e-01 1.67e-05  -5.7 9.81e-03    -  1.00e+00 1.00e+00h  1
       8  3.2076320e-03 6.66e-04 6.43e-09  -5.7 3.16e-04    -  1.00e+00 1.00e+00h  1
       9  3.2076317e-03 8.04e-05 2.41e-09  -8.6 1.15e-04    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  3.2076317e-03 2.91e-11 2.21e-14  -8.6 4.46e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 10
    
                                       (scaled)                 (unscaled)
    Objective...............:   3.2076316541892677e-03    3.2076316541892677e-03
    Dual infeasibility......:   2.2139808557097437e-14    2.2139808557097437e-14
    Constraint violation....:   1.0763356438560375e-13    2.9103830456733704e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 11
    Number of objective gradient evaluations             = 11
    Number of equality constraint evaluations            = 11
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 11
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 10
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.010
    Total CPU secs in NLP function evaluations           =      0.015
    
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
          <td>0.476793</td>
          <td>-0.408376</td>
        </tr>
        <tr>
          <th>1</th>
          <td>0.466116</td>
          <td>-0.399796</td>
        </tr>
        <tr>
          <th>2</th>
          <td>0.458028</td>
          <td>-0.394253</td>
        </tr>
        <tr>
          <th>3</th>
          <td>0.449050</td>
          <td>-0.386724</td>
        </tr>
        <tr>
          <th>4</th>
          <td>0.465797</td>
          <td>-0.399956</td>
        </tr>
        <tr>
          <th>5</th>
          <td>0.468683</td>
          <td>-0.401897</td>
        </tr>
        <tr>
          <th>6</th>
          <td>0.482269</td>
          <td>-0.412374</td>
        </tr>
        <tr>
          <th>7</th>
          <td>0.481936</td>
          <td>-0.412269</td>
        </tr>
        <tr>
          <th>8</th>
          <td>0.478817</td>
          <td>-0.409642</td>
        </tr>
        <tr>
          <th>9</th>
          <td>0.398945</td>
          <td>-0.347210</td>
        </tr>
      </tbody>
    </table>
    </div>

