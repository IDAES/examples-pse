Parameter Estimation Using the NRTL State Block
-----------------------------------------------

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
    
        
        # Fix the state variables on the state block
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

    2020-07-21 19:09:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:09:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
    Total CPU secs in NLP function evaluations           =      0.022
    
    EXIT: Optimal Solution Found.
    

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

    2020-07-21 19:10:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  5.4705532e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  7.7359849e-02 1.42e+03 4.73e-01  -1.0 1.37e+04    -  9.87e-01 1.00e+00f  1
       2  6.1740435e-03 9.79e+00 1.78e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  6.2018817e-03 1.39e-01 4.37e-05  -2.5 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  4.7446564e-03 4.14e+03 3.34e-02  -3.8 7.62e-01    -  9.13e-01 1.00e+00h  1
       5  5.6002884e-03 1.07e+03 7.79e-03  -3.8 3.96e-01    -  1.00e+00 1.00e+00h  1
       6  4.6585713e-03 2.15e+01 1.55e-03  -3.8 1.28e-01    -  1.00e+00 1.00e+00h  1
       7  4.6373143e-03 1.90e-01 2.99e-06  -3.8 3.43e-02    -  1.00e+00 1.00e+00h  1
       8  4.6358151e-03 1.94e-01 8.37e-06  -5.7 5.77e-03    -  1.00e+00 1.00e+00h  1
       9  4.6357935e-03 6.20e-05 8.87e-10  -5.7 9.81e-05    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.6357933e-03 2.74e-05 1.20e-09  -8.6 6.86e-05    -  1.00e+00 1.00e+00h  1
      11  4.6357933e-03 1.46e-11 2.19e-14  -8.6 1.37e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 11
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.6357933119725320e-03    4.6357933119725320e-03
    Dual infeasibility......:   2.1852115100346825e-14    2.1852115100346825e-14
    Constraint violation....:   5.6843418860808015e-14    1.4551915228366852e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 12
    Number of objective gradient evaluations             = 12
    Number of equality constraint evaluations            = 12
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 12
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 11
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.031
    Total CPU secs in NLP function evaluations           =      0.044
    
    EXIT: Optimal Solution Found.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  4.4603578e+01 2.79e+00 3.68e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  5.2855825e-02 1.12e+03 3.42e-01  -1.0 1.22e+04    -  1.00e+00 1.00e+00f  1
       2  5.2496916e-03 4.83e+01 2.37e-02  -1.7 3.73e+02    -  1.00e+00 1.00e+00h  1
       3  5.2217163e-03 1.63e-02 4.39e-05  -2.5 4.14e-01    -  1.00e+00 1.00e+00h  1
       4  4.9126666e-03 1.86e+02 1.62e-03  -3.8 1.60e-01    -  1.00e+00 1.00e+00h  1
       5  4.7074423e-03 6.25e-02 6.66e-05  -3.8 6.59e-03  -2.0 1.00e+00 1.00e+00h  1
       6  4.6697861e-03 5.16e-01 3.07e-05  -5.7 8.89e-03  -2.5 1.00e+00 1.00e+00h  1
       7  4.5563520e-03 5.24e+00 8.19e-05  -5.7 2.84e-02  -3.0 1.00e+00 1.00e+00h  1
       8  4.2393738e-03 5.23e+01 8.07e-04  -5.7 9.03e-02  -3.4 1.00e+00 1.00e+00h  1
       9  4.1362828e-03 1.08e+02 2.74e-03  -5.7 1.34e-01    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.0665198e-03 2.02e+01 3.84e-04  -5.7 5.68e-02    -  1.00e+00 1.00e+00h  1
      11  4.0628480e-03 6.33e-01 8.63e-06  -5.7 9.85e-03    -  1.00e+00 1.00e+00h  1
      12  4.0627772e-03 3.93e-04 4.66e-09  -5.7 2.42e-04    -  1.00e+00 1.00e+00h  1
      13  4.0627769e-03 4.20e-05 1.41e-09  -8.6 8.43e-05    -  1.00e+00 1.00e+00h  1
      14  4.0627769e-03 1.46e-11 2.15e-14  -8.6 2.18e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 14
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.0627769309870060e-03    4.0627769309870060e-03
    Dual infeasibility......:   2.1483925725301556e-14    2.1483925725301556e-14
    Constraint violation....:   5.6843418860808015e-14    1.4551915228366852e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 15
    Number of objective gradient evaluations             = 15
    Number of equality constraint evaluations            = 15
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 15
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 14
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.016
    Total CPU secs in NLP function evaluations           =      0.026
    
    EXIT: Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  5.2515380e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  7.4474581e-02 1.41e+03 4.90e-01  -1.0 1.37e+04    -  9.88e-01 1.00e+00f  1
       2  5.5802755e-03 1.25e+01 1.74e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  5.6212258e-03 1.05e-03 1.54e-05  -2.5 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  4.9119984e-03 7.78e+02 7.24e-03  -3.8 3.30e-01    -  9.89e-01 1.00e+00h  1
       5  4.4223108e-03 3.20e-01 1.15e-03  -3.8 3.69e-02    -  1.00e+00 1.00e+00h  1
       6  4.3881659e-03 4.19e+00 1.27e-04  -5.7 2.69e-02    -  1.00e+00 1.00e+00h  1
       7  4.3873743e-03 8.81e-02 1.83e-06  -5.7 3.77e-03    -  1.00e+00 1.00e+00h  1
       8  4.3873637e-03 1.43e-05 2.21e-10  -5.7 4.70e-05    -  1.00e+00 1.00e+00h  1
       9  4.3873635e-03 3.66e-05 1.35e-09  -8.6 7.90e-05    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.3873635e-03 7.28e-11 2.18e-14  -8.6 1.88e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 10
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.3873635039169594e-03    4.3873635039169594e-03
    Dual infeasibility......:   2.1766601399404484e-14    2.1766601399404484e-14
    Constraint violation....:   2.6908391096400935e-13    7.2759576141834246e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 11
    Number of objective gradient evaluations             = 11
    Number of equality constraint evaluations            = 11
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 11
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 10
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.019
    Total CPU secs in NLP function evaluations           =      0.022
    
    EXIT: Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  5.8585149e+01 2.91e+00 3.84e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  7.8047361e-02 1.21e+03 4.06e-01  -1.0 1.27e+04    -  1.00e+00 1.00e+00f  1
       2  6.7458067e-03 2.76e+02 6.36e-02  -1.7 4.05e+02    -  1.00e+00 1.00e+00h  1
       3  6.4144967e-03 5.99e+01 6.38e-04  -1.7 4.87e-01    -  1.00e+00 1.00e+00h  1
       4  6.4095596e-03 1.84e+01 9.80e-05  -1.7 4.93e-02    -  1.00e+00 1.00e+00h  1
       5  6.3924435e-03 8.30e-01 9.22e-05  -2.5 1.10e-02    -  1.00e+00 1.00e+00h  1
       6  6.2119090e-03 6.81e+01 6.77e-04  -3.8 9.50e-02    -  1.00e+00 1.00e+00h  1
       7  6.1675357e-03 1.04e-03 3.21e-05  -3.8 3.20e-03  -2.0 1.00e+00 1.00e+00h  1
       8  6.1507621e-03 2.33e-01 1.96e-05  -5.7 5.75e-03  -2.5 1.00e+00 1.00e+00h  1
       9  6.0949532e-03 2.59e+00 3.17e-05  -5.7 1.92e-02  -3.0 1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  5.8302245e-03 5.09e+01 5.35e-04  -5.7 8.55e-02  -3.4 1.00e+00 1.00e+00h  1
      11  5.5758105e-03 9.90e+00 1.87e-04  -5.7 4.01e-02  -3.0 1.00e+00 1.00e+00h  1
      12  4.9211986e-03 1.68e+02 2.65e-03  -5.7 1.61e-01  -3.5 1.00e+00 1.00e+00h  1
      13  4.7276539e-03 3.53e+01 1.61e-03  -5.7 8.15e-02    -  1.00e+00 1.00e+00h  1
      14  4.6902143e-03 7.19e+00 1.68e-04  -5.7 3.46e-02    -  1.00e+00 1.00e+00h  1
      15  4.6887638e-03 1.12e-01 1.92e-06  -5.7 4.19e-03    -  1.00e+00 1.00e+00h  1
      16  4.6887467e-03 1.42e-05 2.09e-10  -5.7 4.64e-05    -  1.00e+00 1.00e+00h  1
      17  4.6887465e-03 3.04e-05 1.07e-09  -8.6 7.23e-05    -  1.00e+00 1.00e+00h  1
      18  4.6887465e-03 1.46e-10 2.18e-14  -8.6 1.54e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 18
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.6887464777415136e-03    4.6887464777415136e-03
    Dual infeasibility......:   2.1817085217360705e-14    2.1817085217360705e-14
    Constraint violation....:   5.3816782192801870e-13    1.4551915228366849e-10
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 19
    Number of objective gradient evaluations             = 19
    Number of equality constraint evaluations            = 19
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 19
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 18
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.032
    Total CPU secs in NLP function evaluations           =      0.022
    
    EXIT: Optimal Solution Found.
    2020-07-21 19:10:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  4.4569574e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  4.9507652e-02 1.42e+03 4.76e-01  -1.0 1.37e+04    -  9.88e-01 1.00e+00f  1
       2  5.0326281e-03 6.17e+00 1.76e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  5.0028363e-03 3.01e+00 2.14e-04  -3.8 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  4.1001740e-03 2.73e+03 2.24e-02  -3.8 1.25e+00  -4.0 1.00e+00 5.00e-01h  2
       5  4.4716487e-03 6.88e+02 4.58e-03  -3.8 3.15e-01    -  1.00e+00 1.00e+00h  1
       6  4.0703386e-03 1.71e+01 5.72e-04  -3.8 5.35e-02    -  1.00e+00 1.00e+00h  1
       7  4.0593308e-03 4.28e-01 7.19e-06  -3.8 8.91e-03    -  1.00e+00 1.00e+00h  1
       8  4.0593225e-03 3.99e-04 3.53e-09  -3.8 2.48e-04    -  1.00e+00 1.00e+00h  1
       9  4.0576498e-03 3.14e-01 1.09e-05  -5.7 7.28e-03    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.0576184e-03 1.71e-04 1.91e-09  -5.7 1.62e-04    -  1.00e+00 1.00e+00h  1
      11  4.0576182e-03 4.40e-05 1.58e-09  -8.6 8.63e-05    -  1.00e+00 1.00e+00h  1
      12  4.0576182e-03 5.82e-11 2.21e-14  -8.6 2.31e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 12
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.0576181698097537e-03    4.0576181698097537e-03
    Dual infeasibility......:   2.2140424874482600e-14    2.2140424874482600e-14
    Constraint violation....:   2.1526712877120750e-13    5.8207660913467407e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 15
    Number of objective gradient evaluations             = 13
    Number of equality constraint evaluations            = 15
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 13
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 12
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.018
    Total CPU secs in NLP function evaluations           =      0.019
    
    EXIT: Optimal Solution Found.
    2020-07-21 19:10:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:10:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  5.5968915e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  8.3499230e-02 1.42e+03 4.65e-01  -1.0 1.37e+04    -  9.87e-01 1.00e+00f  1
       2  6.3876746e-03 6.81e+00 1.84e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  6.3063636e-03 8.99e+00 4.17e-04  -3.8 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  6.2790242e-03 4.98e-02 2.96e-05  -3.8 2.95e-03  -2.0 1.00e+00 1.00e+00h  1
       5  6.2404164e-03 5.61e-01 3.11e-05  -5.7 9.06e-03  -2.5 1.00e+00 1.00e+00h  1
       6  6.1027925e-03 6.74e+00 8.05e-05  -5.7 3.15e-02  -3.0 1.00e+00 1.00e+00h  1
       7  5.4098942e-03 1.69e+02 1.93e-03  -5.7 1.58e-01  -3.4 1.00e+00 1.00e+00h  1
       8  5.2926211e-03 8.23e+02 1.84e-02  -5.7 1.35e+00    -  1.00e+00 2.50e-01h  3
       9  4.7946246e-03 1.94e+02 3.65e-03  -5.7 1.74e-01    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.7244615e-03 1.60e+01 1.67e-04  -5.7 4.82e-02    -  1.00e+00 1.00e+00h  1
      11  4.7213497e-03 5.87e-02 2.99e-07  -5.7 2.78e-03    -  1.00e+00 1.00e+00h  1
      12  4.7213386e-03 4.70e-08 2.58e-12  -5.7 3.20e-06    -  1.00e+00 1.00e+00h  1
      13  4.7213384e-03 2.47e-05 1.15e-09  -8.6 6.53e-05    -  1.00e+00 1.00e+00h  1
      14  4.7213384e-03 4.37e-11 2.19e-14  -8.6 1.22e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 14
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.7213384447481084e-03    4.7213384447481084e-03
    Dual infeasibility......:   2.1894078466648466e-14    2.1894078466648466e-14
    Constraint violation....:   1.6145034657840562e-13    4.3655745685100555e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 18
    Number of objective gradient evaluations             = 15
    Number of equality constraint evaluations            = 18
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 15
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 14
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.026
    Total CPU secs in NLP function evaluations           =      0.038
    
    EXIT: Optimal Solution Found.
    2020-07-21 19:11:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  5.5100702e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  7.5547650e-02 1.42e+03 4.66e-01  -1.0 1.37e+04    -  9.87e-01 1.00e+00f  1
       2  6.3100561e-03 7.56e+00 1.83e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  6.2321627e-03 7.94e+00 3.86e-04  -3.8 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  6.2077541e-03 4.63e-02 2.83e-05  -3.8 2.82e-03  -2.0 1.00e+00 1.00e+00h  1
       5  6.1720315e-03 5.19e-01 2.98e-05  -5.7 8.70e-03  -2.5 1.00e+00 1.00e+00h  1
       6  6.0452302e-03 6.17e+00 7.35e-05  -5.7 3.01e-02  -3.0 1.00e+00 1.00e+00h  1
       7  5.4084876e-03 1.51e+02 1.70e-03  -5.7 1.50e-01  -3.4 1.00e+00 1.00e+00h  1
       8  4.8508448e-03 6.71e+02 1.23e-02  -5.7 1.69e+01    -  3.17e-01 1.71e-02h  5
       9  4.7167617e-03 1.33e+02 1.51e-03  -5.7 1.42e-01    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.6831665e-03 4.53e+00 2.31e-05  -5.7 2.47e-02    -  1.00e+00 1.00e+00h  1
      11  4.6822686e-03 1.74e-04 2.03e-08  -5.7 2.15e-04    -  1.00e+00 1.00e+00h  1
      12  4.6822682e-03 2.56e-05 1.17e-09  -8.6 6.64e-05    -  1.00e+00 1.00e+00h  1
      13  4.6822682e-03 4.37e-11 2.14e-14  -8.6 1.26e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.6822681752350477e-03    4.6822681752350477e-03
    Dual infeasibility......:   2.1426204137369681e-14    2.1426204137369681e-14
    Constraint violation....:   1.6145034657840562e-13    4.3655745685100555e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 19
    Number of objective gradient evaluations             = 14
    Number of equality constraint evaluations            = 19
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 14
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 13
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.024
    Total CPU secs in NLP function evaluations           =      0.018
    
    EXIT: Optimal Solution Found.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  5.7104224e+01 2.56e+00 3.36e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  6.6612962e-02 9.44e+02 2.78e-01  -1.0 1.12e+04    -  1.00e+00 1.00e+00f  1
       2  6.8344329e-03 3.44e+01 1.99e-02  -1.7 3.12e+02    -  1.00e+00 1.00e+00h  1
       3  6.7767919e-03 6.39e-01 8.54e-05  -2.5 2.93e-01    -  1.00e+00 1.00e+00h  1
       4  5.4421417e-03 2.25e+03 1.74e-02  -3.8 5.56e-01    -  9.48e-01 1.00e+00h  1
       5  5.4827743e-03 4.59e+02 2.79e-03  -3.8 2.47e-01    -  1.00e+00 1.00e+00h  1
       6  5.4380978e-03 3.22e+02 1.03e-02  -3.8 2.44e-01    -  1.00e+00 1.00e+00h  1
       7  4.8905320e-03 9.97e+01 2.36e-03  -3.8 1.29e-01    -  1.00e+00 1.00e+00h  1
       8  4.8469130e-03 1.11e+01 1.80e-04  -3.8 4.17e-02    -  1.00e+00 1.00e+00h  1
       9  4.8456412e-03 1.09e-01 1.40e-06  -3.8 4.05e-03    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.8443091e-03 1.56e-01 6.04e-06  -5.7 5.20e-03    -  1.00e+00 1.00e+00h  1
      11  4.8442880e-03 3.75e-05 4.84e-10  -5.7 7.64e-05    -  1.00e+00 1.00e+00h  1
      12  4.8442878e-03 2.29e-05 9.01e-10  -8.6 6.30e-05    -  1.00e+00 1.00e+00h  1
      13  4.8442878e-03 1.46e-11 2.18e-14  -8.6 1.13e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.8442878123546871e-03    4.8442878123546871e-03
    Dual infeasibility......:   2.1781404500396822e-14    2.1781404500396822e-14
    Constraint violation....:   5.6843418860808015e-14    1.4551915228366852e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 14
    Number of objective gradient evaluations             = 14
    Number of equality constraint evaluations            = 14
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 14
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 13
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.019
    Total CPU secs in NLP function evaluations           =      0.032
    
    EXIT: Optimal Solution Found.
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  6.5161124e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  1.0065104e-01 1.40e+03 5.16e-01  -1.0 1.37e+04    -  9.88e-01 1.00e+00f  1
       2  6.2330428e-03 4.23e+01 2.32e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  6.2511330e-03 9.00e-02 3.72e-05  -2.5 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  5.7502888e-03 4.04e+02 3.69e-03  -3.8 2.36e-01    -  1.00e+00 1.00e+00h  1
       5  5.2191588e-03 3.45e-01 2.21e-04  -3.8 1.23e-02  -2.0 1.00e+00 1.00e+00h  1
       6  4.9225004e-03 4.67e+02 8.63e-03  -5.7 5.45e-01    -  9.18e-01 5.00e-01h  2
       7  4.7732896e-03 8.75e+01 1.11e-03  -5.7 1.16e-01    -  1.00e+00 1.00e+00h  1
       8  4.7492582e-03 3.70e+00 3.77e-05  -5.7 2.31e-02    -  1.00e+00 1.00e+00h  1
       9  4.7485993e-03 3.60e-03 2.57e-08  -5.7 6.93e-04    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.7485985e-03 3.15e-05 1.12e-09  -8.6 7.37e-05    -  1.00e+00 1.00e+00h  1
      11  4.7485985e-03 2.91e-11 2.21e-14  -8.6 1.63e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 11
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.7485984592461322e-03    4.7485984592461322e-03
    Dual infeasibility......:   2.2093717116170349e-14    2.2093717116170349e-14
    Constraint violation....:   1.0763356438560375e-13    2.9103830456733704e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 14
    Number of objective gradient evaluations             = 12
    Number of equality constraint evaluations            = 14
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 12
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 11
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.018
    Total CPU secs in NLP function evaluations           =      0.019
    
    EXIT: Optimal Solution Found.
    2020-07-21 19:11:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 19:11:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 19:11:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 19:11:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 19:11:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 19:11:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 19:11:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 19:11:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  6.3680621e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  8.6348135e-02 1.41e+03 4.85e-01  -1.0 1.37e+04    -  9.88e-01 1.00e+00f  1
       2  7.0910810e-03 1.39e+01 1.73e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  7.0928807e-03 6.59e-01 8.97e-05  -2.5 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  5.5329371e-03 1.50e+04 1.13e-01  -3.8 2.95e+00    -  6.58e-01 5.00e-01h  2
       5  5.6298669e-03 1.45e+04 1.10e-01  -3.8 1.71e+00    -  8.44e-01 3.12e-02h  6
       6  5.6415195e-03 1.45e+04 1.09e-01  -3.8 2.06e+00    -  1.00e+00 3.91e-03h  9
       7  5.6444433e-03 1.45e+04 1.09e-01  -3.8 2.10e+00    -  1.00e+00 9.77e-04h 11
       8  5.6451739e-03 1.45e+04 1.09e-01  -3.8 2.11e+00    -  1.00e+00 2.44e-04h 13
       9  5.6455393e-03 1.45e+04 1.09e-01  -3.8 2.11e+00    -  1.00e+00 1.22e-04h 14
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  1.7707006e-02 4.08e+03 1.19e-01  -3.8 2.11e+00    -  1.00e+00 1.00e+00h  1
      11  1.8350416e-02 4.08e+03 2.09e-01  -3.8 1.58e+03  -4.0 2.40e-01 2.42e-03h  8
      12  1.8988464e-02 4.07e+03 4.13e+01  -3.8 1.58e+03  -2.7 1.00e+00 1.04e-02h  7
      13  5.8811173e+00 1.96e+02 4.63e+04  -3.8 1.56e+03  -2.2 1.00e+00 1.00e+00h  1
      14  5.9573610e+00 8.09e+01 1.58e+04  -3.8 1.58e+03    -  1.00e+00 1.00e+00h  1
      15  4.5844658e-01 4.30e+02 8.25e+01  -3.8 2.57e+03    -  1.00e+00 1.00e+00f  1
      16  2.2281190e+00 4.23e+02 8.07e+01  -3.8 2.19e+03  -0.9 7.47e-01 2.14e-02h  6
      17  3.9771035e+00 4.18e+02 7.95e+01  -3.8 2.15e+03  -0.5 1.00e+00 1.56e-02h  7
      18  5.0328515e+00 4.15e+02 7.70e+01  -3.8 2.12e+03  -1.0 1.00e+00 3.12e-02h  6
      19  5.5892434e+00 4.12e+02 7.64e+01  -3.8 2.05e+03  -0.5 1.00e+00 7.81e-03h  8
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      20  5.7706725e+00 4.10e+02 7.58e+01  -3.8 2.04e+03  -1.0 1.00e+00 7.81e-03h  8
      21  5.8283128e+00 4.10e+02 7.57e+01  -3.8 2.02e+03  -0.6 1.00e+00 9.77e-04h 11
      22  5.8662407e+00 4.09e+02 7.55e+01  -3.8 2.02e+03  -1.1 1.00e+00 1.95e-03h 10
      23  1.2947762e+02 3.08e+02 6.12e+06  -3.8 2.02e+03  -0.6 1.00e+00 1.00e+00H  1
      24  1.2816033e+02 3.29e+01 6.31e+05  -3.8 1.41e+01   0.7 1.00e+00 1.00e+00H  1
      25  1.1531196e+02 9.98e+01 3.89e+05  -3.8 1.36e+03    -  1.00e+00 2.50e-01f  3
      26  1.1267720e+02 1.25e+02 3.20e+05  -3.8 6.69e+01    -  1.00e+00 1.25e-01f  4
      27  1.1021718e+02 1.49e+02 2.61e+05  -3.8 5.36e+01    -  1.00e+00 1.25e-01f  4
      28  1.0785395e+02 1.74e+02 2.13e+05  -3.8 4.41e+01    -  1.00e+00 1.25e-01f  4
      29  8.9975243e+01 3.21e+03 8.80e+05  -3.8 3.69e+01    -  1.00e+00 1.00e+00f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      30  9.0398618e+01 9.26e+01 7.81e+03  -3.8 1.69e-01   0.2 1.00e+00 1.00e+00h  1
      31  8.3022276e+01 6.04e+02 4.14e+02  -3.8 3.09e-01  -0.3 1.00e+00 1.00e+00f  1
      32  2.4682408e+01 4.68e+04 1.07e+04  -3.8 3.95e+01  -0.7 1.46e-01 1.07e-01f  1
      33  2.4687702e+01 4.68e+04 1.07e+04  -3.8 1.84e+03  -1.2 9.56e-02 1.33e-04h  1
      34  6.3643096e+01 4.62e+04 1.05e+04  -3.8 6.21e+05    -  2.57e-04 1.21e-02h  4
      35r 6.3643096e+01 4.62e+04 9.99e+02   2.2 0.00e+00  -1.7 0.00e+00 3.57e-07R  5
      36r 7.1691771e+02 9.94e+04 1.08e+09   2.2 1.67e+06    -  7.50e-06 1.05e-05f  1
      37r 7.0906982e+02 9.93e+04 1.08e+09   2.2 4.35e+01   6.0 4.63e-03 2.56e-03f  1
      38r 7.0688124e+02 9.93e+04 1.08e+09   2.2 5.41e+01   5.5 6.53e-03 7.52e-04f  1
      39r 6.8973933e+02 9.89e+04 1.07e+09   2.2 6.44e+01   5.0 6.84e-03 6.09e-03f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      40r 6.8763316e+02 9.88e+04 1.07e+09   2.2 4.01e+01   5.5 1.82e-02 8.51e-04h  1
      41r 6.5976514e+02 9.81e+04 1.05e+09   2.2 4.60e+01   5.0 1.41e-02 1.21e-02f  1
      42r 6.5121089e+02 9.78e+04 1.05e+09   2.2 2.60e+01   5.4 5.04e-02 4.84e-03h  1
      43r 5.9558572e+02 9.58e+04 1.02e+09   2.2 2.54e+01   4.9 6.04e-02 2.99e-02h  1
      44r 5.7326269e+02 9.43e+04 9.91e+08   2.2 1.38e+01   5.4 2.58e-01 2.69e-02h  1
      45r 5.5814538e+02 9.30e+04 9.67e+08   2.2 1.05e+01   5.8 4.43e-01 2.43e-02h  1
      46r 5.0989430e+02 8.33e+04 8.23e+08   2.2 8.74e+00   5.3 1.84e-01 1.50e-01h  1
      47r 4.8192146e+02 7.89e+04 7.34e+08   2.2 3.84e+00   5.7 2.72e-01 1.08e-01h  1
      48r 4.2342051e+02 6.70e+04 5.23e+08   2.2 3.94e+00   5.3 1.52e-01 2.86e-01h  1
      49r 4.2049985e+02 6.56e+04 4.91e+08   2.2 2.38e+00   5.7 3.04e-01 6.02e-02h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      50r 5.0331875e+02 5.55e+04 2.50e+08   2.2 2.27e+00   5.2 1.36e-01 4.83e-01f  1
      51r 5.0203319e+02 5.55e+04 2.49e+08   2.2 1.75e+00   5.6 2.16e-01 6.13e-03h  1
      52r 4.9773892e+02 5.50e+04 2.37e+08   2.2 2.40e+00   5.2 9.83e-02 4.91e-02f  1
      53r 4.9110206e+02 5.47e+04 2.31e+08   2.2 2.05e+00   5.6 1.05e-01 2.43e-02h  1
      54r 4.9119651e+02 5.38e+04 2.08e+08   2.2 2.54e+00   5.1 1.33e-01 9.86e-02f  1
      55r 4.8570955e+02 5.34e+04 2.00e+08   2.2 2.20e+00   5.5 8.06e-02 3.93e-02h  1
      56r 5.3490745e+02 5.21e+04 1.67e+08   2.2 2.37e+00   5.1 1.99e-01 1.65e-01f  1
      57r 5.3351272e+02 5.20e+04 1.66e+08   2.2 1.75e+00   6.4 5.99e-02 4.73e-03h  1
      58r 5.2758534e+02 5.07e+04 1.35e+08   2.2 1.85e+00   5.9 3.65e-02 1.87e-01f  1
      59r 5.2864267e+02 5.07e+04 1.33e+08   2.2 1.60e+00   7.2 3.05e-02 9.42e-03h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      60r 5.3065412e+02 5.06e+04 1.31e+08   2.2 1.58e+00   7.7 4.34e-02 1.68e-02h  1
      61r 5.6156030e+02 4.98e+04 2.20e+08   2.2 1.56e+00   7.2 6.82e-02 1.38e-01f  1
      62r 5.6108753e+02 4.97e+04 1.13e+08   2.2 1.09e+01   6.7 7.63e-03 2.16e-02h  1
      63r 5.6142069e+02 4.97e+04 1.13e+08   2.2 1.29e+00   7.1 1.85e-01 9.92e-04h  1
      64r 6.3681457e+02 4.92e+04 9.46e+07   2.2 1.24e+00   6.7 1.50e-01 1.18e-01f  1
      65r 6.4065726e+02 4.91e+04 9.39e+07   2.2 9.88e-01   7.1 7.08e-02 7.63e-03h  1
      66r 7.0047901e+02 4.91e+04 9.38e+07   2.2 1.30e+02   6.6 1.33e-03 1.03e-03h  1
      67r 7.3881618e+02 4.90e+04 9.01e+07   2.2 1.37e+00   7.0 3.80e-01 3.94e-02h  1
      68r 7.6627609e+02 4.88e+04 8.64e+07   2.2 8.21e-01   7.5 6.37e-02 3.73e-02h  1
      69r 7.8558315e+02 4.87e+04 8.36e+07   2.2 7.90e-01   7.9 6.50e-02 2.81e-02h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      70r 8.4562914e+02 4.85e+04 7.86e+07   2.2 9.74e-01   7.4 4.33e-02 5.62e-02h  1
      71r 8.8115215e+02 4.83e+04 7.38e+07   2.2 6.70e-01   7.8 3.24e-03 4.83e-02f  1
      72r 8.8494730e+02 4.83e+04 7.36e+07   2.2 2.38e+00   7.4 1.57e-03 3.12e-03f  1
      73r 8.9885434e+02 4.82e+04 7.17e+07   2.2 6.23e-01   7.8 1.94e-02 2.28e-02h  1
      74r 9.3047892e+02 4.81e+04 6.98e+07   2.2 3.64e+00   7.3 7.57e-03 2.71e-02h  1
      75r 9.3102931e+02 4.81e+04 6.97e+07   2.2 6.08e-01   7.7 7.32e-02 1.42e-03h  1
      76r 9.4306016e+02 4.80e+04 6.61e+07   2.2 6.35e-01   8.2 9.21e-03 3.51e-02h  1
      77r 9.5207352e+02 4.80e+04 6.44e+07   2.2 6.39e-01   7.7 9.02e-02 2.51e-02h  1
      78r 9.5663048e+02 4.79e+04 6.32e+07   2.2 6.42e-01   8.1 6.69e-02 1.58e-02h  1
      79r 9.8799123e+02 4.77e+04 5.73e+07   2.2 8.09e-01   7.6 1.40e-02 8.27e-02f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      80r 9.8918022e+02 4.77e+04 5.71e+07   2.2 6.04e-01   8.1 1.73e-01 4.73e-03h  1
      81r 1.0234640e+03 4.76e+04 4.60e+08   2.2 3.97e+00   7.6 9.87e-03 3.78e-02h  1
      82r 1.0244705e+03 4.75e+04 1.75e+09   2.2 6.46e-01   8.9 1.86e-02 7.70e-03h  1
      83r 1.0269940e+03 4.75e+04 1.99e+09   2.2 5.89e-01   8.4 5.71e-02 1.91e-02f  1
      84r 1.0356507e+03 4.74e+04 2.15e+09   2.2 6.11e-01   8.0 7.90e-02 5.13e-02f  1
      85r 1.0374554e+03 4.73e+04 2.18e+09   2.2 5.69e-01   8.4 1.24e-01 1.52e-02h  1
      86r 1.0529516e+03 4.72e+04 2.21e+09   2.2 5.80e-01   7.9 1.43e-01 7.79e-02h  1
      87r 1.0579798e+03 4.71e+04 2.27e+09   2.2 5.19e-01   8.3 4.80e-01 5.02e-02h  1
      88r 1.0650820e+03 4.70e+04 2.21e+09   2.2 2.40e+00   7.9 1.57e-02 3.52e-02h  1
      89r 1.0701871e+03 4.69e+04 2.20e+09   2.2 4.89e-01   8.3 1.41e-01 6.23e-02f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      90r 1.0704407e+03 4.69e+04 2.19e+09   2.2 3.88e-01   8.7 3.26e-01 5.21e-03h  1
      91r 1.0774630e+03 4.67e+04 2.15e+09   2.2 6.59e-01   8.2 9.92e-03 1.06e-01h  1
      92r 1.0777047e+03 4.67e+04 2.14e+09   2.2 3.63e-01   8.7 2.10e-01 9.56e-03h  1
      93r 1.0777926e+03 4.66e+04 2.13e+09   2.2 3.08e-01   9.1 4.76e-03 4.15e-03h  1
      94r 1.0837850e+03 4.63e+04 2.36e+09   2.2 4.91e-01   8.6 1.07e-04 2.04e-01h  1
      95r 1.0836489e+03 4.63e+04 2.34e+09   2.2 2.05e-01   9.0 2.88e-01 9.99e-03h  1
      96r 1.0812554e+03 4.61e+04 2.50e+11   2.2 1.56e+00   8.6 1.02e-03 1.36e-01h  1
      97r 1.0812112e+03 4.61e+04 2.49e+11   2.2 3.14e-01   9.9 9.65e-03 2.55e-03h  1
      98r 1.0810610e+03 4.61e+04 2.47e+11   2.2 1.96e-01  10.3 1.01e-02 7.34e-03h  1
      99r 1.0810592e+03 4.61e+04 2.47e+11   2.2 3.76e-01   9.8 8.58e-01 1.03e-04h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     100r 1.0749950e+03 4.58e+04 1.44e+11   2.2 1.96e-01  10.3 1.24e-02 3.07e-01h  1
     101r 1.0715200e+03 4.57e+04 1.12e+11   2.2 7.98e-01   9.8 2.09e-01 8.27e-02h  1
     102r 1.0714872e+03 4.57e+04 1.12e+11   2.2 2.56e-01  10.2 1.00e+00 9.65e-04h  1
     103r 1.0645068e+03 4.56e+04 1.06e+11   2.2 2.01e-01  10.6 1.79e-01 2.24e-01h  1
     104r 1.0506682e+03 4.54e+04 7.69e+10   2.2 1.86e-01  10.2 9.90e-01 4.22e-01f  1
     105r 1.0505419e+03 4.54e+04 7.63e+10   2.2 1.78e-01  10.6 4.10e-01 7.77e-03h  1
     106r 1.0478230e+03 4.53e+04 6.52e+10   2.2 1.65e-01  10.1 8.32e-01 2.51e-01h  1
     107r 1.0476453e+03 4.53e+04 6.39e+10   2.2 1.66e-01  10.5 6.36e-01 2.03e-02h  1
     108r 1.0476508e+03 4.53e+04 6.36e+10   2.2 1.75e-01  10.1 8.65e-01 4.63e-03h  1
     109r 1.0447295e+03 4.52e+04 6.16e+10   2.2 1.63e-01  10.5 7.03e-01 5.12e-01H  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     110r 1.0500118e+03 4.52e+04 5.41e+10   2.2 2.66e-01  10.0 8.87e-01 1.19e-01f  1
     111r 1.0505347e+03 4.51e+04 5.24e+10   2.2 1.97e-01  10.4 1.00e+00 3.04e-02h  1
     112r 1.0511396e+03 4.51e+04 5.20e+10   2.2 3.56e-01  10.0 8.04e-01 7.80e-03h  1
     113r 1.0708552e+03 4.51e+04 6.57e+11   2.2 2.00e-01  10.4 1.00e+00 7.68e-01f  1
     114r 1.0714352e+03 4.51e+04 6.55e+11   2.2 4.09e-01   9.9 6.13e-01 4.17e-03h  2
     115r 1.0716176e+03 4.51e+04 6.52e+11   2.2 1.83e-01  10.3 9.31e-01 3.30e-03h  2
     116r 1.1018234e+03 4.51e+04 5.78e+11   2.2 6.62e-01   9.8 1.10e-01 1.21e-01f  1
     117r 1.1030031e+03 4.51e+04 5.71e+11   2.2 2.12e-01  10.3 1.00e+00 1.29e-02h  1
     118r 1.1032414e+03 4.51e+04 5.71e+11   2.2 2.41e+00   9.8 8.89e-02 2.54e-04h  2
     119r 1.1063382e+03 4.51e+04 5.52e+11   2.2 2.16e-01  10.2 7.05e-02 3.38e-02f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     120r 1.1064966e+03 4.51e+04 5.51e+11   2.2 8.30e-01   9.7 3.37e-03 1.33e-03h  2
     121r 1.1080895e+03 4.51e+04 5.35e+11   2.2 2.31e-01  10.2 5.00e-03 3.00e-02h  1
     122r 1.1084942e+03 4.51e+04 5.34e+11   2.2 1.59e+00   9.7 1.43e-04 1.86e-03f  2
     123r 1.1094219e+03 4.51e+04 5.25e+11   2.2 2.47e-01  10.1 2.11e-02 1.66e-02h  1
     124r 1.1095775e+03 4.51e+04 5.25e+11   2.2 1.43e+01   9.6 5.17e-04 8.26e-05h  3
     125r 1.1101603e+03 4.51e+04 5.20e+11   2.2 2.55e-01  10.1 6.87e-01 9.93e-03f  2
     126r 1.1110417e+03 4.51e+04 5.12e+11   2.2 2.70e-01  10.5 7.51e-03 1.49e-02f  1
     127r 1.1622066e+03 4.51e+04 2.64e+11   2.2 2.70e-01  10.0 2.61e-04 6.95e-01f  1
     128r 1.1623651e+03 4.51e+04 2.62e+11   2.2 5.60e-01  10.4 2.14e-02 8.25e-03h  1
     129r 1.1623655e+03 4.51e+04 2.62e+11   2.2 5.42e-01  10.9 9.02e-04 4.00e-04h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     130r 1.1633587e+03 4.51e+04 2.46e+11   2.2 5.58e-01  10.4 1.25e-03 6.60e-02f  1
     131r 1.1633923e+03 4.51e+04 2.45e+11   2.2 5.12e-01  10.8 4.53e-02 1.67e-03h  1
     132r 1.1633956e+03 4.51e+04 2.45e+11   2.2 1.75e+00  10.3 3.61e-03 8.10e-05h  1
     133r 1.1634035e+03 4.51e+04 2.45e+11   2.2 5.13e-01  10.8 1.84e-01 3.45e-04h  1
     134r 1.1682801e+03 4.51e+04 2.26e+11   2.2 1.74e+00  10.3 1.18e-02 8.30e-02f  1
     135r 1.1683637e+03 4.51e+04 2.26e+11   2.2 8.02e-01  10.7 2.66e-03 1.08e-03h  1
     136r 1.1683721e+03 4.51e+04 2.26e+11   2.2 4.71e-01  11.1 2.94e-05 1.38e-04h  1
     137r 1.1809007e+03 4.51e+04 1.97e+11   2.2 9.06e-01  10.7 2.25e-01 1.39e-01h  1
     138r 1.2326511e+03 4.51e+04 1.47e+11   2.2 4.59e-01  10.2 5.06e-02 3.22e-01f  1
     139r 1.2332547e+03 4.51e+04 1.47e+11   2.2 3.93e-01   9.7 1.98e-03 3.52e-03h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     140r 1.2332606e+03 4.51e+04 1.47e+11   2.2 4.01e-01  10.1 2.51e-02 5.98e-05h  1
     141r 1.2332288e+03 4.51e+04 1.47e+11   2.2 4.71e-01   9.7 6.85e-04 1.81e-03h  1
     142r 1.2230516e+03 4.51e+04 6.88e+10   2.2 4.70e-01  10.1 7.49e-06 1.00e+00f  1
     143r 1.2229458e+03 4.51e+04 6.34e+10   2.2 2.21e-02  10.5 1.09e-01 8.27e-02h  1
     144r 1.2229417e+03 4.51e+04 6.33e+10   2.2 4.65e-02  10.0 9.22e-01 2.32e-03h  1
     145r 1.2216410e+03 4.51e+04 2.69e+10   2.2 2.11e-02  10.5 1.60e-01 1.00e+00f  1
     146r 1.2203021e+03 4.51e+04 1.20e+10   2.2 2.72e-02  10.0 9.84e-01 8.92e-01f  1
     147r 1.2177897e+03 4.51e+04 5.34e+09   2.2 8.50e-02   9.5 1.00e+00 1.00e+00f  1
     148r 1.2152909e+03 4.51e+04 2.36e+09   2.2 1.93e-01   9.0 1.00e+00 1.00e+00f  1
     149r 1.2162988e+03 4.51e+04 2.18e+09   2.2 3.61e-01   8.6 2.26e-01 6.10e-02h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     150r 1.2163040e+03 4.51e+04 2.17e+09   2.2 1.89e-01   9.0 4.81e-03 3.32e-03h  1
     151r 1.2156715e+03 4.51e+04 5.17e+09   2.2 5.01e-01   8.5 5.99e-01 3.18e-01f  2
     152r 1.2156963e+03 4.51e+04 4.35e+09   2.2 2.41e-02   9.8 1.00e+00 1.75e-01h  1
     153r 1.2153965e+03 4.51e+04 2.58e+09   2.2 5.50e-02   9.4 1.00e+00 7.20e-01f  1
     154r 1.2150628e+03 4.51e+04 1.33e+09   2.2 1.62e-02   9.8 1.00e+00 1.00e+00f  1
     155r 1.2148137e+03 4.51e+04 5.24e+08   2.2 1.40e-02   9.3 1.64e-01 1.00e+00f  1
     156r 1.2145996e+03 4.51e+04 7.07e+07   2.2 4.36e-03   8.8 1.00e+00 1.00e+00f  1
     157r 1.2141544e+03 4.51e+04 7.80e+05   2.2 3.99e-04   8.3 1.00e+00 1.00e+00f  1
     158r 1.2132011e+03 4.51e+04 5.80e+04   2.2 7.78e-04   7.9 1.00e+00 1.00e+00f  1
     159r 1.2112218e+03 4.51e+04 3.64e+04   2.2 1.46e-03   7.4 1.00e+00 1.00e+00f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     160r 1.2071203e+03 4.51e+04 2.19e+04   2.2 2.65e-03   6.9 1.00e+00 1.00e+00f  1
     161r 1.2052793e+03 4.51e+04 2.65e+03   1.5 9.62e-04   6.4 1.00e+00 1.00e+00f  1
     162r 1.1987948e+03 4.51e+04 3.64e+03   1.5 3.95e-03   6.0 1.00e+00 1.00e+00f  1
     163r 1.1857268e+03 4.51e+04 1.72e+03   1.5 5.62e-03   5.5 1.00e+00 1.00e+00f  1
     164r 1.1532596e+03 4.51e+04 1.29e+03   1.5 1.26e-02   5.0 1.00e+00 1.00e+00f  1
     165r 1.0780079e+03 4.50e+04 1.21e+03   1.5 3.56e-02   4.5 1.00e+00 1.00e+00f  1
     166r 8.9395409e+02 4.50e+04 1.20e+03   1.5 1.06e-01   4.1 1.00e+00 1.00e+00f  1
     167r 6.1535729e+02 4.50e+04 1.16e+03   1.5 3.06e-01   3.6 1.00e+00 6.06e-01f  1
     168r 5.5170374e+02 4.50e+04 1.16e+03   0.8 1.18e-01   4.0 1.00e+00 5.51e-01f  1
     169r 5.1472965e+02 4.50e+04 1.20e+03   0.8 4.45e-02   4.4 1.00e+00 1.00e+00f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     170r 4.5301207e+02 4.49e+04 1.25e+03   0.8 1.47e-01   4.0 1.00e+00 4.79e-01f  1
     171r 4.1744081e+02 4.49e+04 1.29e+03   0.8 5.47e-02   4.4 1.00e+00 8.74e-01f  1
     172r 3.6441865e+02 4.49e+04 1.28e+03   0.8 2.12e-01   3.9 1.00e+00 5.59e-01f  1
     173r 3.4795862e+02 4.49e+04 1.28e+03   0.8 7.37e-01   3.4 4.69e-01 8.46e-02f  1
     174r 3.2470699e+02 4.49e+04 1.24e+03   0.8 1.64e+00   2.9 1.54e-01 1.91e-01f  1
     175r 3.4644042e+02 4.48e+04 1.27e+03   0.8 6.50e-01   3.4 8.55e-01 2.83e-01f  1
     176r 3.5603098e+02 4.48e+04 8.29e+03   0.8 2.49e-01   3.8 1.00e+00 5.09e-01f  1
     177r 3.6952634e+02 4.48e+04 8.27e+03   0.8 7.57e-01   3.3 8.09e-01 2.84e-01f  1
     178r 3.7786810e+02 4.47e+04 7.08e+03   0.8 2.79e-01   3.7 1.00e+00 3.29e-01f  1
     179r 3.8789944e+02 4.47e+04 6.42e+03   0.8 2.06e+00   3.3 1.12e-01 1.89e-01f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     180r 4.0881970e+02 4.47e+04 5.95e+03   0.8 1.78e+00   2.8 4.51e-01 1.12e-01f  1
     181r 4.3151463e+02 4.46e+04 5.31e+03   0.8 6.59e-01   3.2 1.00e+00 2.29e-01f  1
     182r 4.8901280e+02 4.46e+04 4.73e+03   0.8 1.90e+00   2.7 2.73e-01 1.85e-01f  1
     183r 5.1937777e+02 4.45e+04 4.16e+03   0.8 7.18e-01   3.2 1.00e+00 1.85e-01f  1
     184r 6.2662901e+02 4.44e+04 2.13e+04   0.8 2.07e+00   2.7 4.10e-01 1.81e-01f  1
     185r 6.4696846e+02 4.44e+04 2.15e+04   0.8 1.30e-01   4.0 1.00e+00 6.02e-01f  1
     186r 6.4756866e+02 4.44e+04 3.16e+03   0.8 2.46e-02   5.3 8.10e-01 5.72e-01f  1
     187r 6.5182926e+02 4.44e+04 2.52e+04   0.8 1.74e-02   4.9 1.00e+00 1.00e+00f  1
     188r 6.5396605e+02 4.44e+04 2.36e+04   0.8 6.33e-03   5.3 1.00e+00 1.00e+00f  1
     189r 6.5782165e+02 4.44e+04 2.02e+04   0.8 2.81e-02   4.8 1.00e+00 1.00e+00f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     190r 6.7526981e+02 4.44e+04 4.80e+03   0.8 4.79e-02   4.3 1.00e+00 9.64e-01f  1
     191r 7.0053119e+02 4.43e+04 2.75e+03   0.8 2.18e-01   3.9 1.00e+00 6.03e-01f  1
     192r 7.0428560e+02 4.43e+04 2.01e+03   0.8 2.54e-01   4.3 2.09e-01 9.54e-02f  1
     193r 7.1548182e+02 4.43e+04 5.78e+03   0.8 4.69e-01   3.8 5.15e-01 1.57e-01f  1
     194r 7.8720297e+02 4.43e+04 4.44e+03   0.8 6.39e-01   3.3 1.00e+00 3.58e-01f  1
     195r 8.0733650e+02 4.43e+04 2.34e+03   0.8 1.79e-01   3.8 1.00e+00 4.56e-01f  1
     196r 8.4622363e+02 4.42e+04 1.55e+03   0.8 5.28e-01   3.3 1.00e+00 3.40e-01f  1
     197r 8.8519227e+02 4.41e+04 1.19e+03   0.8 1.57e+00   2.8 1.00e+00 2.28e-01f  1
     198r 9.2236617e+02 4.41e+04 1.01e+03   0.8 5.88e-01   3.2 1.00e+00 5.18e-01f  1
     199r 9.1922287e+02 4.40e+04 1.02e+03   0.8 1.78e+00   2.8 5.85e-01 1.04e-01f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     200r 9.2366201e+02 4.39e+04 1.02e+03   0.8 6.65e-01   3.2 9.17e-01 3.99e-01f  1
     201r 9.0710015e+02 4.39e+04 1.02e+03   0.8 2.00e+00   2.7 6.29e-01 1.56e-01f  1
     202r 9.0460720e+02 4.38e+04 1.02e+03   0.8 1.44e+00   3.1 4.11e-01 1.54e-01f  1
     203r 8.7583868e+02 4.37e+04 1.02e+03   0.8 2.24e+00   2.7 5.25e-01 2.15e-01f  1
     204r 8.1296283e+02 4.36e+04 1.02e+03   0.8 1.31e+01   2.2 9.45e-02 3.81e-02f  1
     205r 8.0312228e+02 4.35e+04 1.02e+03   0.8 2.51e+00   2.6 3.62e-01 1.85e-01f  1
     206r 8.0146299e+02 4.34e+04 1.02e+03   0.8 9.43e-01   3.0 1.00e+00 4.94e-01f  1
     207r 7.8386209e+02 4.32e+04 1.02e+03   0.8 2.82e+00   2.6 8.95e-01 2.14e-01f  1
     208r 7.7888180e+02 4.32e+04 1.01e+03   0.8 1.06e+00   3.0 1.00e+00 2.34e-01f  1
     209r 7.6167609e+02 4.30e+04 1.01e+03   0.8 3.17e+00   2.5 1.00e+00 1.96e-01f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     210r 7.5467413e+02 4.29e+04 1.01e+03   0.8 1.19e+00   2.9 1.00e+00 3.39e-01f  1
     211r 7.3934325e+02 4.27e+04 1.01e+03   0.8 3.56e+00   2.5 1.00e+00 1.73e-01f  1
     212r 7.1627533e+02 4.24e+04 1.01e+03   0.8 1.07e+01   2.0 1.00e+00 1.23e-01f  1
     213r 7.0175397e+02 4.21e+04 1.01e+03   0.8 4.00e+00   2.4 1.00e+00 1.99e-01f  1
     214r 6.3706391e+02 4.02e+04 1.43e+03   0.8 1.20e+01   1.9 8.36e-01 7.37e-01f  1
     215r 6.1768383e+02 3.94e+04 1.28e+03   0.8 3.59e+01   1.4 7.54e-02 8.43e-02f  1
     216r 6.1939169e+02 3.72e+04 1.01e+03   0.8 1.35e+01   1.9 1.00e+00 5.84e-01f  1
     217r 6.4468882e+02 3.59e+04 1.01e+03   0.8 5.05e+00   2.3 1.00e+00 7.42e-01f  1
     218r 6.4644350e+02 3.48e+04 1.01e+03   0.8 1.52e+01   1.8 1.00e+00 2.54e-01f  1
     219r 6.6042912e+02 3.32e+04 1.00e+03   0.8 5.67e+00   2.2 1.00e+00 1.00e+00f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     220r 6.6586984e+02 3.07e+04 1.00e+03   0.8 1.70e+01   1.8 8.83e-01 5.53e-01f  1
     221r 6.6400989e+02 2.90e+04 1.01e+03   0.8 6.38e+00   2.2 1.00e+00 1.00e+00f  1
     222r 6.4380384e+02 2.37e+04 3.16e+05   0.8 1.92e+01   1.7 7.23e-01 1.00e+00f  1
     223r 6.4065848e+02 2.36e+04 4.56e+04   0.8 8.99e-01   3.0 1.00e+00 7.96e-01f  1
     224r 6.4848229e+02 2.34e+04 9.63e+05   0.8 2.90e+00   2.6 7.19e-02 2.19e-01f  1
     225r 6.3765147e+02 2.32e+04 4.30e+05   0.8 1.03e+00   3.0 7.06e-01 1.00e+00f  1
     226r 6.3734304e+02 2.31e+04 4.16e+05   0.8 3.07e+00   2.5 5.32e-01 3.45e-02f  1
     227r 6.1816303e+02 2.15e+04 5.16e+06   0.8 9.19e+00   2.0 1.00e+00 6.68e-01f  1
     228r 6.1621293e+02 2.15e+04 3.19e+06   0.8 8.72e-02   4.3 6.09e-01 3.83e-01f  1
     229r 6.1509314e+02 2.15e+04 1.23e+06   0.8 2.79e-01   3.8 4.31e-01 4.62e-01f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     230r 6.1195120e+02 2.15e+04 1.55e+06   0.8 2.82e-02   5.1 1.00e+00 8.82e-01h  1
     231r 6.1370097e+02 2.15e+04 7.14e+05   0.8 4.05e-02   4.7 1.00e+00 1.00e+00f  1
     232r 6.1460337e+02 2.15e+04 4.33e+05   0.8 4.63e-02   5.1 1.00e+00 6.54e-01h  1
     233r 6.1487357e+02 2.15e+04 3.84e+05   0.8 9.73e-02   5.5 6.10e-01 1.32e-01h  1
     234r 6.1547680e+02 2.15e+04 3.16e+05   0.8 2.22e-01   5.0 2.54e-01 2.14e-01h  1
     235r 6.1514611e+02 2.15e+04 2.91e+05   0.8 8.96e-01   4.5 2.58e-02 8.53e-02f  1
     236r 6.1477430e+02 2.15e+04 2.80e+05   0.8 7.48e-01   5.0 7.46e-02 3.93e-02h  1
     237r 6.1474651e+02 2.15e+04 2.79e+05   0.8 9.23e-01   5.4 2.27e-02 2.44e-03h  1
     238r 6.1434726e+02 2.15e+04 2.67e+05   0.8 4.76e-01   5.8 1.32e-02 4.73e-02h  1
     239r 6.1429355e+02 2.15e+04 2.54e+05   0.8 2.26e-01   6.3 1.31e-02 5.96e-02h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     240r 6.1440612e+02 2.15e+04 2.35e+05   0.8 9.56e-02   5.8 9.23e-02 9.09e-02f  1
     241r 6.1452299e+02 2.15e+04 2.90e+05   0.8 4.95e-02   5.3 1.23e-02 2.27e-01f  1
     242r 6.1452055e+02 2.15e+04 2.88e+05   0.8 1.33e-01   6.6 1.42e-01 5.96e-03h  1
     243r 6.1449481e+02 2.15e+04 2.57e+05   0.8 1.68e-01   6.1 2.02e-01 9.84e-02h  1
     244r 6.1422464e+02 2.15e+04 2.69e+05   0.8 1.82e-01   6.6 1.78e-02 9.78e-02H  1
     245r 6.1421220e+02 2.15e+04 2.69e+05   0.8 2.16e+00   6.1 1.40e-02 1.58e-03h  1
     246r 6.1420545e+02 2.15e+04 2.68e+05   0.8 7.21e-01   6.5 6.01e-04 1.54e-03h  1
     247r 6.1420486e+02 2.15e+04 2.68e+05   0.8 5.34e+00   6.0 2.87e-03 7.25e-05h  1
     248r 6.1416538e+02 2.15e+04 2.66e+05   0.8 7.86e-01   6.5 5.38e-04 8.02e-03f  1
     249r 6.1416141e+02 2.15e+04 2.66e+05   0.8 9.32e-01   6.9 3.33e-02 7.06e-04h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     250r 6.1409969e+02 2.15e+04 2.63e+05   0.8 9.83e-01   6.4 3.56e-03 1.04e-02h  1
     251r 6.1407542e+02 2.15e+04 8.52e+05   0.8 9.83e-01   6.8 1.02e-01 4.11e-03h  1
     252r 6.1400538e+02 2.15e+04 8.34e+05   0.8 9.28e-01   6.4 9.31e-03 1.39e-02h  1
     253r 6.1398392e+02 2.15e+04 1.18e+06   0.8 9.65e-01   6.8 7.25e-02 3.97e-03h  1
     254r 6.1356771e+02 2.15e+04 1.02e+06   0.8 9.16e-01   6.3 3.01e-02 9.95e-02h  1
     255r 6.1342557e+02 2.15e+04 9.33e+05   0.8 8.39e-01   5.8 9.50e-02 6.52e-02h  1
     256r 6.1325729e+02 2.15e+04 9.32e+05   0.8 8.08e-01   6.3 1.54e-04 4.70e-04H  1
     257r 6.1325551e+02 2.15e+04 9.24e+05   0.8 8.04e-01   5.8 1.15e-02 2.94e-04h  3
     258r 6.1322699e+02 2.15e+04 9.22e+05   0.8 8.06e-01   6.2 3.00e-02 3.97e-03h  1
     259r 6.1300914e+02 2.15e+04 9.05e+05   0.8 1.01e+00   5.7 1.68e-02 3.10e-02h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     260r 6.1285508e+02 2.15e+04 8.78e+05   0.8 9.71e-01   6.2 1.10e-01 2.06e-02h  1
     261r 6.1206020e+02 2.15e+04 7.73e+05   0.8 8.66e-01   5.7 1.17e-01 1.90e-01h  1
     262r 6.1141418e+02 2.15e+04 5.42e+05   0.8 7.61e-01   6.1 5.02e-01 1.66e-01h  1
     263r 6.1025640e+02 2.16e+04 2.99e+05   0.8 6.52e-01   5.6 6.43e-01 6.65e-01h  1
     264r 6.0970073e+02 2.15e+04 1.92e+05   0.8 2.27e-01   5.2 3.45e-01 3.94e-01h  1
     265r 6.0969601e+02 2.15e+04 1.92e+05   0.8 2.12e+00   5.6 3.14e-03 1.17e-03h  1
     266r 6.0948787e+02 2.15e+04 1.88e+05   0.8 8.17e-01   6.0 2.75e-02 4.73e-02h  1
     267r 6.0947455e+02 2.15e+04 1.87e+05   0.8 7.21e-01   5.5 3.57e-04 4.37e-03h  1
     268r 6.0947415e+02 2.15e+04 1.87e+05   0.8 8.08e-01   6.0 1.02e-01 1.04e-04h  1
     269r 6.0926502e+02 2.15e+04 1.81e+05   0.8 7.68e-01   5.5 2.81e-03 4.38e-02h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     270r 6.0885308e+02 2.15e+04 1.53e+05   0.8 9.83e-01   5.0 2.36e-01 1.62e-01h  1
     271r 6.0863272e+02 2.15e+04 1.40e+05   0.8 7.76e-01   5.4 2.72e-01 1.02e-01h  1
     272r 6.0867346e+02 2.16e+04 6.34e+04   0.8 7.23e-01   5.0 1.00e+00 8.61e-01h  1
     273r 6.1033541e+02 2.15e+04 2.54e+04   0.8 1.10e-01   4.5 1.00e+00 6.32e-01f  1
     274r 6.1301815e+02 2.14e+04 2.96e+03   0.8 1.01e-01   4.0 1.00e+00 9.02e-01f  1
     275r 6.1276999e+02 2.14e+04 1.01e+03   0.8 3.02e-01   3.5 1.00e+00 7.83e-01f  1
     276r 6.1271461e+02 2.14e+04 1.01e+03   0.8 1.13e-01   4.0 1.00e+00 1.00e+00f  1
     277r 6.1108024e+02 2.13e+04 1.01e+03   0.8 3.40e-01   3.5 1.00e+00 1.00e+00f  1
     278r 6.1035278e+02 2.13e+04 1.01e+03   0.8 1.02e+00   3.0 1.00e+00 1.34e-01f  1
     279r 6.0399287e+02 2.10e+04 3.44e+07   0.8 3.04e+00   2.5 5.59e-01 2.89e-01f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     280r 6.0119263e+02 2.09e+04 2.24e+07   0.8 1.14e+00   2.9 1.00e+00 3.92e-01f  1
     281r 6.0032086e+02 2.08e+04 1.47e+07   0.8 4.29e-01   3.4 1.00e+00 3.95e-01f  1
     282r 5.9235961e+02 2.05e+04 7.26e+06   0.8 1.28e+00   2.9 7.89e-01 1.00e+00f  1
     283r 5.9164724e+02 2.05e+04 3.61e+06   0.8 6.02e-02   4.2 1.00e+00 1.00e+00f  1
     284r 5.9361028e+02 2.04e+04 3.29e+06   0.8 2.79e-01   3.7 1.87e-01 2.30e-01f  1
     285r 5.9335553e+02 2.04e+04 3.20e+06   0.8 7.09e-01   4.2 3.50e-02 2.86e-02h  1
     286r 5.9297638e+02 2.04e+04 3.11e+06   0.8 9.07e-01   4.6 1.13e-02 3.06e-02h  1
     287r 5.9302455e+02 2.05e+04 2.63e+06   0.8 1.14e+00   4.1 4.75e-02 2.33e-01h  1
     288r 5.9276024e+02 2.05e+04 2.48e+06   0.8 7.21e-01   4.5 7.87e-02 6.41e-02h  1
     289r 5.9319922e+02 2.05e+04 2.34e+06   0.8 7.03e-01   4.1 1.36e-02 7.53e-02h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     290r 5.9270943e+02 2.05e+04 2.01e+06   0.8 9.45e-01   4.5 5.08e-02 1.68e-01h  1
     291r 5.9435832e+02 2.05e+04 2.00e+06   0.8 8.09e-01   4.0 1.61e-01 9.61e-02h  2
     292r 5.9418723e+02 2.05e+04 1.89e+06   0.8 1.22e+00   4.4 3.00e-02 6.26e-02h  1
     293r 5.9456885e+02 2.05e+04 1.86e+06   0.8 4.66e+00   4.9 1.09e-02 1.32e-02h  1
     294r 5.9516283e+02 2.05e+04 1.83e+06   0.8 1.50e+00   5.3 3.66e-02 1.86e-02h  1
     295r 5.9506092e+02 2.05e+04 1.77e+06   0.8 1.18e+00   5.7 1.27e-02 3.34e-02h  1
     296r 5.9423599e+02 2.05e+04 1.59e+06   0.8 1.06e+00   5.2 6.27e-04 1.10e-01h  1
     297r 5.9422249e+02 2.05e+04 1.59e+06   0.8 7.14e-01   5.7 3.20e-02 1.66e-03h  1
     298r 5.9401474e+02 2.05e+04 1.56e+06   0.8 1.08e+00   5.2 1.24e-02 2.12e-02h  1
     299r 5.9397460e+02 2.05e+04 1.55e+06   0.8 1.06e+00   5.6 7.19e-02 3.68e-03h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     300r 5.9316114e+02 2.05e+04 1.47e+06   0.8 1.07e+00   5.1 3.25e-02 5.53e-02h  1
     301r 5.9219548e+02 2.05e+04 1.37e+06   0.8 9.88e-01   5.6 1.02e-01 7.41e-02h  1
     302r 5.9137188e+02 2.05e+04 1.30e+06   0.8 8.36e-01   6.0 2.41e-02 5.79e-02h  1
     303r 5.9129956e+02 2.05e+04 1.29e+06   0.8 6.49e-01   6.4 2.59e-03 4.34e-03h  1
     304r 5.9115157e+02 2.05e+04 1.28e+06   0.8 5.93e-01   6.8 2.61e-01 1.10e-02h  1
     305r 5.8964514e+02 2.05e+04 1.26e+06   0.8 1.02e+01   6.4 8.73e-04 1.55e-02h  1
     306r 5.8970965e+02 2.05e+04 1.26e+06   0.8 8.52e-01   5.9 2.30e-02 2.04e-03H  1
     307r 5.8825525e+02 2.04e+04 1.09e+06   0.8 9.75e-01   6.3 5.18e-04 1.90e-01h  1
     308r 5.8820886e+02 2.04e+04 1.16e+06   0.8 7.69e-01   7.6 2.06e-02 1.19e-02h  1
     309r 5.8817682e+02 2.04e+04 1.85e+06   0.8 6.72e-01   8.1 1.30e-04 1.13e-02h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     310r 5.8816726e+02 2.04e+04 3.32e+06   0.8 4.80e-01   8.5 1.17e-03 6.80e-03h  1
     311r 5.8816399e+02 2.04e+04 4.41e+06   0.8 4.41e-01   8.9 8.19e-03 2.87e-03h  1
     312r 5.8819969e+02 2.04e+04 1.47e+07   0.8 7.81e-01   8.4 6.90e-03 1.81e-02h  1
     313r 5.8831563e+02 2.04e+04 1.20e+07   0.8 6.68e-01   8.0 5.01e-03 1.85e-02h  1
     314r 5.8879132e+02 2.04e+04 3.57e+07   0.8 8.82e-01   7.5 4.59e-03 4.72e-02h  1
     315r 5.8878945e+02 2.04e+04 3.57e+07   0.8 3.04e-01   8.8 5.24e-02 2.51e-03h  1
     316r 5.8878017e+02 2.04e+04 3.49e+07   0.8 3.04e-01   9.2 2.87e-02 1.26e-02h  1
     317r 5.8877085e+02 2.04e+04 3.48e+07   0.8 3.00e-01   8.8 7.76e-04 1.76e-02h  1
     318r 5.8876745e+02 2.04e+04 3.48e+07   0.8 2.95e-01   9.2 8.55e-03 6.53e-03h  1
     319r 5.8876324e+02 2.04e+04 3.48e+07   0.8 2.94e-01   8.7 5.52e-04 1.15e-02h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     320r 5.8876151e+02 2.04e+04 3.48e+07   0.8 2.90e-01   9.1 8.75e-03 3.71e-03h  1
     321r 5.8875769e+02 2.04e+04 3.53e+07   0.8 2.89e-01   8.7 6.36e-03 1.67e-02h  1
     322r 5.8875406e+02 2.04e+04 3.52e+07   0.8 2.85e-01   9.1 2.96e-02 9.98e-03h  1
     323r 5.8875660e+02 2.04e+04 4.33e+07   0.8 2.82e-01   8.6 1.84e-03 2.04e-02h  1
     324r 5.8875456e+02 2.04e+04 4.26e+07   0.8 2.76e-01   9.0 6.09e-02 1.42e-02h  1
     325r 5.8875135e+02 2.04e+04 1.84e+08   0.8 2.73e-01   9.5 2.23e-01 1.66e-02h  1
     326r 5.8875005e+02 2.04e+04 1.81e+08   0.8 2.68e-01   9.0 1.36e-01 3.85e-02h  1
     327r 5.8879870e+02 2.04e+04 3.44e+08   0.8 2.59e-01   8.5 1.29e-02 4.13e-02H  1
     328r 5.8880355e+02 2.04e+04 3.33e+08   0.8 2.50e-01   8.0 3.68e-03 2.99e-02h  1
     329r 5.8880196e+02 2.04e+04 2.91e+08   0.8 2.42e-01   8.5 1.55e-01 1.25e-01h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     330r 5.8879363e+02 2.04e+04 2.84e+08   0.8 3.50e-01   8.0 1.39e-02 2.25e-02h  1
     331r 5.8878758e+02 2.04e+04 2.82e+08   0.8 4.53e-01   8.4 4.25e-02 9.34e-03h  1
     332r 5.8878767e+02 2.04e+04 2.82e+08   0.8 2.05e-01   8.8 5.12e-01 6.48e-04h  1
     333r 5.8874127e+02 2.04e+04 2.74e+08   0.8 1.41e+00   8.4 1.92e-02 2.64e-02h  1
     334r 5.8873842e+02 2.04e+04 2.73e+08   0.8 8.38e-01   8.8 1.84e-01 3.39e-03h  1
     335r 5.8872531e+02 2.04e+04 2.69e+08   0.8 8.32e-01   8.3 2.59e-02 1.68e-02h  1
     336r 5.8872438e+02 2.04e+04 2.68e+08   0.8 8.18e-01   7.8 4.07e-01 7.54e-04h  1
     337r 5.8869765e+02 2.04e+04 2.55e+08   0.8 8.20e-01   8.3 2.55e-02 4.99e-02h  1
     338r 5.8877488e+02 2.04e+04 2.52e+08   0.8 1.39e+00   8.7 2.78e-03 1.26e-02h  1
     339r 5.8878663e+02 2.04e+04 2.51e+08   0.8 2.74e-01   8.2 7.85e-03 3.56e-03h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     340r 5.8971730e+02 2.04e+04 6.81e+08   0.8 3.82e+00   7.7 4.36e-03 3.41e-02f  1
     341r 5.8960647e+02 2.04e+04 6.66e+08   0.8 3.49e-01   8.2 1.44e-01 1.90e-02h  1
     342r 5.8960638e+02 2.04e+04 6.66e+08   0.8 1.73e-01   8.6 1.01e-01 6.34e-04h  1
     343r 5.8946978e+02 2.04e+04 6.58e+08   0.8 6.05e-01   8.1 8.79e-04 1.21e-02h  1
     344r 5.8938455e+02 2.04e+04 4.99e+08   0.8 2.00e-01   8.5 4.29e-01 1.48e-01h  1
     345r 5.8946041e+02 2.04e+04 4.42e+08   0.8 2.91e-01   9.0 1.95e-03 5.24e-02h  1
     346r 5.8959637e+02 2.04e+04 4.12e+08   0.8 2.75e-01   9.4 2.59e-05 4.37e-02h  1
     347r 5.8973991e+02 2.04e+04 4.44e+08   0.8 3.64e-01   9.8 3.24e-03 3.36e-02h  1
     348r 5.8974767e+02 2.04e+04 4.44e+08   0.8 3.89e-01  10.2 2.44e-03 1.69e-03h  1
     349r 5.8980114e+02 2.04e+04 4.50e+08   0.8 6.03e-01   9.8 1.07e-02 8.89e-03h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     350r 5.8982124e+02 2.04e+04 4.54e+08   0.8 4.13e-01  10.2 4.86e-03 4.26e-03h  1
     351r 5.8990462e+02 2.04e+04 4.71e+08   0.8 5.83e-01   9.7 1.80e-04 1.34e-02h  1
     352r 5.8990533e+02 2.04e+04 4.71e+08   0.8 3.19e+00   9.2 1.05e-03 2.41e-05h  1
     353r 5.8990545e+02 2.04e+04 4.71e+08   0.8 4.67e-01   9.7 8.85e-03 2.30e-05h  1
     354r 5.8994873e+02 2.04e+04 4.71e+08   0.8 5.11e+00   9.2 3.57e-05 8.55e-04f  1
     355r 5.8994934e+02 2.04e+04 4.71e+08   0.8 4.47e-01   9.6 4.34e-02 1.06e-04h  1
     356r 5.8995815e+02 2.04e+04 4.71e+08   0.8 2.90e+00   9.1 1.89e-03 2.14e-04h  1
     357r 5.9007237e+02 2.04e+04 4.87e+08   0.8 4.17e-01   9.6 2.77e-02 1.95e-02h  1
     358r 5.9007347e+02 2.04e+04 4.87e+08   0.8 4.57e-01   9.1 4.50e-02 1.79e-04h  1
     359r 5.9007372e+02 2.04e+04 4.87e+08   0.8 7.40e-01   9.5 7.51e-05 6.65e-05h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     360r 5.9057257e+02 2.04e+04 5.48e+08   0.8 5.70e-01   9.0 8.58e-06 6.01e-02h  1
     361r 5.9057086e+02 2.04e+04 5.47e+08   0.8 9.53e-01   8.6 5.22e-02 1.04e-03h  1
     362r 5.9059746e+02 2.04e+04 5.44e+08   0.8 9.13e-01   8.1 8.42e-03 7.42e-03h  1
     363r 5.9049199e+02 2.04e+04 5.21e+08   0.8 9.54e-01   8.5 9.76e-02 5.01e-02h  1
     364r 5.9034117e+02 2.04e+04 5.12e+08   0.8 8.53e-01   8.0 1.63e-02 1.84e-02h  1
     365r 5.9028722e+02 2.04e+04 5.06e+08   0.8 8.13e-01   8.4 1.03e-01 1.33e-02h  1
     366r 5.9024781e+02 2.04e+04 5.06e+08   0.8 1.12e+02   8.0 5.28e-05 9.93e-05f  2
     367r 5.9024669e+02 2.04e+04 5.05e+08   0.8 7.56e-01   8.4 1.09e-04 2.63e-04h  1
     368r 5.9022646e+02 2.04e+04 5.01e+08   0.8 7.88e-01   8.8 4.95e-02 1.05e-02h  1
     369r 5.9021944e+02 2.04e+04 5.00e+08   0.8 1.00e+00   8.3 1.41e-03 1.21e-03h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     370r 5.9019496e+02 2.04e+04 4.94e+08   0.8 7.37e-01   8.8 2.28e-03 1.89e-02h  1
     371r 5.9019470e+02 2.04e+04 4.94e+08   0.8 6.84e-01   9.2 1.57e-02 5.39e-04h  1
     372r 5.9011320e+02 2.04e+04 5.17e+08   0.8 5.68e-01   8.7 1.32e-03 4.43e-02h  1
     373r 5.9012434e+02 2.04e+04 5.20e+08   0.8 6.28e-01   9.1 1.45e-02 1.28e-02h  1
     374r 5.9018273e+02 2.04e+04 5.79e+08   0.8 3.97e-01   9.6 8.80e-05 2.75e-02h  1
     375r 5.9037999e+02 2.04e+04 3.91e+10   0.8 9.42e+00   9.1 2.16e-05 1.09e-02f  1
     376r 5.9044214e+02 2.04e+04 3.86e+10   0.8 3.04e-01  10.4 5.10e-02 1.29e-02h  1
     377r 5.9044156e+02 2.04e+04 3.85e+10   0.8 5.27e-01   9.9 4.94e-04 1.63e-03h  1
     378r 5.9044229e+02 2.04e+04 3.85e+10   0.8 3.71e-01  10.4 1.44e-05 1.84e-04h  1
     379r 5.9044229e+02 2.04e+04 3.85e+10   0.8 8.81e-01   9.9 0.00e+00 2.43e-07R  2
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     380r 5.9044230e+02 2.04e+04 2.23e+06   0.8 6.62e-08  10.3 9.90e-01 1.00e+00f  1
     381r 5.9044233e+02 2.04e+04 2.20e+04   0.8 1.47e-07   9.8 9.90e-01 1.00e+00f  1
     382r 5.9044240e+02 2.04e+04 1.03e+03   0.8 4.40e-07   9.4 9.90e-01 1.00e+00f  1
     383r 5.9044263e+02 2.04e+04 1.03e+03   0.8 1.32e-06   8.9 1.00e+00 1.00e+00f  1
     384r 5.9044332e+02 2.04e+04 1.03e+03   0.8 3.96e-06   8.4 1.00e+00 1.00e+00f  1
     385r 5.9044537e+02 2.04e+04 1.03e+03   0.8 1.19e-05   7.9 1.00e+00 1.00e+00f  1
     386r 5.9045151e+02 2.04e+04 1.03e+03   0.8 3.56e-05   7.5 1.00e+00 1.00e+00f  1
     387r 5.9046967e+02 2.04e+04 1.03e+03   0.8 1.07e-04   7.0 1.00e+00 1.00e+00f  1
     388r 5.9052200e+02 2.04e+04 1.02e+03   0.8 3.20e-04   6.5 1.00e+00 1.00e+00f  1
     389r 5.9066546e+02 2.04e+04 1.02e+03   0.8 9.57e-04   6.0 1.00e+00 1.00e+00f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     390r 5.9103472e+02 2.04e+04 1.02e+03   0.8 2.87e-03   5.6 1.00e+00 1.00e+00f  1
     391r 5.9194135e+02 2.04e+04 1.02e+03   0.8 8.61e-03   5.1 1.00e+00 1.00e+00f  1
     392r 5.9379242e+02 2.04e+04 1.02e+03   0.8 2.59e-02   4.6 1.00e+00 1.00e+00f  1
     393r 5.9383145e+02 2.04e+04 1.16e+06   0.8 1.60e-01   4.1 8.68e-01 4.80e-01f  1
     394r 5.9384723e+02 2.04e+04 8.94e+05   0.8 3.06e-02   5.4 1.16e-01 2.27e-01f  1
     395r 5.9387095e+02 2.04e+04 3.56e+05   0.8 7.84e-02   5.0 7.02e-01 7.01e-01f  1
     396r 5.9213744e+02 2.04e+04 2.02e+05   0.8 9.81e-02   4.5 4.03e-01 1.00e+00f  1
     397r 5.8925556e+02 2.04e+04 1.78e+05   0.8 2.27e-01   4.0 1.41e-01 3.30e-01f  1
     398r 5.8813815e+02 2.04e+04 1.03e+05   0.8 5.17e-02   4.4 1.00e+00 1.00e+00f  1
     399r 5.9000330e+02 2.04e+04 6.04e+04   0.8 3.89e-01   4.0 7.62e-01 8.87e-01f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     400r 5.9033005e+02 2.04e+04 4.63e+04   0.8 9.72e-02   5.3 1.26e-01 5.38e-01f  1
     401r 5.9068107e+02 2.04e+04 7.03e+04   0.8 1.84e-01   5.7 1.33e-01 7.74e-01f  1
     402r 5.9045655e+02 2.04e+04 5.50e+04   0.8 2.06e-01   5.2 7.63e-02 2.29e-01f  1
     403r 5.9009232e+02 2.04e+04 1.16e+05   0.8 1.77e-01   5.7 4.72e-02 1.55e-01f  1
     404r 5.9021669e+02 2.04e+04 4.01e+05   0.8 5.50e-01   5.2 8.67e-02 2.74e-02h  1
     405r 5.9019540e+02 2.04e+04 1.68e+05   0.8 9.57e-02   7.4 2.59e-01 4.90e-02h  1
     406r 5.9018472e+02 2.04e+04 2.34e+05   0.8 9.31e-02   6.9 9.79e-01 2.36e-02f  3
     407r 5.9017664e+02 2.04e+04 1.67e+05   0.8 9.14e-02   6.5 3.82e-01 1.95e-02f  5
     408r 5.9016841e+02 2.04e+04 1.84e+05   0.8 8.94e-02   6.9 1.00e+00 1.96e-02f  5
     409r 5.8995117e+02 2.04e+04 2.10e+05   0.8 8.79e-02   6.4 4.28e-01 6.13e-01h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     410r 5.9018246e+02 2.04e+04 1.94e+05   0.8 2.79e-01   5.9 2.52e-01 8.98e-02h  1
     411r 5.9020635e+02 2.04e+04 1.92e+05   0.8 5.59e-01   6.4 1.64e-02 2.02e-02h  1
     412r 5.8995878e+02 2.04e+04 2.87e+05   0.8 5.85e-01   5.9 5.20e-04 2.58e-02f  1
     413r 5.8991248e+02 2.04e+04 2.93e+05   0.8 6.43e-01   6.3 1.49e-02 8.61e-03h  1
     414r 5.8991284e+02 2.04e+04 2.93e+05   0.8 4.52e+00   5.8 1.20e-04 2.50e-05h  1
     415r 5.8943436e+02 2.04e+04 2.19e+05   0.8 5.02e-01   6.3 1.72e-03 1.76e-01f  1
     416r 5.8900740e+02 2.04e+04 3.24e+05   0.8 1.03e+00   6.7 2.84e-03 7.61e-02f  1
     417r 5.8896372e+02 2.04e+04 3.50e+05   0.8 8.94e-01   6.2 4.01e-02 2.32e-02f  1
     418r 5.8896409e+02 2.04e+04 3.51e+05   0.8 8.38e-01   5.7 2.58e-03 7.13e-04h  1
     419r 5.8889340e+02 2.04e+04 3.85e+05   0.8 8.82e-01   6.2 1.05e-01 3.70e-02f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     420r 5.8888616e+02 2.04e+04 3.86e+05   0.8 8.67e-01   5.7 4.17e-02 3.08e-02h  1
     421r 5.8893599e+02 2.04e+04 3.84e+05   0.8 1.06e+00   5.2 1.47e-02 6.77e-03h  1
     422r 5.8892050e+02 2.04e+04 3.84e+05   0.8 8.51e-01   5.6 2.77e-02 3.46e-02h  1
     423r 5.8897030e+02 2.04e+04 3.79e+05   0.8 7.37e-01   5.2 1.05e-02 1.70e-02h  1
     424r 5.8895179e+02 2.04e+04 3.76e+05   0.8 8.05e-01   5.6 3.54e-01 4.56e-02h  1
     425r 5.8921335e+02 2.04e+04 3.40e+05   0.8 7.40e-01   5.1 3.12e-02 1.28e-01h  1
     426r 5.8918650e+02 2.04e+04 3.35e+05   0.8 6.87e-01   5.5 2.05e-01 4.98e-02h  1
     427r 5.8926164e+02 2.04e+04 2.88e+05   0.8 6.27e-01   5.1 3.40e-01 1.75e-01h  1
     428r 5.8920130e+02 2.04e+04 2.49e+05   0.8 5.37e-01   5.5 6.98e-01 3.12e-01h  1
     429r 5.8855524e+02 2.04e+04 9.73e+04   0.8 3.72e-01   5.0 4.02e-01 7.18e-01h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     430r 5.8687024e+02 2.04e+04 3.67e+04   0.8 1.46e-01   4.5 1.00e+00 1.00e+00f  1
     431r 5.8658585e+02 2.03e+04 6.78e+05   0.8 1.24e-01   4.1 1.43e-01 1.35e-01F  1
     432r 5.8658741e+02 2.03e+04 9.19e+05   0.8 9.28e-02   6.3 3.37e-01 2.52e-02h  1
     433r 5.8658636e+02 2.03e+04 1.48e+06   0.8 6.02e-02   6.7 4.50e-01 2.70e-01f  1
     434r 5.8611174e+02 2.03e+04 9.79e+06   0.8 6.62e-02   8.0 1.89e-02 9.40e-01f  1
     435r 5.8610500e+02 2.03e+04 9.77e+06   0.8 6.66e-01   7.6 1.19e-02 1.77e-03h  1
     436r 5.8607952e+02 2.03e+04 9.70e+06   0.8 6.60e-01   8.0 1.28e-02 6.76e-03h  1
     437r 5.8607811e+02 2.03e+04 9.70e+06   0.8 8.97e-01   7.5 8.56e-03 2.77e-04h  1
     438r 5.8606085e+02 2.03e+04 9.67e+06   0.8 9.59e-01   7.9 7.45e-02 3.09e-03h  1
     439r 5.8595163e+02 2.03e+04 9.48e+06   0.8 9.68e-01   7.5 1.18e-02 1.99e-02h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     440r 5.8594377e+02 2.03e+04 9.46e+06   0.8 9.94e-01   7.0 2.06e-01 1.55e-03h  1
     441r 5.8598745e+02 2.03e+04 8.86e+06   0.8 8.94e-01   6.5 5.97e-02 3.41e-02h  2
     442r 5.8607503e+02 2.03e+04 8.43e+06   0.8 1.01e+00   6.0 5.92e-02 4.94e-02h  1
     443r 5.8605666e+02 2.03e+04 8.38e+06   0.8 9.40e-01   6.5 1.04e-01 6.50e-03h  1
     444r 5.8598864e+02 2.03e+04 7.76e+06   0.8 9.96e-01   6.0 1.09e-01 7.33e-02h  1
     445r 5.8588831e+02 2.03e+04 7.47e+06   0.8 8.77e-01   6.4 2.16e-01 3.80e-02h  1
     446r 5.8539271e+02 2.04e+04 3.47e+06   0.8 8.91e-01   5.9 3.41e-01 5.22e-01h  1
     447r 5.8509480e+02 2.04e+04 1.54e+06   0.8 3.83e-01   6.4 1.00e+00 5.49e-01h  1
     448r 5.8527021e+02 2.03e+04 2.85e+05   0.8 1.80e-01   5.9 5.40e-01 9.42e-01h  1
     449r 5.8518347e+02 2.03e+04 3.30e+03   0.8 1.11e-02   5.4 1.00e+00 1.00e+00f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     450r 5.8484811e+02 2.03e+04 1.42e+03   0.8 1.71e-02   4.9 1.00e+00 1.00e+00f  1
     451r 5.8374084e+02 2.03e+04 1.34e+03   0.8 4.82e-02   4.4 1.00e+00 9.49e-01f  1
     452r 5.8294348e+02 2.03e+04 1.03e+03   0.8 1.11e-01   4.0 1.00e+00 2.59e-01f  1
     453r 5.8066356e+02 2.03e+04 1.03e+03   0.8 3.34e-01   3.5 1.00e+00 4.15e-01f  1
     454r 5.8101411e+02 2.03e+04 1.06e+03   0.8 1.25e-01   3.9 1.00e+00 9.47e-01f  1
     455r 5.7969740e+02 2.02e+04 1.12e+03   0.8 3.76e-01   3.4 1.00e+00 3.70e-01f  1
     456r 5.7909851e+02 2.02e+04 5.56e+04   0.8 1.41e-01   3.9 1.00e+00 1.00e+00f  1
     457r 5.7749638e+02 2.02e+04 3.51e+04   0.8 4.23e-01   3.4 1.00e+00 5.45e-01f  1
     458r 5.7652608e+02 2.01e+04 1.70e+04   0.8 1.59e-01   3.8 1.00e+00 1.00e+00f  1
     459r 5.7456336e+02 2.00e+04 1.14e+04   0.8 4.75e-01   3.3 1.00e+00 4.04e-01f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     460r 5.7060125e+02 2.00e+04 1.00e+04   0.8 4.15e+00   2.9 8.17e-02 8.51e-02f  1
     461r 5.6388777e+02 1.99e+04 4.39e+07   0.8 8.71e-01   3.3 2.33e-02 7.95e-01f  1
     462r 5.6063945e+02 1.98e+04 6.46e+06   0.8 2.05e-01   3.7 8.96e-01 7.94e-01f  1
     463r 5.5977400e+02 1.98e+04 5.97e+06   0.8 1.03e+01   3.2 4.51e-02 5.86e-02f  1
     464r 5.5948269e+02 1.98e+04 5.78e+06   0.8 1.60e+00   3.7 4.03e-02 2.96e-02f  1
     465r 5.5587372e+02 1.98e+04 5.34e+06   0.8 1.29e+00   3.2 1.45e-02 2.27e-01f  1
     466r 5.5290174e+02 1.98e+04 3.27e+06   0.8 8.55e-01   3.6 1.58e-01 5.07e-01f  1
     467r 5.5245118e+02 1.98e+04 2.66e+06   0.8 3.45e-01   4.0 4.74e-01 2.26e-01f  1
     468r 5.5099105e+02 1.97e+04 2.07e+06   0.8 2.85e-01   3.6 3.46e-01 2.92e-01f  1
     469r 5.5039873e+02 1.97e+04 2.02e+06   0.8 1.38e+00   3.1 9.51e-02 2.87e-02f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     470r 5.4995463e+02 1.97e+04 1.84e+06   0.8 6.84e-01   3.5 5.17e-02 9.65e-02f  1
     471r 5.4968221e+02 1.97e+04 1.76e+06   0.8 8.16e-01   3.9 2.56e-02 4.76e-02f  1
     472r 5.4942131e+02 1.97e+04 1.69e+06   0.8 8.34e-01   4.4 2.22e-02 4.26e-02f  1
     473r 5.4930065e+02 1.97e+04 3.44e+06   0.8 1.81e+00   4.8 1.53e-02 2.71e-02f  1
     474r 5.4929969e+02 1.97e+04 3.39e+06   0.8 2.49e-01   6.1 4.37e-02 1.48e-02h  1
     475r 5.4917218e+02 1.97e+04 3.29e+06   0.8 3.94e-01   5.6 4.49e-02 2.79e-02h  1
     476r 5.4896228e+02 1.97e+04 3.18e+06   0.8 3.60e-01   6.1 1.59e-02 3.39e-02h  1
     477r 5.4884028e+02 1.97e+04 3.06e+06   0.8 3.07e-01   6.5 7.02e-02 3.80e-02h  1
     478r 5.4879508e+02 1.97e+04 1.12e+07   0.8 1.78e-01   6.0 9.26e-03 1.18e-01f  1
     479r 5.4877569e+02 1.97e+04 1.11e+07   0.8 1.45e-01   6.4 9.89e-02 5.42e-03h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     480r 5.4785283e+02 1.97e+04 1.05e+07   0.8 1.96e-01   6.0 5.15e-02 7.18e-02h  1
     481r 5.4769502e+02 1.97e+04 1.02e+07   0.8 3.24e-01   6.4 3.65e-02 2.62e-02h  1
     482r 5.4758344e+02 1.97e+04 1.02e+07   0.8 8.29e-01   5.9 4.59e-02 7.05e-03h  1
     483r 5.4736567e+02 1.97e+04 9.93e+06   0.8 8.64e-01   6.3 6.57e-02 2.22e-02h  1
     484r 5.4702727e+02 1.97e+04 9.54e+06   0.8 9.70e-01   6.8 1.90e-01 3.83e-02h  1
     485r 5.4688107e+02 1.97e+04 9.36e+06   0.8 9.21e-01   7.2 4.33e-02 1.89e-02h  1
     486r 5.4674906e+02 1.97e+04 9.18e+06   0.8 8.79e-01   7.6 8.96e-02 1.83e-02h  1
     487r 5.4674548e+02 1.97e+04 9.17e+06   0.8 1.17e+00   7.1 1.01e-02 1.57e-03h  1
     488r 5.4670806e+02 1.97e+04 1.87e+07   0.8 8.12e-01   7.6 6.19e-01 6.25e-03h  1
     489r 5.4666624e+02 1.97e+04 1.84e+07   0.8 6.32e-01   7.1 1.64e-02 1.58e-02h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     490r 5.4691970e+02 1.97e+04 1.83e+07   0.8 4.50e+00   6.6 4.25e-03 5.70e-03h  1
     491r 5.4713011e+02 1.97e+04 1.67e+07   0.8 6.08e-01   6.1 8.64e-02 4.10e-02h  1
     492r 5.4683170e+02 1.97e+04 1.64e+07   0.8 7.28e-01   6.6 1.23e-02 9.72e-02h  1
     493r 5.4664280e+02 1.97e+04 1.25e+07   0.8 6.60e-01   6.1 2.29e-01 4.45e-02h  1
     494r 5.4565144e+02 1.97e+04 1.21e+07   0.8 6.38e-01   6.5 2.34e-02 3.08e-01h  1
     495r 5.4552205e+02 1.97e+04 1.10e+07   0.8 4.46e-01   6.9 9.46e-02 4.88e-02h  1
     496r 5.4547640e+02 1.97e+04 1.09e+07   0.8 9.18e+00   6.5 6.91e-03 1.85e-03h  1
     497r 5.4520385e+02 1.97e+04 9.60e+06   0.8 9.96e-01   6.9 1.30e-01 7.54e-02h  1
     498r 5.4514353e+02 1.97e+04 8.50e+06   0.8 6.32e-01   8.2 3.78e-04 1.83e-02h  1
     499r 5.4491723e+02 1.97e+04 7.70e+06   0.8 8.09e-01   7.7 2.57e-02 5.83e-02h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     500r 5.4475456e+02 1.97e+04 7.55e+06   0.8 9.72e-01   7.3 6.88e-04 6.03e-02h  1
     501r 5.4475318e+02 1.97e+04 7.45e+06   0.8 8.80e-01   6.8 1.04e-02 1.11e-03h  1
     502r 5.4476565e+02 1.97e+04 6.98e+06   0.8 8.88e-01   6.3 4.58e-02 9.93e-03h  1
     503r 5.4475158e+02 1.97e+04 6.87e+06   0.8 9.00e-01   6.7 1.56e-02 1.15e-02h  1
     504r 5.4481625e+02 1.97e+04 6.11e+06   0.8 8.83e-01   6.2 8.67e-02 4.01e-02h  1
     505r 5.4481581e+02 1.97e+04 6.10e+06   0.8 7.69e-01   6.7 4.21e-04 4.61e-04h  1
     506r 5.4486411e+02 1.97e+04 6.18e+06   0.8 8.25e-01   6.2 2.23e-04 2.33e-02h  1
     507r 5.4486385e+02 1.97e+04 5.17e+06   0.8 6.70e-01   6.6 1.71e-01 4.47e-04h  1
     508r 5.4494355e+02 1.97e+04 5.00e+06   0.8 8.19e-01   6.1 2.68e-02 3.27e-02h  1
     509r 5.4494236e+02 1.97e+04 4.99e+06   0.8 8.40e-01   6.6 5.52e-01 1.88e-03h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     510r 5.4495034e+02 1.97e+04 4.97e+06   0.8 8.17e-01   6.1 6.14e-02 2.77e-03h  1
     511r 5.4490953e+02 1.97e+04 4.46e+06   0.8 8.44e-01   6.5 1.80e-01 1.03e-01h  1
     512r 5.4514140e+02 1.97e+04 4.14e+06   0.8 7.29e-01   6.0 1.27e-01 7.20e-02h  1
     513r 5.4512478e+02 1.97e+04 3.63e+06   0.8 7.13e-01   6.5 4.73e-01 1.23e-01h  1
     514r 5.4707546e+02 1.99e+04 1.42e+06   0.8 5.80e-01   6.0 9.90e-03 5.56e-01h  1
     515r 5.4703964e+02 1.97e+04 8.89e+05   0.8 2.64e-01   6.4 8.09e-01 1.00e+00h  1
     516r 5.4716995e+02 1.97e+04 1.27e+06   0.8 3.55e-03   5.9 9.93e-01 1.00e+00f  1
     517r 5.4718921e+02 1.97e+04 2.29e+05   0.8 1.32e-03   7.3 1.00e+00 1.00e+00f  1
     518r 5.4712426e+02 1.97e+04 1.01e+05   0.8 1.80e-03   6.8 1.00e+00 1.00e+00f  1
     519r 5.4708596e+02 1.97e+04 2.06e+04   0.8 1.73e-03   6.3 1.00e+00 1.00e+00f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     520r 5.4708372e+02 1.97e+04 3.61e+03   0.8 1.66e-03   5.8 1.00e+00 1.00e+00f  1
     521r 5.4708306e+02 1.97e+04 1.03e+03   0.8 4.47e-03   5.4 1.00e+00 1.00e+00f  1
     522r 5.4707901e+02 1.97e+04 1.03e+03   0.8 1.34e-02   4.9 1.00e+00 1.00e+00f  1
     523r 5.4679937e+02 1.97e+04 1.03e+03   0.8 4.02e-02   4.4 1.00e+00 1.00e+00f  1
     524r 5.4460810e+02 1.97e+04 1.03e+03   0.8 1.21e-01   3.9 1.00e+00 1.00e+00f  1
     525r 5.4055684e+02 1.96e+04 1.34e+03   0.8 3.62e-01   3.5 1.00e+00 5.11e-01f  1
     526r 5.3759233e+02 1.96e+04 1.03e+03   0.8 1.36e-01   3.9 1.00e+00 1.00e+00f  1
     527r 5.3385802e+02 1.95e+04 1.39e+03   0.8 4.06e-01   3.4 1.00e+00 4.07e-01f  1
     528r 5.3051978e+02 1.95e+04 2.19e+03   0.8 1.52e-01   3.8 1.00e+00 1.00e+00f  1
     529r 5.2607317e+02 1.95e+04 6.53e+03   0.8 4.57e-01   3.4 1.00e+00 4.21e-01f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     530r 5.1949274e+02 1.94e+04 4.15e+03   0.8 1.37e+00   2.9 1.00e+00 1.73e-01f  1
     531r 5.0713656e+02 1.92e+04 2.23e+04   0.8 5.13e-01   3.3 1.00e+00 1.00e+00f  1
     532r 5.0564299e+02 1.92e+04 5.27e+04   0.8 2.41e-02   4.6 1.00e+00 1.00e+00f  1
     533r 5.0379210e+02 1.92e+04 3.00e+04   0.8 7.21e-02   4.2 1.00e+00 1.00e+00f  1
     534r 4.9812617e+02 1.92e+04 1.41e+04   0.8 2.16e-01   3.7 1.00e+00 1.00e+00f  1
     535r 4.9022744e+02 1.91e+04 9.02e+03   0.8 6.48e-01   3.2 1.00e+00 3.99e-01f  1
     536r 4.7547943e+02 1.90e+04 7.49e+03   0.8 1.94e+00   2.7 1.00e+00 1.77e-01f  1
     537r 4.4815366e+02 1.88e+04 1.53e+04   0.8 7.27e-01   3.1 1.00e+00 1.00e+00f  1
     538r 4.3288911e+02 1.87e+04 4.68e+03   0.8 2.73e-01   3.6 1.00e+00 1.00e+00f  1
     539r 4.2799565e+02 1.87e+04 2.23e+03   0.8 1.02e-01   4.0 1.00e+00 1.00e+00f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     540r 4.1416762e+02 1.86e+04 1.02e+03   0.8 3.06e-01   3.5 1.00e+00 1.00e+00f  1
     541r 3.7581962e+02 1.84e+04 1.01e+03   0.8 9.18e-01   3.0 7.49e-01 7.89e-01f  1
     542r 3.5112653e+02 1.83e+04 1.81e+03   0.8 3.56e+00   2.6 1.36e-01 1.12e-01f  1
     543r 3.2662696e+02 1.82e+04 1.22e+03   0.8 1.04e+00   3.0 3.71e-01 4.29e-01f  1
     544r 3.0380420e+02 1.81e+04 1.02e+03   0.8 3.88e-01   3.4 6.70e-01 1.00e+00f  1
     545r 2.7153430e+02 1.79e+04 4.86e+03   0.8 1.17e+00   2.9 5.69e-01 7.03e-01f  1
     546r 2.5866949e+02 1.78e+04 1.61e+03   0.8 4.38e-01   3.4 5.94e-01 1.00e+00f  1
     547r 2.4861176e+02 1.77e+04 3.64e+03   0.8 1.32e+00   2.9 6.49e-01 3.06e-01f  1
     548r 2.4067369e+02 1.76e+04 2.59e+03   0.8 4.94e-01   3.3 1.00e+00 9.29e-01f  1
     549r 2.3845090e+02 1.75e+04 1.26e+03   0.8 1.85e-01   3.7 1.00e+00 1.00e+00f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     550r 2.3463779e+02 1.74e+04 1.02e+03   0.8 5.56e-01   3.3 1.00e+00 9.60e-01f  1
     551r 2.2649785e+02 1.70e+04 1.33e+03   0.8 1.66e+00   2.8 4.65e-01 7.60e-01f  1
     552r 2.2565985e+02 1.70e+04 1.02e+03   0.8 6.23e-01   3.2 5.67e-01 2.80e-01f  1
     553r 2.2480919e+02 1.69e+04 1.02e+03   0.8 2.33e-01   3.6 8.48e-01 1.00e+00f  1
     554r 2.2150603e+02 1.67e+04 1.02e+03   0.8 6.99e-01   3.2 1.00e+00 1.00e+00f  1
     555r 2.1976855e+02 1.67e+04 1.02e+03   0.8 2.62e-01   3.6 1.00e+00 1.00e+00f  1
     556r 2.1589897e+02 1.65e+04 1.02e+03   0.8 7.86e-01   3.1 1.00e+00 1.00e+00f  1
     557r 2.1453058e+02 1.64e+04 1.02e+03   0.8 2.37e+00   2.6 5.38e-01 1.16e-01f  1
     558r 2.1045440e+02 1.62e+04 3.43e+03   0.8 8.85e-01   3.1 1.00e+00 1.00e+00f  1
     559r 2.0838963e+02 1.61e+04 3.37e+03   0.8 2.65e+00   2.6 1.89e-01 8.94e-02f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     560r 2.0251623e+02 1.58e+04 5.60e+03   0.8 9.94e-01   3.0 9.37e-01 1.00e+00f  1
     561r 2.0134335e+02 1.58e+04 2.63e+03   0.8 6.69e-01   3.4 4.87e-01 5.78e-01f  1
     562r 2.0132945e+02 1.58e+04 2.32e+03   0.8 1.37e-01   4.8 2.09e-01 1.18e-01f  1
     563r 2.0127040e+02 1.58e+04 2.42e+03   0.8 2.10e-01   4.3 4.52e-02 2.83e-01f  1
     564r 2.0093206e+02 1.57e+04 1.02e+03   0.8 1.57e-01   3.8 1.00e+00 4.85e-01f  1
     565r 1.9860775e+02 1.56e+04 1.02e+03   0.8 4.72e-01   3.3 1.00e+00 1.00e+00f  1
     566r 1.9773683e+02 1.56e+04 1.02e+03   0.8 1.77e-01   3.8 1.00e+00 1.00e+00f  1
     567r 1.9504075e+02 1.54e+04 1.02e+03   0.8 5.31e-01   3.3 1.00e+00 1.00e+00f  1
     568r 1.9235812e+02 1.53e+04 1.02e+03   0.8 2.86e+00   2.8 8.22e-02 3.47e-01f  1
     569r 1.9195607e+02 1.53e+04 1.02e+03   0.8 5.98e-01   3.2 9.82e-01 1.50e-01f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     570r 1.8572584e+02 1.49e+04 1.51e+03   0.8 1.79e+00   2.8 1.00e+00 6.99e-01f  1
     571r 1.6104726e+02 1.35e+04 2.59e+03   0.8 5.37e+00   2.3 1.00e+00 1.00e+00f  1
     572r 1.5803598e+02 1.35e+04 3.92e+03   0.8 4.51e-02   4.5 5.62e-01 5.79e-01f  1
     573r 1.5586266e+02 1.35e+04 3.98e+03   0.8 1.72e-02   4.9 6.33e-01 1.00e+00f  1
     574r 1.5515725e+02 1.35e+04 7.70e+04   0.8 4.40e-02   4.5 1.00e+00 1.00e+00f  1
     575r 1.5503745e+02 1.35e+04 5.78e+04   0.8 9.28e-02   4.9 1.23e-01 1.21e-01f  1
     576r 1.5478747e+02 1.35e+04 9.78e+05   0.8 5.71e-02   5.3 6.95e-01 6.26e-01f  1
     577r 1.5431114e+02 1.35e+04 5.68e+05   0.8 8.08e-02   4.8 7.02e-01 4.89e-01f  1
     578r 1.5322430e+02 1.35e+04 3.91e+05   0.8 6.86e-02   4.4 5.99e-01 4.23e-01f  1
     579r 1.5254488e+02 1.35e+04 1.72e+05   0.8 2.86e-02   4.8 8.87e-01 1.00e+00f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     580r 1.5094304e+02 1.34e+04 9.27e+04   0.8 1.26e-01   4.3 7.43e-01 1.00e+00f  1
     581r 1.5012516e+02 1.34e+04 4.71e+04   0.8 9.09e-02   4.7 1.00e+00 1.00e+00f  1
     582r 1.4947150e+02 1.34e+04 4.81e+04   0.8 5.16e-01   5.2 5.47e-02 1.94e-01f  1
     583r 1.4942708e+02 1.34e+04 6.48e+04   0.8 1.87e-01   6.5 1.68e-02 1.38e-01h  1
     584r 1.4942321e+02 1.34e+04 7.57e+04   0.8 1.63e-01   6.0 4.60e-01 1.08e-01h  1
     585r 1.4961967e+02 1.34e+04 5.24e+04   0.8 1.50e-01   5.5 3.20e-01 9.40e-01f  1
     586r 1.4998152e+02 1.34e+04 7.77e+03   0.8 1.61e-02   5.1 1.00e+00 8.78e-01f  1
     587r 1.5101058e+02 1.34e+04 1.01e+03   0.8 2.69e-02   4.6 1.00e+00 1.00e+00f  1
     588r 1.5252141e+02 1.34e+04 1.01e+03   0.8 8.07e-02   4.1 1.00e+00 9.91e-01f  1
     589r 1.5152645e+02 1.34e+04 1.01e+03   0.8 2.42e-01   3.6 1.00e+00 8.51e-01f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     590r 1.4952427e+02 1.32e+04 1.02e+03   0.8 7.34e-01   3.1 8.70e-01 1.00e+00f  1
     591r 1.4792016e+02 1.30e+04 1.02e+03   0.8 2.19e+00   2.7 1.00e+00 3.47e-01f  1
     592r 1.4699839e+02 1.28e+04 1.02e+03   0.8 6.53e+00   2.2 1.00e+00 7.14e-02f  1
    Error in an AMPL evaluation. Run with "halt_on_ampl_error yes" to see details.
    Warning: Cutting back alpha due to evaluation error
     593r 1.3578310e+02 2.71e+04 1.02e+03   0.8 1.95e+01   1.7 4.09e-01 2.69e-01f  2
     594r 1.2930650e+02 1.91e+04 1.01e+03   0.8 9.18e-01   3.0 1.00e+00 3.01e-01f  1
     595r 1.2303881e+02 1.66e+04 1.01e+03   0.8 2.75e+00   2.6 5.10e-01 1.95e-01f  1
     596r 1.1934604e+02 9.74e+03 1.01e+03   0.8 8.23e+00   2.1 2.52e-01 6.86e-01f  1
     597r 1.2369514e+02 8.99e+03 1.01e+03   0.8 2.45e+01   1.6 7.06e-01 1.14e-01f  1
     598r 2.0613999e+02 4.77e+03 1.00e+03   0.8 7.33e+01   1.1 4.06e-01 2.67e-01f  1
     599r 2.0837274e+02 4.76e+03 1.27e+03   0.8 2.17e+02   0.7 1.53e-01 2.62e-03f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     600r 2.0916571e+02 4.75e+03 1.18e+03   0.8 7.66e+03    -  3.17e-01 1.67e-03f  1
    Error in an AMPL evaluation. Run with "halt_on_ampl_error yes" to see details.
    Warning: Cutting back alpha due to evaluation error
    Error in an AMPL evaluation. Run with "halt_on_ampl_error yes" to see details.
    Warning: Cutting back alpha due to evaluation error
     601  1.7565396e+02 5.01e+03 1.03e+02  -3.8 1.13e+05    -  1.00e+00 1.25e-01f  4
     602  1.6334003e+02 8.27e+03 8.85e+01  -3.8 1.46e+04    -  1.00e+00 1.74e-01f  3
     603  1.5979687e+02 8.02e+03 8.06e+01  -3.8 1.19e+04    -  1.00e+00 6.07e-02h  5
     604  1.5718164e+02 7.76e+03 7.40e+01  -3.8 1.12e+04    -  1.00e+00 5.75e-02h  5
     605  1.5405058e+02 8.83e+03 4.39e+01  -3.8 1.05e+04    -  1.00e+00 2.16e-01h  3
    Error in an AMPL evaluation. Run with "halt_on_ampl_error yes" to see details.
    Warning: SOC step rejected due to evaluation error
     606  1.5594387e+02 8.53e+03 4.37e+01  -3.8 9.21e+03    -  1.00e+00 8.14e-02h  4
     607  1.5695578e+02 8.42e+03 4.33e+01  -3.8 8.78e+03    -  1.00e+00 1.81e-02h  6
     608  1.5709627e+02 8.41e+03 4.33e+01  -3.8 7.83e+03    -  1.00e+00 1.42e-03h  9
     609  1.5716972e+02 8.40e+03 4.36e+01  -3.8 7.30e+03    -  1.00e+00 4.70e-04h 10
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     610  1.5717441e+02 8.40e+03 4.51e+01  -3.8 7.30e+03    -  1.00e+00 2.15e-05h 14
     611  2.1836460e+02 4.01e+04 2.36e+06  -3.8 7.30e+03    -  1.00e+00 1.44e-01w  1
    Error in an AMPL evaluation. Run with "halt_on_ampl_error yes" to see details.
    Warning: Cutting back alpha due to evaluation error
     612r 1.5717441e+02 8.40e+03 9.99e+02   1.1 0.00e+00    -  0.00e+00 2.14e-09R 26
     613r 1.5277612e+02 4.90e+03 9.84e+02   1.1 1.28e+04    -  4.47e-02 9.92e-04f  1
     614  1.5336892e+02 4.83e+03 1.24e+02  -3.8 8.93e+03    -  1.00e+00 1.69e-02h  6
     615  1.5453155e+02 4.78e+03 1.22e+02  -3.8 8.01e+03    -  1.00e+00 1.72e-02h  6
     616  1.5583526e+02 4.77e+03 1.21e+02  -3.8 7.79e+03    -  1.00e+00 1.10e-02h  6
     617  1.5845671e+02 4.85e+03 1.21e+02  -3.8 7.56e+03    -  1.00e+00 1.37e-02h  5
     618  1.6108665e+02 4.94e+03 1.22e+02  -3.8 7.28e+03    -  1.00e+00 9.59e-03h  5
    Error in an AMPL evaluation. Run with "halt_on_ampl_error yes" to see details.
    Warning: SOC step rejected due to evaluation error
     619  1.6242464e+02 4.97e+03 1.30e+02  -3.8 6.79e+03  -2.2 3.71e-01 1.77e-03h  6
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     620  1.6492728e+02 5.06e+03 1.82e+02  -3.8 9.71e+03    -  1.00e+00 9.36e-03h  5
     621  1.6746141e+02 5.15e+03 3.98e+02  -3.8 9.54e+03    -  5.71e-01 8.65e-03h  5
     622  1.6988643e+02 5.22e+03 1.49e+03  -3.8 1.10e+04    -  1.00e+00 1.24e-02h  5
     623  1.7099322e+02 5.20e+03 3.83e+03  -3.8 1.23e+04    -  1.00e+00 1.13e-02h  6
     624  2.4318744e+02 3.96e+04 4.53e+04  -3.8 1.24e+04    -  1.00e+00 4.35e-01w  1
     625  2.4471916e+02 3.96e+04 7.92e+06  -3.8 1.37e+04  -2.7 2.14e-04 3.89e-02w  1
     626  2.4471660e+02 3.96e+04 7.92e+06  -3.8 2.59e+05  -3.1 7.11e-03 7.50e-06w  1
     627  1.7210272e+02 5.16e+03 8.43e+03  -3.8 1.22e+06  -3.6 1.00e+00 1.36e-02h  5
     628  1.7325525e+02 5.13e+03 1.86e+04  -3.8 1.22e+04    -  1.00e+00 1.38e-02h  6
     629  1.7445808e+02 5.09e+03 4.18e+04  -3.8 1.20e+04    -  1.00e+00 1.35e-02h  6
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     630  1.7574414e+02 5.06e+03 9.67e+04  -3.8 1.12e+04    -  1.00e+00 1.32e-02h  6
     631  1.7853479e+02 5.05e+03 1.25e+05  -3.8 1.88e+05    -  2.01e-01 3.14e-03h  7
     632  1.8129876e+02 5.04e+03 2.15e+05  -3.8 1.95e+05    -  4.83e-01 2.88e-03h  7
     633  1.8436518e+02 5.03e+03 2.44e+05  -3.8 3.51e+05    -  1.15e-01 1.49e-03h  7
     634  1.8453738e+02 4.91e+03 9.04e+04  -3.8 9.72e+03  -4.1 8.39e-01 2.38e-02h  5
     635  1.8489509e+02 4.89e+03 1.99e+05  -3.8 6.30e+04  -4.6 1.00e+00 1.20e-02h  6
     636  1.8542657e+02 4.92e+03 1.62e+06  -3.8 2.28e+05  -5.0 1.00e+00 3.86e-03h  6
     637  3.4306445e+02 4.88e+04 2.30e+08  -3.8 1.39e+05  -5.5 6.97e-01 2.15e-01w  1
    Error in an AMPL evaluation. Run with "halt_on_ampl_error yes" to see details.
    Warning: Cutting back alpha due to evaluation error
     638  1.8807376e+02 4.94e+03 5.67e+06  -3.8 4.56e+04  -2.4 6.97e-01 6.71e-03h  5
     639  1.8838170e+02 4.84e+03 2.64e+06  -3.8 1.54e+04  -2.9 1.00e+00 2.13e-02h  5
    Error in an AMPL evaluation. Run with "halt_on_ampl_error yes" to see details.
    Warning: SOC step rejected due to evaluation error
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     640  1.8887583e+02 4.86e+03 1.14e+07  -3.8 1.15e+05  -3.3 1.00e+00 7.33e-03h  6
     641  1.2886483e+03 1.14e+04 1.14e+07  -3.8 2.29e+07    -  1.49e-03 1.49e-03s 13
     642  1.2868846e+03 1.14e+04 1.14e+07  -3.8 1.24e+05  -3.8 1.85e-03 1.85e-03s 13
     643r 1.2868846e+03 1.14e+04 1.00e+03   0.9 0.00e+00  -4.3 0.00e+00 0.00e+00R  1
     644r 1.1539738e+03 9.06e+03 9.98e+02   0.9 8.35e+03    -  1.10e-02 1.03e-03f  1
     645r 1.0161957e+03 5.87e+03 9.88e+02   0.9 6.14e+03    -  6.81e-03 9.41e-03f  1
     646r 9.9094389e+02 5.77e+03 9.82e+02   0.9 3.58e+03    -  2.64e-02 1.58e-03f  1
     647r 7.1650319e+02 5.05e+03 9.53e+02   0.9 3.41e+02    -  4.84e-02 2.82e-02f  1
     648r 5.2601214e+02 4.08e+03 8.54e+02   0.9 2.45e+00   2.0 2.81e-01 1.06e-01f  1
     649r 4.0249004e+02 3.85e+03 8.13e+02   0.9 2.16e+00   1.5 6.18e-02 4.62e-02f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     650  4.0066710e+02 3.76e+03 1.47e+02  -3.8 1.27e+05    -  3.99e-01 2.32e-02f  4
    Error in an AMPL evaluation. Run with "halt_on_ampl_error yes" to see details.
    Warning: SOC step rejected due to evaluation error
     651  4.0030821e+02 3.74e+03 1.46e+02  -3.8 2.16e+05    -  4.53e-01 5.90e-03h  5
     652  4.0001003e+02 3.73e+03 1.89e+02  -3.8 1.65e+05    -  5.09e-01 3.68e-03h  6
     653  3.9875208e+02 3.68e+03 3.21e+02  -3.8 9.51e+04    -  6.07e-01 1.28e-02f  5
     654  3.9786400e+02 3.65e+03 6.34e+02  -3.8 6.95e+04    -  6.30e-01 8.42e-03f  6
     655  3.8364583e+02 4.66e+04 2.17e+04  -3.8 6.43e+04    -  5.70e-01 2.82e-01f  1
     656  3.8266543e+02 4.62e+04 2.18e+04  -3.8 1.60e+04  -4.8 2.08e-03 1.02e-02h  1
     657  3.8266100e+02 4.62e+04 3.18e+06  -3.8 1.17e+06  -5.2 1.25e-01 1.64e-05h  1
     658r 3.8266100e+02 4.62e+04 1.00e+03   1.5 0.00e+00    -  0.00e+00 3.01e-07R 11
     659r 3.7596893e+02 2.38e+03 1.03e+04   1.5 8.04e+03    -  3.93e-04 3.61e-03f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     660  3.7602153e+02 2.38e+03 2.07e+03  -3.8 3.45e+05    -  3.30e-01 1.60e-04h  1
     661r 3.7602153e+02 2.38e+03 1.00e+03   0.7 0.00e+00    -  0.00e+00 2.86e-07R  8
     662r 3.7342782e+02 1.62e+03 9.95e+02   0.7 1.01e+04    -  5.09e-03 1.68e-03f  1
     663  3.7342908e+02 1.62e+03 1.91e+04  -3.8 5.77e+05    -  4.13e-01 2.17e-05h  1
     664  3.7229945e+02 1.62e+03 2.24e+05  -3.8 1.13e+07    -  9.86e-04 6.27e-05f  1
     665  3.7192835e+02 1.72e+03 2.33e+05  -3.8 4.06e+07    -  2.88e-04 2.12e-04f  3
     666r 3.7192835e+02 1.72e+03 1.00e+03   0.4 0.00e+00    -  0.00e+00 3.30e-07R  8
     667r 3.7330798e+02 9.48e+02 9.91e+02   0.4 9.34e+03    -  1.10e-02 1.52e-03f  1
     668r 3.7330798e+02 9.48e+02 9.99e+02  -0.2 0.00e+00    -  0.00e+00 4.89e-07R  5
     669r 3.6591936e+02 2.23e+02 9.96e+02  -0.2 2.69e+03    -  6.40e-03 1.11e-03f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     670  3.6578799e+02 5.18e+02 1.47e+02  -3.8 4.69e+05    -  1.28e-01 1.38e-03f  1
     671  3.6579241e+02 5.18e+02 4.67e+05  -3.8 1.24e+05    -  1.89e-01 3.65e-05h  1
     672r 3.6579241e+02 5.18e+02 9.99e+02  -0.5 0.00e+00    -  0.00e+00 3.72e-07R 11
     673r 3.5125908e+02 1.84e+02 9.98e+02  -0.5 2.12e+03    -  4.21e-03 9.93e-04f  1
     674r 3.5125908e+02 1.84e+02 9.99e+02  -0.6 0.00e+00    -  0.00e+00 4.44e-07R  4
     675r 3.3956110e+02 2.09e+02 9.97e+02  -0.6 2.67e+03    -  4.36e-03 1.82e-03f  1
     676r 3.3725294e+02 1.01e+03 9.93e+02  -0.6 2.69e+03    -  2.96e-03 4.32e-03f  1
     677r 3.3653726e+02 1.35e+03 9.92e+02  -0.6 5.00e+03    -  1.04e-03 1.29e-03f  1
     678r 3.3799625e+02 1.37e+03 9.88e+02  -0.6 1.87e+03    -  1.16e-02 3.64e-03f  1
     679r 3.4230325e+02 1.33e+03 9.70e+02  -0.6 1.84e+03    -  7.63e-02 1.77e-02f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     680r 3.3124375e+02 1.27e+03 9.25e+02  -0.6 1.80e+03    -  2.21e-01 4.77e-02f  1
     681r 3.1989580e+02 9.46e+02 6.88e+02  -0.6 2.71e+00   0.0 4.82e-01 2.56e-01f  1
     682r 3.3118688e+02 6.09e+02 4.69e+02  -0.6 2.56e+00   0.4 1.43e-02 3.67e-01f  1
     683r 3.3154323e+02 5.41e+02 4.62e+02  -0.6 1.18e+01   0.9 4.54e-02 1.13e-01f  1
     684r 3.3145272e+02 5.32e+02 3.89e+02  -0.6 9.23e-01   2.2 2.95e-01 1.68e-02f  1
     685r 3.3278322e+02 4.66e+02 3.33e+02  -0.6 3.35e-01   1.7 3.86e-01 1.25e-01f  1
     686r 3.4797857e+02 9.36e+02 1.87e+02  -0.6 9.10e-01   1.2 6.75e-01 5.76e-01f  1
     687r 3.5598616e+02 1.48e+03 2.92e+02  -0.6 1.14e+00   1.7 1.00e+00 2.50e-01f  1
     688r 3.5748900e+02 1.27e+03 1.08e+03  -0.6 3.26e-01   1.2 1.00e+00 1.45e-01f  1
     689r 3.6794662e+02 2.21e+02 1.54e+02  -0.6 3.64e-01   0.7 1.00e+00 8.16e-01f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     690r 3.7829610e+02 4.50e+01 2.02e+01  -0.6 1.29e+00   0.2 1.00e+00 1.00e+00f  1
     691r 3.8010555e+02 1.55e+01 1.15e+01  -0.6 7.00e-01   0.6 1.00e+00 1.00e+00f  1
     692r 3.7865578e+02 4.48e+00 4.25e+00  -1.3 3.30e-01   1.1 1.00e+00 1.00e+00f  1
     693r 3.7724778e+02 2.31e+02 1.56e+02  -1.3 2.29e+00   0.6 1.00e+00 1.00e+00f  1
     694r 3.7520189e+02 2.65e+02 7.39e+01  -1.3 2.88e-01   1.9 1.00e+00 6.14e-01f  1
     695r 3.7458037e+02 2.80e+02 3.21e+02  -1.3 1.72e+01   1.4 1.03e-01 7.61e-03f  1
     696r 3.7178433e+02 2.78e+02 3.88e+02  -1.3 4.05e+03    -  1.89e-01 2.61e-02f  1
     697r 3.5983038e+02 8.03e+02 9.25e+02  -1.3 3.84e+03    -  7.29e-01 2.76e-01f  1
     698r 3.7368756e+02 1.98e+03 4.98e+02  -1.3 2.48e+00   1.0 1.00e+00 5.33e-01f  1
     699r 3.7366577e+02 1.95e+03 9.92e+02  -1.3 4.46e-01   2.3 1.00e+00 1.42e-02f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     700r 3.7492049e+02 1.97e+03 1.00e+03  -1.3 4.82e+04    -  4.68e-02 7.44e-03f  1
     701r 3.7475842e+02 1.89e+03 6.19e+02  -1.3 2.71e+03    -  4.55e-01 3.87e-02f  1
     702r 3.7407467e+02 1.40e+03 4.31e+02  -1.3 2.64e+03    -  3.24e-01 2.50e-01f  1
     703r 3.7440981e+02 1.04e+03 4.30e+02  -1.3 2.29e+03    -  6.15e-01 2.52e-01f  1
     704r 3.7700310e+02 3.24e+02 2.26e+03  -1.3 1.80e+03    -  6.20e-01 1.00e+00f  1
     705r 3.7874315e+02 6.04e+01 1.53e+01  -1.3 1.53e+02    -  1.00e+00 1.00e+00h  1
     706r 3.7866215e+02 5.71e-01 2.04e-01  -1.3 5.49e+01    -  1.00e+00 1.00e+00h  1
     707r 3.7655565e+02 9.69e-01 4.44e+00  -3.1 5.28e+01    -  8.38e-01 9.83e-01f  1
     708r 3.7662711e+02 2.25e-01 7.37e-01  -3.1 1.11e-02   1.8 1.00e+00 1.00e+00f  1
     709r 3.7694739e+02 2.25e-01 7.32e-01  -3.1 3.30e-02   1.3 1.00e+00 1.00e+00f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     710r 3.7779429e+02 2.92e-01 1.50e+00  -3.1 9.68e-02   0.9 1.00e+00 1.00e+00f  1
     711r 3.8001688e+02 1.86e+00 1.42e+00  -3.1 2.68e-01   0.4 1.00e+00 1.00e+00f  1
     712r 3.8512632e+02 8.18e+00 3.20e+01  -3.1 6.44e-01  -0.1 7.51e-01 1.00e+00f  1
     713r 3.8593206e+02 7.39e+00 3.42e+01  -3.1 1.57e+00  -0.6 2.48e-01 1.07e-01f  1
     714r 4.0040099e+02 3.53e+01 1.61e+02  -3.1 2.72e+00  -1.0 8.15e-02 1.00e+00f  1
     715r 4.0291413e+02 2.97e+01 1.45e+02  -3.1 3.96e+00  -1.5 1.21e-01 1.65e-01f  1
     716r 4.6115860e+02 1.15e+03 2.84e+02  -3.1 6.53e+00  -2.0 4.17e-02 1.00e+00f  1
     717r 4.5590299e+02 8.17e+02 1.98e+02  -3.1 7.64e+00  -2.5 2.99e-01 3.13e-01h  1
     718r 4.6610529e+02 2.83e+02 1.05e+02  -3.1 1.08e+01  -2.9 8.49e-01 1.00e+00f  1
     719r 4.5424124e+02 5.69e+01 1.67e+01  -3.1 1.68e+01  -3.4 9.91e-01 1.00e+00f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     720r 4.5232505e+02 1.29e+00 1.80e-01  -3.1 4.05e-02  -0.3 1.00e+00 1.00e+00h  1
     721r 4.5152396e+02 2.25e-01 4.09e-02  -3.1 1.07e-01  -0.8 1.00e+00 1.00e+00h  1
     722r 4.5132968e+02 2.25e-01 1.63e+01  -4.6 3.60e-02  -0.3 9.57e-01 8.62e-01f  1
     723r 4.6969038e+02 1.08e+02 8.69e+01  -4.6 1.27e+00  -0.8 1.55e-01 1.00e+00f  1
     724r 4.9510933e+02 3.79e+00 8.60e+00  -4.6 2.85e+00  -1.3 1.00e+00 1.00e+00f  1
     725r 5.0986947e+02 1.66e+00 1.31e+01  -4.6 4.90e+00  -1.8 1.00e+00 7.40e-01f  1
     726r 5.1387670e+02 1.14e+00 1.36e+02  -4.6 6.88e+00  -2.2 1.00e+00 3.31e-01f  1
     727r 5.2155721e+02 1.52e+00 4.51e+01  -4.6 1.12e+01  -2.7 1.00e+00 6.92e-01f  1
     728r 5.2923607e+02 1.44e+00 3.17e+01  -4.6 1.57e+01  -3.2 1.00e+00 1.00e+00f  1
     729r 5.3061167e+02 4.40e-01 1.72e+01  -4.6 3.50e+01  -3.7 1.00e+00 1.00e+00f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     730r 5.3133342e+02 9.20e-01 1.13e+02  -4.6 1.06e+02  -4.2 3.48e-01 1.00e+00f  1
     731r 5.3555963e+02 8.22e+00 2.12e+02  -4.6 3.37e+02  -4.6 1.65e-01 1.00e+00f  1
     732r 5.3787027e+02 8.31e+00 1.89e+02  -4.6 9.28e+02  -5.1 5.07e-01 1.30e-01f  1
     733r 5.3892449e+02 8.30e+00 1.88e+02  -4.6 3.18e+03  -5.6 8.37e-02 1.41e-02f  1
     734r 5.4059004e+02 8.21e+00 2.26e+02  -4.6 1.18e+03  -5.2 1.00e+00 8.87e-02f  1
     735r 5.5545430e+02 2.74e+01 1.90e+02  -4.6 3.57e+03  -5.6 1.08e-01 1.59e-01f  1
     736r 5.6925287e+02 6.14e+01 2.29e+02  -4.6 1.41e+03  -5.2 1.97e-02 7.48e-01f  1
     737r 5.6925415e+02 6.14e+01 2.34e+02  -4.6 4.05e+03  -5.7 5.60e-02 1.35e-05f  1
     738r 5.8199594e+02 5.99e+01 1.67e+02  -4.6 1.51e+03  -5.3 5.18e-02 2.96e-01f  1
     739r 5.8525522e+02 6.47e+01 2.14e+02  -4.6 4.57e+03  -5.7 2.34e-02 8.15e-02f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     740r 5.8486147e+02 6.21e+01 7.55e+02  -4.6 1.70e+03  -5.3 1.00e+00 5.01e-02f  1
     741r 5.8403551e+02 6.76e+01 7.03e+02  -4.6 5.16e+03  -5.8 8.74e-02 7.20e-02f  1
     742r 5.7284823e+02 6.60e+01 6.23e+02  -4.6 1.92e+03  -5.4 1.00e+00 2.25e-01f  1
     743r 5.7284109e+02 6.60e+01 6.36e+02  -4.6 5.83e+03  -5.8 6.91e-03 1.43e-04f  1
     744r 5.7251043e+02 6.52e+01 1.13e+03  -4.6 2.16e+03  -5.4 6.00e-01 1.39e-02f  1
     745r 5.7239154e+02 6.72e+01 1.28e+03  -4.6 6.58e+03  -5.9 1.83e-01 4.30e-02f  1
     746r 5.6838946e+02 6.49e+01 1.10e+03  -4.6 2.44e+03  -5.5 5.29e-01 1.03e-01f  1
     747r 5.6410908e+02 7.71e+01 8.92e+02  -4.6 7.44e+03  -5.9 9.99e-06 9.34e-02f  1
     748r 5.6310714e+02 7.60e+01 9.10e+02  -4.6 2.75e+03  -5.5 1.76e-01 1.57e-02f  1
     749r 5.2541783e+02 8.14e+02 7.05e+02  -4.6 8.43e+03  -6.0 2.86e-02 2.25e-01f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     750r 5.2541680e+02 8.14e+02 7.05e+02  -4.6 2.97e+03  -5.6 2.88e-01 2.03e-05f  1
     751r 5.2371810e+02 1.59e+03 5.95e+02  -4.6 1.11e+03  -5.1 6.54e-01 1.57e-01f  1
     752r 5.2341362e+02 1.51e+03 5.52e+02  -4.6 3.35e+03  -5.6 1.98e-05 7.15e-02f  1
     753r 5.2340599e+02 1.50e+03 5.59e+02  -4.6 1.25e+03  -5.2 4.00e-01 5.96e-03f  1
     754r 5.4526760e+02 2.89e+04 1.53e+03  -4.6 1.99e+04  -5.7 2.08e-06 3.79e-02f  1
     755r 5.5000928e+02 2.15e+04 9.80e+02  -4.6 4.18e+03  -5.2 1.01e-01 2.88e-01f  1
     756r 5.5268210e+02 1.98e+04 1.47e+03  -4.6 7.51e+02  -3.9 9.98e-03 8.09e-02H  1
     757r 5.5317319e+02 1.95e+04 1.42e+03  -4.6 3.68e+03  -4.4 1.35e-02 2.18e-02h  1
     758r 5.5205856e+02 1.91e+04 1.38e+03  -4.6 6.93e+03  -4.9 3.30e-02 2.27e-02f  1
     759r 5.5205856e+02 1.91e+04 1.38e+03  -4.6 5.29e+01   3.0 0.00e+00 4.72e-07R  5
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     760r 5.5208860e+02 1.79e+04 9.02e+02  -4.6 2.55e+00   2.5 3.36e-01 3.04e-01f  1
     761r 5.5089972e+02 1.49e+04 6.53e+02  -4.6 5.98e+00   2.0 2.75e-01 3.21e-01f  1
     762r 5.4827421e+02 6.42e+03 4.47e+02  -4.6 9.90e+00   1.5 1.77e-01 5.54e-01f  1
     763r 5.4763810e+02 5.36e+03 3.98e+02  -4.6 5.91e+00   1.1 8.04e-02 1.17e-01f  1
     764r 5.4592847e+02 4.45e+03 3.50e+02  -4.6 4.97e+00   0.6 1.34e-01 2.12e-01f  1
     765r 5.3940493e+02 2.46e+04 9.98e+02  -4.6 4.83e+01   0.1 6.19e-02 6.12e-02F  1
     766r 5.3471639e+02 1.43e+04 6.25e+02  -4.6 2.90e+00  -0.4 3.06e-01 4.07e-01f  1
     767r 5.3471604e+02 1.43e+04 6.25e+02  -4.6 1.14e+00  -0.8 4.49e-02 6.19e-05f  1
     768r 5.3438130e+02 1.20e+04 5.37e+02  -4.6 1.12e+00  -1.3 3.65e-01 1.60e-01f  1
     769r 5.3343181e+02 8.03e+03 3.84e+02  -4.6 8.39e-01  -1.8 1.74e-01 3.25e-01f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     770r 5.3343140e+02 8.03e+03 4.57e+02  -4.6 1.84e+00  -2.3 3.62e-01 7.94e-05f  1
     771r 5.3173650e+02 5.89e+03 3.59e+02  -4.6 5.45e+00  -2.8 3.15e-01 2.85e-01f  1
     772r 5.2821619e+02 3.41e+03 1.33e+02  -4.6 1.40e+01  -3.2 5.78e-01 5.93e-01f  1
     773r 5.2754316e+02 3.04e+03 2.40e+02  -4.6 5.62e+01  -3.7 1.39e-02 1.35e-01f  1
     774r 5.2843982e+02 6.00e+02 1.97e+02  -4.6 1.40e+02  -4.2 2.75e-01 7.93e-01f  1
     775r 5.2863353e+02 4.93e+02 2.36e+02  -4.6 4.07e+02  -4.7 6.83e-01 1.78e-01f  1
     776r 5.2925446e+02 4.05e+02 1.98e+02  -4.6 1.17e+03  -5.1 5.30e-02 1.78e-01f  1
     777r 5.2926626e+02 4.05e+02 3.76e+02  -4.6 3.32e+03  -5.6 1.29e-01 2.31e-04f  1
     778r 5.2932713e+02 4.04e+02 3.92e+02  -4.6 1.32e+04  -6.1 6.90e-03 2.19e-03f  1
     779r 5.3184233e+02 3.86e+02 8.80e+02  -4.6 3.74e+03  -5.7 2.82e-01 4.44e-02f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     780r 5.3208249e+02 3.86e+02 8.79e+02  -4.6 2.27e+05  -6.1 4.82e-06 2.89e-04f  1
     781r 5.3281422e+02 3.71e+02 8.43e+02  -4.6 5.05e+03  -5.7 1.34e-02 4.10e-02f  1
     782r 5.3281929e+02 3.71e+02 8.43e+02  -4.6 1.57e+03  -5.3 1.45e-01 2.57e-04f  1
     783r 5.3218253e+02 3.59e+02 8.15e+02  -4.6 4.75e+03  -5.8 2.49e-01 3.29e-02f  1
     784r 5.2315118e+02 1.48e+03 7.55e+02  -4.6 2.16e+04  -6.2 1.65e-02 7.45e-02f  1
     785r 5.2533534e+02 1.46e+03 7.25e+02  -4.6 6.83e+03  -5.8 6.44e-02 3.87e-02f  1
     786r 5.2467210e+02 1.46e+03 7.24e+02  -4.6 1.69e+04  -6.3 6.66e-02 1.33e-03f  1
     787r 5.1793885e+02 1.55e+03 9.01e+02  -4.6 6.09e+03  -5.9 4.15e-01 5.79e-02f  1
     788r 5.1837998e+02 1.58e+03 9.00e+02  -4.6 8.53e+03  -3.6 3.68e-04 3.68e-04s 14
     789r 5.1837998e+02 1.58e+03 9.00e+02  -4.6 6.85e+02  -3.2 0.00e+00 0.00e+00R  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     790r 5.1989780e+02 1.10e+03 5.98e+02  -4.6 9.68e-01   1.7 4.09e-03 3.17e-01f  1
     791r 5.2000789e+02 1.08e+03 5.92e+02  -4.6 7.02e-01   1.2 2.29e-01 1.16e-02f  1
     792r 5.2258495e+02 8.99e+02 4.89e+02  -4.6 7.00e-01   0.8 2.02e-01 1.73e-01f  1
     793r 5.2308706e+02 8.67e+02 9.01e+02  -4.6 5.80e-01   0.3 4.66e-01 3.48e-02f  1
     794r 5.2589886e+02 7.20e+02 6.85e+02  -4.6 5.59e-01  -0.2 1.08e-01 1.78e-01f  1
     795r 5.2714789e+02 6.47e+02 9.21e+02  -4.6 4.59e-01  -0.7 6.81e-01 1.00e-01f  1
     796r 5.2944209e+02 5.08e+02 5.55e+02  -4.6 4.11e-01  -1.1 5.58e-02 2.16e-01f  1
     797r 5.3332542e+02 2.56e+02 3.75e+02  -4.6 5.33e-01  -1.6 8.20e-01 4.94e-01f  1
     798r 5.3494769e+02 1.23e+02 1.73e+02  -4.6 1.85e+00  -2.1 4.96e-01 5.18e-01f  1
     799r 5.3501733e+02 1.08e+02 1.50e+02  -4.6 4.39e+00  -2.6 2.72e-01 1.27e-01f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     800r 5.3490514e+02 9.27e+01 1.24e+02  -4.6 8.89e+00  -3.0 2.02e-01 1.41e-01f  1
     801r 5.3454870e+02 6.13e+01 2.03e+02  -4.6 2.66e+01  -3.5 1.71e-01 3.44e-01f  1
     802r 5.3429646e+02 5.10e+01 1.18e+02  -4.6 7.96e+01  -4.0 4.54e-01 1.69e-01f  1
     803r 5.3049120e+02 7.95e+00 2.24e+02  -4.6 2.39e+02  -4.5 4.45e-01 8.66e-01f  1
     804r 5.2922361e+02 6.65e+00 1.30e+02  -4.6 7.19e+02  -5.0 3.48e-01 1.66e-01f  1
     805r 5.2878727e+02 6.34e+00 7.57e+01  -4.6 2.17e+03  -5.4 3.01e-01 4.86e-02f  1
     806r 5.2754081e+02 6.22e+00 1.41e+02  -4.6 6.62e+03  -5.9 2.63e-01 2.92e-02f  1
     807r 5.2756760e+02 6.20e+00 6.29e+02  -4.6 1.42e+03  -5.5 3.72e-01 2.44e-03f  1
     808r 5.2719046e+02 1.82e+01 5.07e+02  -4.6 4.95e+03  -6.0 9.64e-06 4.91e-02f  1
     809r 5.2719263e+02 1.90e+01 8.56e+02  -4.6 2.37e+03  -5.5 7.51e-01 4.05e-02f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     810r 5.2761210e+02 2.48e+01 8.41e+02  -4.6 7.71e+03  -6.0 2.37e-02 1.90e-02f  1
     811r 5.2827497e+02 2.73e+01 8.11e+02  -4.6 2.69e+03  -5.6 7.26e-02 4.02e-02f  1
     812r 5.2827425e+02 3.03e+01 8.13e+02  -4.6 9.44e+03  -6.1 5.95e-02 7.24e-03f  1
     813r 5.3812332e+02 1.19e+03 1.01e+03  -4.6 2.46e+03  -5.6 4.32e-01 1.75e-01f  1
     814r 5.4182336e+02 1.17e+03 1.02e+03  -4.6 7.32e+03  -6.1 5.97e-02 4.00e-02f  1
     815r 5.4715072e+02 1.09e+03 9.81e+02  -4.6 2.27e+03  -5.7 1.18e-01 9.04e-02f  1
     816r 5.4495049e+02 1.07e+03 9.20e+02  -4.6 6.80e+03  -6.2 4.56e-06 2.75e-02f  1
     817r 5.4072232e+02 1.04e+03 8.50e+02  -4.6 2.55e+03  -5.7 6.55e-01 7.66e-02f  1
     818r 5.3603458e+02 1.01e+03 8.23e+02  -4.6 7.61e+03  -6.2 1.59e-03 3.14e-02f  1
     819r 5.2166018e+02 8.45e+02 5.28e+02  -4.6 2.87e+03  -5.8 1.91e-01 3.58e-01f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     820r 5.2163271e+02 8.45e+02 5.28e+02  -4.6 8.60e+03  -6.3 8.23e-02 2.05e-04f  1
     821r 5.1564342e+02 6.25e+02 3.67e+02  -4.6 3.24e+03  -5.8 3.33e-01 3.05e-01f  1
     822r 5.1559545e+02 6.25e+02 3.67e+02  -4.6 9.72e+03  -6.3 1.04e-01 3.41e-04f  1
     823r 5.1218496e+02 4.80e+02 2.82e+02  -4.6 3.65e+03  -5.9 6.36e-03 2.32e-01f  1
     824r 5.1139365e+02 4.77e+02 2.80e+02  -4.6 1.33e+04  -6.4 1.30e-01 6.81e-03f  1
     825r 5.1108191e+02 4.52e+02 2.65e+02  -4.6 4.13e+03  -5.9 3.17e-01 5.32e-02f  1
     826r 5.0833182e+02 4.39e+02 2.58e+02  -4.6 2.11e+04  -6.4 1.89e-05 2.77e-02f  1
     827r 5.0834721e+02 4.38e+02 2.57e+02  -4.6 4.66e+03  -6.0 3.60e-02 2.83e-03f  1
     828r 5.0843090e+02 4.08e+02 2.39e+02  -4.6 1.74e+03  -5.6 3.31e-01 6.82e-02f  1
     829r 5.2180210e+02 4.46e+02 1.26e+02  -4.6 5.25e+03  -6.0 1.17e-01 5.81e-01f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     830r 5.2625437e+02 3.33e+02 9.93e+01  -4.6 1.98e+03  -5.6 2.36e-01 2.66e-01f  1
     831r 5.2661029e+02 3.31e+02 7.81e+01  -4.6 9.56e+03  -6.1 1.13e-01 6.84e-03f  1
     832r 5.3117418e+02 2.69e+02 8.37e+01  -4.6 2.23e+03  -5.7 3.20e-01 2.05e-01f  1
     833r 5.3173351e+02 2.68e+02 8.85e+01  -4.6 1.09e+04  -6.1 3.67e-02 8.01e-03f  1
     834r 5.3356364e+02 2.55e+02 8.43e+01  -4.6 2.50e+03  -5.7 4.77e-02 4.78e-02f  1
     835r 5.4158527e+02 3.28e+02 6.99e+01  -4.6 1.47e+04  -6.2 1.79e-02 6.67e-02f  1
     836r 5.4490447e+02 3.07e+02 6.69e+01  -4.6 2.83e+03  -5.8 8.07e-02 6.96e-02f  1
     837r 5.4545895e+02 3.06e+02 8.91e+01  -4.6 1.54e+04  -6.2 6.43e-02 3.55e-03f  1
     838r 5.5005488e+02 2.84e+02 7.80e+01  -4.6 3.19e+03  -5.8 3.68e-02 8.13e-02f  1
     839r 5.5050625e+02 2.83e+02 1.22e+02  -4.6 2.37e+04  -6.3 3.54e-02 2.35e-03f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     840r 5.5239097e+02 2.76e+02 1.49e+02  -4.6 3.59e+03  -5.9 4.16e-01 2.82e-02f  1
     841r 5.5549891e+02 2.84e+02 1.64e+02  -4.6 4.87e+04  -6.4 2.04e-02 1.17e-02f  1
     842r 5.5586121e+02 2.83e+02 1.68e+02  -4.6 3.89e+03  -5.9 6.30e-02 2.81e-03f  1
     843r 5.6058996e+02 2.67e+02 1.50e+02  -4.6 1.45e+03  -5.5 1.25e-01 5.86e-02f  1
     844r 5.7132731e+02 2.54e+02 1.40e+02  -4.6 4.50e+03  -6.0 8.14e-02 7.21e-02f  1
     845r 5.7328590e+02 2.48e+02 1.41e+02  -4.6 1.63e+03  -5.6 2.60e-01 2.33e-02f  1
     846r 6.5142000e+02 4.74e+02 5.86e+01  -4.6 5.23e+03  -6.0 1.70e-01 4.33e-01f  1
     847r 6.5634067e+02 4.44e+02 5.61e+01  -4.6 1.86e+03  -5.6 1.18e-01 6.48e-02f  1
     848r 6.7329810e+02 4.30e+02 6.93e+01  -4.6 6.40e+03  -6.1 2.40e-01 9.51e-02f  1
     849r 6.7657740e+02 4.13e+02 6.90e+01  -4.6 2.10e+03  -5.7 1.75e-01 4.10e-02f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     850r 6.8056111e+02 4.06e+02 2.75e+02  -4.6 7.53e+03  -6.1 4.95e-01 2.00e-02f  1
     851r 6.8658252e+02 3.79e+02 2.45e+02  -4.6 2.37e+03  -5.7 1.33e-01 7.00e-02f  1
     852r 6.9454285e+02 3.74e+02 2.58e+02  -4.6 9.10e+03  -6.2 9.76e-02 3.54e-02f  1
     853r 6.9495328e+02 3.73e+02 2.04e+02  -4.6 2.60e+03  -5.8 4.27e-01 4.55e-03f  1
     854r 6.9996162e+02 3.69e+02 1.91e+02  -4.6 1.14e+04  -6.2 4.37e-03 2.01e-02f  1
     855r 7.0004717e+02 3.69e+02 1.87e+02  -4.6 2.92e+03  -5.8 9.97e-02 6.97e-04f  1
     856r 7.0686949e+02 3.70e+02 1.76e+02  -4.6 1.48e+04  -6.3 1.58e-02 2.51e-02f  1
     857r 7.1379657e+02 3.61e+02 1.64e+02  -4.6 3.55e+03  -5.9 2.72e-01 6.09e-02f  1
     858r 7.2388877e+02 4.02e+02 2.30e+02  -4.6 1.89e+04  -6.3 1.15e-01 3.90e-02f  1
     859r 7.4338009e+02 4.40e+02 2.31e+02  -4.6 3.97e+03  -5.9 6.06e-01 2.49e-01f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     860r 7.5050669e+02 8.86e+02 2.24e+02  -4.6 3.01e+04  -6.4 9.75e-03 2.26e-02f  1
     861r 7.5051868e+02 8.86e+02 2.89e+02  -4.6 4.37e+03  -6.0 5.35e-01 2.18e-04f  1
     862r 7.5831561e+02 9.03e+02 2.79e+02  -4.6 2.39e+04  -6.4 2.14e-02 3.17e-02f  1
     863r 7.6475020e+02 8.85e+02 2.70e+02  -4.6 5.19e+03  -6.0 3.95e-01 9.36e-02f  1
     864r 7.6679857e+02 8.81e+02 2.69e+02  -4.6 4.53e+04  -6.5 5.50e-03 4.87e-03f  1
     865r 7.6738147e+02 8.75e+02 2.67e+02  -4.6 6.56e+03  -6.1 1.33e-02 6.72e-03f  1
     866r 7.6740200e+02 8.75e+02 7.21e+02  -4.6 2.01e+03  -5.6 1.00e+00 7.88e-04f  1
     867r 7.7217881e+02 8.31e+02 6.90e+02  -4.6 7.96e+03  -6.1 1.01e-01 4.97e-02f  1
     868r 7.7338521e+02 7.97e+02 1.47e+03  -4.6 2.45e+03  -5.7 1.74e-01 4.17e-02f  1
     869r 7.7718118e+02 9.67e+02 1.42e+03  -4.6 2.15e+04  -6.2 4.48e-02 1.73e-02f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     870r 7.7706385e+02 9.68e+02 1.55e+03  -4.6 2.70e+01  -0.3 1.21e-03 3.22e-04h  2
     871r 7.7722124e+02 9.67e+02 1.55e+03  -4.6 4.91e+00  -0.8 1.93e-03 5.13e-03h  1
     872r 7.7722048e+02 9.67e+02 1.55e+03  -4.6 3.87e-01   1.4 8.11e-05 1.10e-04h  1
     873r 7.7719829e+02 9.67e+02 1.55e+03  -4.6 7.44e+00   1.0 7.75e-08 2.14e-04h  1
     874r 7.7721107e+02 9.67e+02 1.55e+03  -4.6 3.99e+01   1.4 2.74e-04 2.42e-05f  1
     875r 7.7721102e+02 9.67e+02 1.87e+03  -4.6 3.55e-01   1.8 4.38e-04 1.20e-05h  1
     876r 7.7720107e+02 9.66e+02 1.87e+03  -4.6 4.19e+00   1.3 7.64e-04 1.64e-04h  1
     877r 7.7712409e+02 9.66e+02 2.29e+03  -4.6 5.27e+00   0.9 1.84e-04 1.03e-03h  1
     878r 7.7712400e+02 9.66e+02 2.03e+03  -4.6 9.15e-02   2.2 1.10e-03 1.61e-05h  1
     879r 7.7713241e+02 9.66e+02 1.88e+03  -4.6 7.24e+00   1.7 2.59e-04 9.88e-05f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     880r 7.7719672e+02 9.66e+02 2.17e+03  -4.6 4.01e+01   1.2 8.75e-05 1.24e-04f  1
     881r 7.7719672e+02 9.66e+02 2.17e+03  -4.6 9.13e-02   2.6 0.00e+00 3.56e-07R  4
     882r 7.7618864e+02 8.39e+02 1.88e+03  -4.6 5.80e-01   2.1 2.91e-03 1.41e-01f  1
     883r 7.7607145e+02 8.09e+02 1.81e+03  -4.6 4.58e-01   2.5 3.54e-02 4.22e-02f  1
     884r 7.7429685e+02 5.24e+02 1.26e+03  -4.6 4.96e-01   2.0 2.07e-01 3.67e-01f  1
     885r 7.7390105e+02 4.46e+02 1.01e+03  -4.6 3.13e-01   2.5 3.77e-01 1.60e-01f  1
     886r 7.7015050e+02 6.63e+01 3.79e+02  -4.6 2.92e-01   2.0 6.35e-01 8.38e-01f  1
     887r 7.6844841e+02 5.20e+01 2.83e+02  -4.6 5.41e-02   1.5 1.33e-01 2.66e-01f  1
     888r 7.6818970e+02 4.67e+01 3.37e+02  -4.6 3.07e-02   1.9 7.04e-02 1.01e-01f  1
     889r 7.6397306e+02 3.35e+01 6.02e+02  -4.6 4.76e-02   1.5 1.30e-01 6.98e-01f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     890r 7.5980044e+02 2.81e+01 3.54e+02  -4.6 5.37e-02   1.0 4.17e-01 4.27e-01f  1
     891r 7.5982748e+02 2.62e+01 2.70e+02  -4.6 7.42e-02   0.5 3.71e-01 7.20e-02f  1
     892r 7.5987316e+02 2.35e+01 4.62e+02  -4.6 4.44e-02   0.0 5.69e-01 1.22e-01f  1
     893r 7.6048201e+02 2.04e+01 2.08e+02  -4.6 8.14e-02  -0.5 4.64e-02 5.75e-01f  1
     894r 7.6247725e+02 4.72e+00 3.33e+02  -4.6 2.32e-01  -0.9 3.67e-01 7.85e-01f  1
     895r 7.6590473e+02 5.98e+00 1.46e+02  -4.6 5.91e-01  -1.4 5.80e-01 5.63e-01f  1
     896r 7.6830242e+02 4.68e+00 4.14e+02  -4.6 1.41e+00  -1.9 8.85e-01 2.17e-01f  1
     897r 7.7541625e+02 2.64e+00 3.78e+02  -4.6 4.18e+00  -2.4 9.91e-01 4.36e-01f  1
     898r 7.7533213e+02 2.47e+00 4.87e+02  -4.6 1.02e+01  -2.8 1.00e+00 6.40e-02f  1
     899r 7.7533475e+02 2.04e+00 1.37e+03  -4.6 2.87e+01  -3.3 9.82e-01 1.72e-01f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     900r 7.7533950e+02 1.99e+00 1.32e+03  -4.6 8.19e+01  -3.8 7.69e-01 2.46e-02f  1
     901r 7.7543096e+02 1.91e+00 1.27e+03  -4.6 1.06e+02  -4.3 1.47e-01 4.15e-02f  1
     902r 7.7636704e+02 1.68e+00 1.13e+03  -4.6 3.20e+02  -4.8 6.09e-01 1.23e-01f  1
     903r 7.7642430e+02 1.67e+00 1.12e+03  -4.6 9.62e+02  -5.2 5.87e-02 3.19e-03f  1
     904r 7.7704512e+02 2.53e+00 1.05e+03  -4.6 2.95e+03  -5.7 1.59e-01 1.53e-02f  1
     905r 7.7766724e+02 4.59e+00 9.87e+02  -4.6 8.73e+03  -6.2 1.51e-02 6.30e-03f  1
     906r 7.8180821e+02 1.88e+02 8.96e+02  -4.6 5.05e+04  -6.7 3.14e-02 1.54e-02f  1
     907r 7.8180848e+02 1.88e+02 1.41e+03  -4.6 9.58e+03  -6.2 1.22e-01 3.32e-06f  1
     908r 7.8555345e+02 4.31e+02 1.04e+03  -4.6 5.92e+04  -6.7 3.43e-01 1.37e-02f  1
     909r 7.8569515e+02 4.31e+02 1.40e+03  -4.6 1.09e+05  -7.2 2.38e-01 5.20e-04f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     910r 7.8744183e+02 4.29e+02 1.61e+03  -4.6 2.29e+05  -7.7 7.73e-01 4.02e-03f  1
     911r 8.0772592e+02 8.80e+02 1.64e+03  -4.6 7.67e+04    -  4.19e-01 2.91e-02f  1
     912r 8.0772647e+02 8.80e+02 1.67e+03  -4.6 7.50e+04    -  2.82e-01 2.08e-06f  1
     913r 8.0772665e+02 8.80e+02 1.74e+03  -4.6 1.07e+04    -  1.00e+00 6.04e-06h  1
     914r 8.1049210e+02 4.84e+02 9.21e+02  -4.6 8.38e+02    -  1.00e+00 4.69e-01f  1
     915r 8.1048189e+02 4.84e+02 9.20e+02  -4.6 5.82e+00   0.2 1.27e-03 8.98e-04h  1
     916r 8.1049419e+02 4.81e+02 1.59e+03  -4.6 7.36e+00  -0.3 1.81e-04 4.96e-03h  1
     917r 8.1054143e+02 4.73e+02 2.86e+03  -4.6 8.11e+00  -0.8 4.35e-03 1.70e-02h  1
     918r 8.1070979e+02 4.49e+02 3.42e+03  -4.6 4.58e+00  -1.3 6.11e-03 5.24e-02h  1
     919r 8.1108329e+02 4.00e+02 3.04e+03  -4.6 1.75e+00  -1.7 5.69e-02 1.10e-01h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     920r 8.1162687e+02 3.18e+02 2.40e+03  -4.6 2.30e+00  -2.2 1.13e-01 2.10e-01h  1
     921r 8.1235075e+02 2.51e+02 1.88e+03  -4.6 5.56e-01  -0.9 2.04e-02 2.19e-01h  1
     922r 8.1237455e+02 2.49e+02 1.86e+03  -4.6 1.76e-01   2.3 7.89e-03 7.97e-03h  1
     923r 8.1247802e+02 2.40e+02 2.01e+03  -4.6 4.47e-02   2.7 2.94e-03 3.62e-02h  1
     924r 8.1279375e+02 2.13e+02 2.10e+03  -4.6 4.42e-02   2.2 3.67e-02 1.13e-01h  1
     925r 8.1298837e+02 1.97e+02 1.99e+03  -4.6 4.15e-02   2.6 7.29e-02 7.77e-02h  1
     926r 8.1315429e+02 1.83e+02 1.84e+03  -4.6 3.92e-02   2.2 4.34e-02 7.10e-02h  1
     927r 8.1331561e+02 1.60e+02 1.60e+03  -4.6 2.62e-02   2.6 1.63e-02 1.28e-01h  1
     928r 8.1338585e+02 1.47e+02 1.47e+03  -4.6 1.84e-02   2.1 1.11e-01 7.97e-02h  1
     929r 8.1343290e+02 1.33e+02 1.33e+03  -4.6 1.70e-02   2.5 1.88e-01 9.75e-02h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     930r 8.1346407e+02 1.29e+02 1.48e+03  -4.6 2.93e-02   2.1 1.27e-01 2.45e-02h  1
     931r 8.1348351e+02 1.28e+02 1.39e+03  -4.6 2.74e-02   2.5 1.98e-01 1.43e-02h  1
     932r 8.1348626e+02 1.27e+02 1.75e+03  -4.6 2.71e-02   2.0 1.06e-01 1.69e-03h  1
     933r 8.1350672e+02 1.26e+02 1.44e+03  -4.6 6.76e-02   1.5 3.70e-02 8.63e-03h  1
     934r 8.1369483e+02 1.10e+02 1.13e+03  -4.6 2.69e-02   2.0 2.98e-02 1.35e-01h  1
     935r 8.1375246e+02 1.07e+02 1.84e+03  -4.6 6.75e-02   1.5 2.59e-03 2.61e-02h  1
     936r 8.1375858e+02 1.06e+02 2.14e+03  -4.6 2.42e-02   2.8 5.31e-02 5.56e-03h  1
     937r 8.1382801e+02 9.93e+01 1.72e+03  -4.6 2.42e-02   2.3 5.35e-02 6.46e-02h  1
     938r 8.1388619e+02 9.52e+01 1.01e+03  -4.6 2.32e-02   1.8 9.56e-02 4.15e-02h  1
     939r 8.1389650e+02 9.43e+01 1.91e+03  -4.6 2.26e-02   2.3 3.54e-01 9.41e-03h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     940r 8.1408565e+02 7.55e+01 1.87e+03  -4.6 2.25e-02   1.8 9.98e-01 2.04e-01h  1
     941r 8.1429516e+02 6.10e+01 1.46e+03  -4.6 1.96e-02   1.3 2.35e-01 1.95e-01h  1
     942r 8.1445689e+02 5.15e+01 7.04e+02  -4.6 1.93e-02   1.7 7.04e-02 1.58e-01h  1
     943r 8.1447479e+02 4.98e+01 5.53e+02  -4.6 1.53e-02   1.3 1.71e-01 3.43e-02h  1
     944r 8.1500779e+02 3.97e+01 8.92e+02  -4.6 5.30e-02   0.8 5.70e-03 2.92e-01f  1
     945r 8.1505268e+02 3.87e+01 4.24e+02  -4.6 3.63e-02   1.2 6.89e-02 2.53e-02h  1
     946r 8.1505569e+02 3.80e+01 6.07e+02  -4.6 4.27e-03   2.5 1.17e-01 1.78e-02h  1
     947r 8.1507722e+02 3.61e+01 5.60e+02  -4.6 9.45e-03   2.1 2.45e-02 4.89e-02h  1
     948r 8.1506430e+02 3.61e+01 5.65e+02  -4.6 1.22e+05    -  7.31e-04 3.96e-05f  1
     949r 8.1507587e+02 3.61e+01 5.68e+02  -4.6 2.38e+04    -  1.15e-03 6.70e-04f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     950r 8.1515087e+02 3.60e+01 5.53e+02  -4.6 7.77e+03    -  6.35e-04 1.81e-03f  1
     951r 8.1522701e+02 3.60e+01 5.52e+02  -4.6 9.27e+03    -  1.82e-03 1.81e-03f  1
     952r 8.1547062e+02 3.58e+01 1.72e+03  -4.6 6.09e+03    -  8.04e-01 5.26e-03f  1
     953r 8.1691065e+02 6.64e+01 1.64e+03  -4.6 5.50e+03    -  1.73e-04 3.78e-02f  1
     954r 8.1814549e+02 8.94e+01 1.55e+03  -4.6 2.45e+03    -  1.95e-02 4.68e-02f  1
     955r 8.1818237e+02 8.90e+01 1.61e+03  -4.6 2.39e+03    -  1.81e-01 3.71e-03h  1
     956r 8.1824706e+02 8.84e+01 1.60e+03  -4.6 6.56e+02    -  1.95e-03 6.88e-03f  1
     957r 8.1831862e+02 8.77e+01 1.59e+03  -4.6 5.94e+02    -  1.46e-02 7.87e-03f  1
     958r 8.1832743e+02 8.76e+01 1.62e+03  -4.6 3.50e+03    -  1.20e-01 1.13e-03h  1
     959r 8.1845015e+02 8.60e+01 1.59e+03  -4.6 7.66e+02    -  1.00e+00 1.84e-02f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     960r 8.1878698e+02 7.90e+01 1.41e+03  -4.6 7.89e+01    -  1.00e+00 8.28e-02h  1
     961r 8.1912423e+02 6.73e+01 9.10e+02  -4.6 6.25e+01    -  1.00e+00 1.53e-01h  1
     962r 8.2057526e+02 2.20e+01 2.30e+02  -4.6 4.13e+01    -  1.00e+00 8.19e-01h  1
     963r 8.2152766e+02 2.86e+01 1.97e+02  -4.6 9.52e+01    -  1.00e+00 8.90e-01h  1
     964r 8.2217928e+02 6.91e+01 3.04e+01  -4.6 1.29e+02    -  1.00e+00 1.00e+00h  1
     965r 8.2246677e+02 2.08e+00 1.58e-01  -4.6 1.99e+01    -  1.00e+00 1.00e+00h  1
     966r 8.2261297e+02 2.25e-01 1.96e-02  -4.6 8.65e-02    -  1.00e+00 1.00e+00h  1
     967r 8.2268881e+02 2.25e-01 4.47e-02  -4.6 1.27e-03    -  1.00e+00 1.00e+00h  1
     968r 8.2270757e+02 2.25e-01 4.60e-02  -4.6 1.41e-03    -  1.00e+00 1.00e+00h  1
     969r 8.2269864e+02 2.25e-01 1.57e-02  -4.6 8.86e-04    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     970r 8.2269998e+02 2.25e-01 4.58e-04  -4.6 1.09e-04    -  1.00e+00 1.00e+00h  1
     971r 8.2272857e+02 2.25e-01 1.37e-01  -6.9 3.83e-02    -  1.00e+00 9.98e-01f  1
     972r 8.2328669e+02 5.62e+02 1.65e+02  -6.9 1.71e+03    -  1.00e+00 2.10e-01f  1
     973r 8.2559199e+02 1.11e+04 3.41e+02  -6.9 1.45e+03    -  1.00e+00 1.00e+00f  1
     974r 8.2602079e+02 1.51e+03 2.62e+02  -6.9 2.44e+02    -  1.00e+00 1.00e+00h  1
     975r 8.2600853e+02 2.84e+00 6.78e+00  -6.9 2.18e+01    -  1.00e+00 1.00e+00h  1
     976r 8.2600947e+02 2.24e-01 1.43e-03  -6.9 5.98e-01    -  1.00e+00 1.00e+00h  1
     977r 8.2600947e+02 2.24e-01 7.28e-10  -6.9 2.52e-04    -  1.00e+00 1.00e+00h  1
     978r 8.2600912e+02 2.24e-01 1.10e-01  -9.0 2.16e-04    -  1.00e+00 9.93e-01f  1
     979r 8.2712388e+02 1.07e+04 6.36e+02  -9.0 7.65e+02    -  1.00e+00 1.00e+00f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     980r 8.2731695e+02 3.31e+03 3.97e+02  -9.0 3.75e+02    -  1.00e+00 1.00e+00h  1
     981r 8.2731592e+02 2.60e+03 3.12e+02  -9.0 7.79e+01    -  2.64e-01 2.13e-01h  1
     982r 8.2731582e+02 2.60e+03 1.21e+03  -9.0 6.11e+01    -  6.24e-03 1.13e-03h  1
     983r 8.2731582e+02 2.60e+03 1.21e+03  -9.0 2.71e+01   1.6 0.00e+00 4.29e-07R  5
     984r 8.2729971e+02 1.95e+03 1.20e+03  -9.0 5.91e+01   1.1 7.08e-03 7.02e-03f  1
     985r 8.2725152e+02 6.74e+00 1.14e+03  -9.0 2.48e+01   1.5 4.84e-09 5.01e-02f  1
     986r 8.2725122e+02 1.51e+00 1.13e+03  -9.0 6.31e-01   1.1 1.14e-01 6.83e-03f  1
     987r 8.2724777e+02 1.30e+00 1.00e+03  -9.0 2.24e-01   1.5 2.97e-03 1.20e-01f  1
     988r 8.2724489e+02 1.14e+00 8.85e+02  -9.0 1.12e-01   1.0 3.10e-02 1.22e-01f  1
     989r 8.2724493e+02 1.14e+00 8.70e+02  -9.0 1.20e+05    -  1.27e-01 1.27e-05f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
     990r 8.2724497e+02 1.14e+00 1.15e+03  -9.0 1.20e+05    -  9.64e-01 1.09e-04f  1
     991r 8.2724513e+02 9.62e-01 1.20e+03  -9.0 1.50e-02   0.5 1.46e-02 1.50e-01f  1
     992r 8.2724513e+02 9.62e-01 1.18e+03  -9.0 1.27e-02   1.9 9.19e-02 2.18e-04f  1
     993r 8.2724631e+02 9.18e-01 1.13e+03  -9.0 1.27e-02   1.4 4.41e-02 4.53e-02f  1
     994r 8.2724923e+02 7.57e-01 9.97e+02  -9.0 1.21e-02   1.8 1.17e-02 1.73e-01f  1
     995r 8.2724978e+02 5.69e-01 6.88e+02  -9.0 1.00e-02   1.3 3.29e-01 2.46e-01f  1
     996r 8.2724705e+02 3.33e-01 5.54e+02  -9.0 7.54e-03   0.9 1.04e-01 4.12e-01f  1
     997r 8.2724610e+02 2.37e-01 9.85e+02  -9.0 4.42e-03   1.3 8.87e-02 2.84e-01f  1
     998r 8.2724493e+02 2.23e-01 9.40e+02  -9.0 3.16e-03   1.7 2.77e-01 4.40e-01f  1
     999r 8.2724446e+02 2.23e-01 8.24e+02  -9.0 2.24e-03   1.2 3.50e-01 2.17e-01f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1000r 8.2724339e+02 2.23e-01 2.78e+02  -9.0 1.92e-03   1.7 7.09e-01 6.75e-01f  1
    1001r 8.2724288e+02 2.23e-01 7.05e+01  -9.0 1.04e-03   2.1 8.72e-01 7.42e-01f  1
    1002r 8.2724327e+02 2.23e-01 1.89e+02  -9.0 8.67e-04   1.6 4.01e-03 6.82e-02f  1
    1003r 8.2725657e+02 2.23e-01 5.04e+02  -9.0 6.61e-04   2.0 5.15e-02 5.58e-01f  1
    1004r 8.2725938e+02 2.23e-01 3.49e+02  -9.0 5.34e-04   1.6 3.38e-01 2.46e-01f  1
    1005r 8.2726416e+02 2.23e-01 2.55e+02  -9.0 9.97e-04   1.1 3.69e-01 5.79e-01f  1
    1006r 8.2726685e+02 2.23e-01 1.33e+02  -9.0 2.96e-03   0.6 7.74e-01 8.44e-01f  1
    1007r 8.2726687e+02 2.23e-01 6.88e+02  -9.0 1.41e-04   1.9 9.99e-01 4.30e-02f  1
    1008r 8.2726770e+02 2.23e-01 4.73e-01  -9.0 4.23e-04   1.5 1.00e+00 1.00e+00f  1
    1009r 8.2726767e+02 2.23e-01 2.16e+00  -9.0 1.27e-03   1.0 1.00e+00 9.91e-01f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1010r 8.2726768e+02 2.23e-01 1.48e+02  -9.0 3.81e-03   0.5 1.00e+00 2.26e-01f  1
    1011r 8.2726768e+02 2.23e-01 3.32e+02  -9.0 1.09e-02   0.0 4.19e-01 1.77e-03f  1
    1012r 8.2726768e+02 2.23e-01 3.26e+02  -9.0 3.43e-02  -0.5 1.96e-02 3.45e-02f  1
    1013r 8.2726769e+02 2.23e-01 2.98e+02  -9.0 9.53e-02  -0.9 3.31e-02 6.67e-02f  1
    1014r 8.2726774e+02 2.23e-01 2.40e+02  -9.0 2.86e-01  -1.4 3.13e-02 1.94e-01f  1
    1015r 8.2726796e+02 2.23e-01 2.16e+02  -9.0 8.56e-01  -1.9 8.26e-02 2.48e-01f  1
    1016r 8.2726872e+02 2.23e-01 2.68e+02  -9.0 2.55e+00  -2.4 1.21e-01 3.00e-01f  1
    1017r 8.2727027e+02 3.33e-01 2.39e+02  -9.0 7.52e+00  -2.8 1.61e-01 2.06e-01f  1
    1018r 8.2727083e+02 3.62e-01 2.33e+02  -9.0 2.15e+01  -3.3 2.57e-02 2.67e-02f  1
    1019r 8.2727088e+02 3.62e-01 2.33e+02  -9.0 5.71e+01  -3.8 2.43e-07 8.34e-04f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1020r 8.2727089e+02 3.62e-01 2.33e+02  -9.0 1.64e+02  -4.3 3.74e-04 4.52e-05f  1
    1021r 8.2727089e+02 3.62e-01 2.33e+02  -9.0 4.97e+02  -4.7 1.12e-04 4.57e-06f  1
    1022r 8.2727093e+02 3.62e-01 2.33e+02  -9.0 1.53e+03  -5.2 5.98e-05 1.47e-04f  1
    1023r 8.2727096e+02 3.62e-01 2.33e+02  -9.0 5.03e+03  -5.7 7.84e-05 1.28e-04f  1
    1024r 8.2727096e+02 3.62e-01 2.33e+02  -9.0 2.41e+04  -6.2 1.16e-04 1.01e-05f  1
    1025r 8.2727097e+02 3.62e-01 2.32e+02  -9.0 5.77e+03  -5.8 7.97e-04 5.40e-06f  1
    1026r 8.2727097e+02 3.62e-01 2.32e+02  -9.0 3.16e+04  -6.2 9.12e-05 5.24e-05f  1
    1027r 8.2727097e+02 3.62e-01 2.22e+02  -9.0 5.82e+03  -5.8 2.44e-02 1.64e-05f  1
    1028r 8.2727099e+02 3.63e-01 2.22e+02  -9.0 1.68e+04  -6.3 1.57e-08 2.28e-04f  1
    1029r 8.2727120e+02 3.72e-01 2.22e+02  -9.0 5.37e+03  -5.9 7.89e-04 9.09e-04f  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
    1030r 8.2727122e+02 3.72e-01 2.21e+02  -9.0 2.30e+04  -6.3 8.10e-04 1.76e-04f  1
    1031r 8.2727122e+02 3.72e-01 7.47e+02  -9.0 1.48e+05    -  1.00e+00 4.38e-06f  1
    1032r 8.2727128e+02 3.73e-01 7.78e+02  -9.0 1.40e+05    -  1.00e+00 1.44e-04f  1
    1033r 8.2727129e+02 3.73e-01 1.00e+03  -9.0 9.05e+04    -  1.00e+00 4.33e-05f  1
    1034r 8.2727484e+02 3.24e-01 1.03e+03  -9.0 4.66e+01    -  1.00e+00 1.55e-01f  1
    1035r 8.2729603e+02 4.84e-01 2.18e+02  -9.0 4.27e+01    -  1.00e+00 1.00e+00f  1
    1036r 8.2732101e+02 8.73e+01 1.01e+01  -9.0 3.50e+01    -  1.00e+00 1.00e+00h  1
    1037r 8.2732035e+02 2.23e-01 1.63e-02  -9.0 9.68e-01    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 1037
    
                                       (scaled)                 (unscaled)
    Objective...............:   8.2732034904976831e+02    8.2732034904976831e+02
    Dual infeasibility......:   2.3895987687923457e+02    2.3895987687923457e+02
    Constraint violation....:   2.2343556576761836e-01    2.2343556576761836e-01
    Complementarity.........:   9.0909579504895325e-10    9.0909579504895325e-10
    Overall NLP error.......:   2.3895987687923457e+02    2.3895987687923457e+02
    
    
    Number of objective function evaluations             = 1614
    Number of objective gradient evaluations             = 107
    Number of equality constraint evaluations            = 1628
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 1053
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 1038
    Total CPU secs in IPOPT (w/o function evaluations)   =      3.126
    Total CPU secs in NLP function evaluations           =      2.509
    
    EXIT: Converged to a point of local infeasibility. Problem may be infeasible.
    WARNING: Loading a SolverResults object with a warning status into
        model=MASTER;
            message from solver=Ipopt 3.13.2\x3a Converged to a locally infeasible
            point. Problem may be infeasible.
    

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
          <td>0.478960</td>
          <td>-0.409798</td>
        </tr>
        <tr>
          <th>1</th>
          <td>0.447785</td>
          <td>-0.385781</td>
        </tr>
        <tr>
          <th>2</th>
          <td>0.460895</td>
          <td>-0.396118</td>
        </tr>
        <tr>
          <th>3</th>
          <td>0.475832</td>
          <td>-0.407662</td>
        </tr>
        <tr>
          <th>4</th>
          <td>0.445738</td>
          <td>-0.384256</td>
        </tr>
        <tr>
          <th>5</th>
          <td>0.485444</td>
          <td>-0.414666</td>
        </tr>
        <tr>
          <th>6</th>
          <td>0.483045</td>
          <td>-0.412851</td>
        </tr>
        <tr>
          <th>7</th>
          <td>0.491341</td>
          <td>-0.419172</td>
        </tr>
        <tr>
          <th>8</th>
          <td>0.475045</td>
          <td>-0.407220</td>
        </tr>
        <tr>
          <th>9</th>
          <td>-0.453060</td>
          <td>5.000000</td>
        </tr>
      </tbody>
    </table>
    </div>

