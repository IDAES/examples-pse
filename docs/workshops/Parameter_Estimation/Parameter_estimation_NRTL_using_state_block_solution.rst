Parameter Estimation Using the NRTL State Block
-----------------------------------------------

In this module, we use Pyomo's ``parmest`` tool in conjunction with
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
   https://idaes-pse.readthedocs.io/en/latest/model\_libraries/core\_library/property\_models/activity\_coefficient.html
-  parmest -
   https://pyomo.readthedocs.io/en/stable/contributed\_packages/parmest/index.html

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
pass the method that returns an initialized model, data, variable\_name,
and the SSE expression to the Estimator method. ``tee=True`` will print
the solver output after solving the parameter estimation problem.

.. code:: ipython3

    # Initialize a parameter estimation object
    pest = parmest.Estimator(NRTL_model, data, variable_name, SSE, tee=True)
    
    # Run parameter estimation using all data
    obj_value, parameters = pest.theta_est()


.. parsed-literal::

    2020-07-21 06:13:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.065
    Total CPU secs in NLP function evaluations           =      0.036
    
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

    2020-07-21 06:13:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:13:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:13:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:13:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:13:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:13:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:13:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:13:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  6.0157501e+01 3.03e+00 4.00e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  7.7117934e-02 1.32e+03 4.19e-01  -1.0 1.32e+04    -  9.94e-01 1.00e+00f  1
       2  7.0324020e-03 1.40e+01 1.61e-02  -1.7 4.39e+02    -  1.00e+00 1.00e+00h  1
       3  7.0108469e-03 8.40e-01 1.08e-04  -2.5 5.70e-01    -  1.00e+00 1.00e+00h  1
       4  5.0928148e-03 5.94e+03 4.49e-02  -3.8 3.65e+00    -  6.03e-01 2.50e-01h  3
       5  6.7847742e-03 1.55e+03 1.18e-02  -3.8 4.89e-01    -  9.52e-01 1.00e+00h  1
       6  5.2065248e-03 1.84e+02 2.89e-03  -3.8 5.09e-01    -  1.00e+00 1.00e+00h  1
       7  5.0011213e-03 7.13e+01 1.87e-03  -3.8 3.10e-01    -  1.00e+00 1.00e+00h  1
       8  4.9817219e-03 8.72e+00 1.68e-04  -3.8 3.74e-02    -  1.00e+00 1.00e+00h  1
       9  4.9808190e-03 7.98e-02 1.27e-06  -3.8 3.49e-03    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.9795440e-03 1.40e-01 6.60e-06  -5.7 4.92e-03    -  1.00e+00 1.00e+00h  1
      11  4.9795252e-03 2.91e-05 4.52e-10  -5.7 6.74e-05    -  1.00e+00 1.00e+00h  1
      12  4.9795250e-03 2.05e-05 9.87e-10  -8.6 5.97e-05    -  1.00e+00 1.00e+00h  1
      13  4.9795250e-03 4.37e-11 2.14e-14  -8.6 9.96e-09    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.9795250235212304e-03    4.9795250235212304e-03
    Dual infeasibility......:   2.1432523782360335e-14    2.1432523782360335e-14
    Constraint violation....:   1.6145034657840562e-13    4.3655745685100555e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 17
    Number of objective gradient evaluations             = 14
    Number of equality constraint evaluations            = 17
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 14
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 13
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.058
    Total CPU secs in NLP function evaluations           =      0.081
    
    EXIT: Optimal Solution Found.
    2020-07-21 06:14:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  4.2215727e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  5.6450829e-02 1.43e+03 4.49e-01  -1.0 1.37e+04    -  9.87e-01 1.00e+00f  1
       2  4.8451957e-03 4.84e+00 1.97e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  4.8137249e-03 3.75e+00 2.79e-04  -3.8 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  3.9400727e-03 3.50e+03 2.81e-02  -3.8 1.41e+00  -4.0 1.00e+00 5.00e-01h  2
       5  4.3599128e-03 9.18e+02 6.93e-03  -3.8 3.69e-01    -  1.00e+00 1.00e+00h  1
       6  3.9264802e-03 8.33e+01 1.22e-03  -3.8 1.07e-01    -  1.00e+00 1.00e+00h  1
       7  3.8762124e-03 6.98e+00 2.10e-04  -3.8 3.62e-02    -  1.00e+00 1.00e+00h  1
       8  3.8768864e-03 2.39e-01 4.13e-06  -3.8 6.16e-03    -  1.00e+00 1.00e+00h  1
       9  3.8752400e-03 2.94e-01 1.14e-05  -5.7 7.04e-03    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  3.8752162e-03 1.49e-04 1.78e-09  -5.7 1.51e-04    -  1.00e+00 1.00e+00h  1
      11  3.8752159e-03 4.28e-05 1.70e-09  -8.6 8.50e-05    -  1.00e+00 1.00e+00h  1
      12  3.8752159e-03 2.91e-11 2.08e-14  -8.6 2.22e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 12
    
                                       (scaled)                 (unscaled)
    Objective...............:   3.8752159063630020e-03    3.8752159063630020e-03
    Dual infeasibility......:   2.0793674331493123e-14    2.0793674331493123e-14
    Constraint violation....:   1.0763356438560375e-13    2.9103830456733704e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 15
    Number of objective gradient evaluations             = 13
    Number of equality constraint evaluations            = 15
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 13
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 12
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.060
    Total CPU secs in NLP function evaluations           =      0.026
    
    EXIT: Optimal Solution Found.
    2020-07-21 06:14:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:14:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:14:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:14:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:14:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:14:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:14:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:14:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  6.6484432e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  1.0137527e-01 1.41e+03 5.02e-01  -1.0 1.37e+04    -  9.88e-01 1.00e+00f  1
       2  6.9526960e-03 3.80e+01 2.11e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  6.9562359e-03 5.41e-01 8.29e-05  -2.5 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  5.8477763e-03 1.58e+03 1.29e-02  -3.8 4.66e-01    -  9.64e-01 1.00e+00h  1
       5  5.3234238e-03 1.90e+02 3.09e-03  -3.8 1.60e-01    -  1.00e+00 1.00e+00h  1
       6  5.0796110e-03 8.19e+01 2.39e-03  -3.8 1.22e-01    -  1.00e+00 1.00e+00h  1
       7  5.0428186e-03 1.34e+01 2.63e-04  -3.8 4.67e-02    -  1.00e+00 1.00e+00h  1
       8  5.0414169e-03 2.26e-01 3.36e-06  -3.8 5.90e-03    -  1.00e+00 1.00e+00h  1
       9  5.0400549e-03 1.58e-01 6.47e-06  -5.7 5.24e-03    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  5.0400339e-03 3.90e-05 5.34e-10  -5.7 7.81e-05    -  1.00e+00 1.00e+00h  1
      11  5.0400337e-03 2.36e-05 9.77e-10  -8.6 6.40e-05    -  1.00e+00 1.00e+00h  1
      12  5.0400337e-03 8.73e-11 2.21e-14  -8.6 1.17e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 12
    
                                       (scaled)                 (unscaled)
    Objective...............:   5.0400337131028365e-03    5.0400337131028365e-03
    Dual infeasibility......:   2.2076995207619794e-14    2.2076995207619794e-14
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.022
    Total CPU secs in NLP function evaluations           =      0.067
    
    EXIT: Optimal Solution Found.
    2020-07-21 06:15:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  4.2203464e+01 3.03e+00 4.00e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  5.4599569e-02 1.32e+03 4.25e-01  -1.0 1.32e+04    -  9.94e-01 1.00e+00f  1
       2  4.7770482e-03 1.49e+01 1.57e-02  -1.7 4.39e+02    -  1.00e+00 1.00e+00h  1
       3  4.8085219e-03 7.85e-02 2.94e-05  -2.5 5.70e-01    -  1.00e+00 1.00e+00h  1
       4  4.4272952e-03 2.75e+02 2.41e-03  -3.8 1.96e-01    -  1.00e+00 1.00e+00h  1
       5  3.8957997e-03 4.79e+02 6.60e-03  -3.8 4.71e-01  -4.0 1.00e+00 5.00e-01h  2
       6  3.8939681e-03 8.29e+01 6.13e-04  -3.8 1.10e-01    -  1.00e+00 1.00e+00h  1
       7  3.8782037e-03 1.06e+00 1.25e-05  -3.8 1.17e-02    -  1.00e+00 1.00e+00h  1
       8  3.8779091e-03 2.31e-04 1.06e-08  -3.8 2.25e-04    -  1.00e+00 1.00e+00h  1
       9  3.8761662e-03 3.56e-01 1.15e-05  -5.7 7.73e-03    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  3.8761365e-03 2.24e-04 2.50e-09  -5.7 1.85e-04    -  1.00e+00 1.00e+00h  1
      11  3.8761362e-03 4.94e-05 1.66e-09  -8.6 9.13e-05    -  1.00e+00 1.00e+00h  1
      12  3.8761362e-03 8.73e-11 2.18e-14  -8.6 2.62e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 12
    
                                       (scaled)                 (unscaled)
    Objective...............:   3.8761362269565279e-03    3.8761362269565279e-03
    Dual infeasibility......:   2.1821869950969789e-14    2.1821869950969789e-14
    Constraint violation....:   3.2290069315681125e-13    8.7311491370201111e-11
    Complementarity.........:   2.5059035596800618e-09    2.5059035596800618e-09
    Overall NLP error.......:   2.5059035596800618e-09    2.5059035596800618e-09
    
    
    Number of objective function evaluations             = 15
    Number of objective gradient evaluations             = 13
    Number of equality constraint evaluations            = 15
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 13
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 12
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.042
    Total CPU secs in NLP function evaluations           =      0.044
    
    EXIT: Optimal Solution Found.
    2020-07-21 06:15:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  4.1292264e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  5.4094285e-02 1.42e+03 4.73e-01  -1.0 1.37e+04    -  9.87e-01 1.00e+00f  1
       2  4.6612705e-03 5.31e+00 1.80e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  4.6412656e-03 2.58e+00 2.03e-04  -3.8 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  3.8907048e-03 3.17e+03 2.75e-02  -3.8 6.72e-01  -4.0 1.00e+00 1.00e+00h  1
       5  4.2981493e-03 8.30e+02 6.67e-03  -3.8 3.49e-01    -  1.00e+00 1.00e+00h  1
       6  3.8761357e-03 6.24e+01 9.56e-04  -3.8 9.31e-02    -  1.00e+00 1.00e+00h  1
       7  3.8390291e-03 4.76e+00 1.19e-04  -3.8 2.97e-02    -  1.00e+00 1.00e+00h  1
       8  3.8395048e-03 1.10e-01 1.57e-06  -3.8 4.16e-03    -  1.00e+00 1.00e+00h  1
       9  3.8377394e-03 3.66e-01 1.21e-05  -5.7 7.84e-03    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  3.8377052e-03 2.41e-04 2.80e-09  -5.7 1.92e-04    -  1.00e+00 1.00e+00h  1
      11  3.8377049e-03 5.15e-05 1.77e-09  -8.6 9.31e-05    -  1.00e+00 1.00e+00h  1
      12  3.8377049e-03 2.91e-11 1.99e-14  -8.6 2.74e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 12
    
                                       (scaled)                 (unscaled)
    Objective...............:   3.8377048694318741e-03    3.8377048694318741e-03
    Dual infeasibility......:   1.9877921978422441e-14    1.9877921978422441e-14
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.027
    Total CPU secs in NLP function evaluations           =      0.044
    
    EXIT: Optimal Solution Found.
    2020-07-21 06:15:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:15:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:15:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:15:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:15:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:15:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:15:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:15:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  6.0009075e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  8.6716236e-02 1.42e+03 4.76e-01  -1.0 1.37e+04    -  9.87e-01 1.00e+00f  1
       2  6.7651989e-03 1.04e+01 1.75e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  6.7793313e-03 4.98e-01 8.10e-05  -2.5 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  5.0264679e-03 7.95e+03 6.19e-02  -3.8 2.12e+00    -  7.36e-01 5.00e-01h  2
       5  6.1110338e-03 4.50e+03 3.60e-02  -3.8 5.97e-01    -  9.53e-01 5.00e-01h  2
       6  6.6613366e-03 1.33e+03 8.60e-03  -3.8 4.39e-01    -  1.00e+00 1.00e+00h  1
       7  5.8206607e-03 1.66e+03 1.79e-02  -3.8 2.39e+00  -4.0 1.00e+00 1.25e-01h  4
       8  5.0637808e-03 3.68e+02 3.54e-03  -3.8 2.29e-01    -  1.00e+00 1.00e+00h  1
       9  4.9278766e-03 1.58e+01 8.64e-05  -3.8 4.53e-02    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.9224936e-03 1.03e-02 8.96e-07  -3.8 1.77e-03    -  1.00e+00 1.00e+00h  1
      11  4.9211838e-03 1.54e-01 7.01e-06  -5.7 5.17e-03    -  1.00e+00 1.00e+00h  1
      12  4.9211640e-03 3.60e-05 5.34e-10  -5.7 7.49e-05    -  1.00e+00 1.00e+00h  1
      13  4.9211638e-03 2.26e-05 1.05e-09  -8.6 6.25e-05    -  1.00e+00 1.00e+00h  1
      14  4.9211638e-03 4.37e-11 2.01e-14  -8.6 1.11e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 14
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.9211637665423800e-03    4.9211637665423800e-03
    Dual infeasibility......:   2.0060520099347981e-14    2.0060520099347981e-14
    Constraint violation....:   1.6145034657840562e-13    4.3655745685100555e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 22
    Number of objective gradient evaluations             = 15
    Number of equality constraint evaluations            = 22
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 15
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 14
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.036
    Total CPU secs in NLP function evaluations           =      0.048
    
    EXIT: Optimal Solution Found.
    2020-07-21 06:16:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  6.4918181e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  9.1781987e-02 1.41e+03 4.98e-01  -1.0 1.37e+04    -  9.88e-01 1.00e+00f  1
       2  6.8918502e-03 1.74e+01 1.76e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  6.9047759e-03 3.63e-01 6.32e-05  -2.5 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  5.1427705e-03 7.32e+03 6.01e-02  -3.8 1.01e+00    -  8.74e-01 1.00e+00h  1
       5  6.2805441e-03 4.14e+03 3.49e-02  -3.8 5.66e-01    -  1.00e+00 5.00e-01h  2
       6  6.7578003e-03 1.21e+03 7.98e-03  -3.8 4.15e-01    -  1.00e+00 1.00e+00h  1
       7  5.1056612e-03 2.30e+00 2.14e-03  -3.8 1.01e-01  -2.0 1.00e+00 1.00e+00h  1
       8  5.0392702e-03 1.66e+00 1.52e-04  -3.8 2.58e-02    -  1.00e+00 1.00e+00h  1
       9  5.0378724e-03 1.32e-02 7.05e-07  -3.8 1.56e-03    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  5.0365432e-03 1.60e-01 6.60e-06  -5.7 5.26e-03    -  1.00e+00 1.00e+00h  1
      11  5.0365201e-03 3.92e-05 5.72e-10  -5.7 7.82e-05    -  1.00e+00 1.00e+00h  1
      12  5.0365199e-03 2.33e-05 9.81e-10  -8.6 6.35e-05    -  1.00e+00 1.00e+00h  1
      13  5.0365199e-03 8.73e-11 2.21e-14  -8.6 1.15e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   5.0365198596233228e-03    5.0365198596233228e-03
    Dual infeasibility......:   2.2131566278838233e-14    2.2131566278838233e-14
    Constraint violation....:   3.2290069315681125e-13    8.7311491370201111e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 15
    Number of objective gradient evaluations             = 14
    Number of equality constraint evaluations            = 15
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 14
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 13
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.042
    Total CPU secs in NLP function evaluations           =      0.043
    
    EXIT: Optimal Solution Found.
    2020-07-21 06:16:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  5.0803979e+01 2.91e+00 3.84e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  5.8533405e-02 1.22e+03 3.77e-01  -1.0 1.27e+04    -  1.00e+00 1.00e+00f  1
       2  6.0770291e-03 4.76e+01 2.35e-02  -1.7 4.05e+02    -  1.00e+00 1.00e+00h  1
       3  6.0278805e-03 2.51e-01 6.91e-05  -2.5 4.87e-01    -  1.00e+00 1.00e+00h  1
       4  5.4571252e-03 4.57e+02 3.46e-03  -3.8 2.51e-01    -  1.00e+00 1.00e+00h  1
       5  4.8655858e-03 4.53e-01 2.49e-04  -3.8 1.33e-02  -2.0 1.00e+00 1.00e+00h  1
       6  4.4018745e-03 5.11e+02 6.39e-04  -3.8 2.89e-01    -  1.00e+00 1.00e+00H  1
       7  4.4939511e-03 1.78e+02 8.60e-04  -3.8 1.53e-01    -  1.00e+00 1.00e+00h  1
       8  4.4738163e-03 2.17e+00 2.01e-05  -3.8 1.76e-02    -  1.00e+00 1.00e+00h  1
       9  4.4732490e-03 1.38e-04 4.22e-09  -3.8 2.14e-04    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.4717895e-03 2.11e-01 8.20e-06  -5.7 6.00e-03    -  1.00e+00 1.00e+00h  1
      11  4.4717664e-03 7.09e-05 8.83e-10  -5.7 1.05e-04    -  1.00e+00 1.00e+00h  1
      12  4.4717662e-03 3.01e-05 1.20e-09  -8.6 7.18e-05    -  1.00e+00 1.00e+00h  1
      13  4.4717662e-03 1.46e-11 2.20e-14  -8.6 1.51e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.4717661526214392e-03    4.4717661526214392e-03
    Dual infeasibility......:   2.1950068576133907e-14    2.1950068576133907e-14
    Constraint violation....:   5.6843418860808015e-14    1.4551915228366852e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 15
    Number of objective gradient evaluations             = 14
    Number of equality constraint evaluations            = 15
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 14
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 13
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.115
    Total CPU secs in NLP function evaluations           =      0.028
    
    EXIT: Optimal Solution Found.
    2020-07-21 06:16:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:16:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:16:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:16:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:16:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:16:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:16:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:16:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  4.6211384e+01 2.91e+00 3.84e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  5.0725517e-02 1.22e+03 3.63e-01  -1.0 1.27e+04    -  1.00e+00 1.00e+00f  1
       2  5.5184453e-03 1.54e+01 1.50e-02  -1.7 4.05e+02    -  1.00e+00 1.00e+00h  1
       3  5.5081301e-03 2.07e-02 2.01e-05  -2.5 4.87e-01    -  1.00e+00 1.00e+00h  1
       4  4.7046708e-03 9.17e+02 7.32e-03  -3.8 3.58e-01    -  9.84e-01 1.00e+00h  1
       5  4.2170533e-03 1.92e+01 1.04e-03  -3.8 5.94e-02    -  1.00e+00 1.00e+00h  1
       6  4.2020413e-03 7.76e-01 1.67e-05  -3.8 1.20e-02    -  1.00e+00 1.00e+00h  1
       7  4.2020127e-03 1.09e-03 9.23e-09  -3.8 4.08e-04    -  1.00e+00 1.00e+00h  1
       8  4.2004920e-03 2.39e-01 9.36e-06  -5.7 6.37e-03    -  1.00e+00 1.00e+00h  1
       9  4.2004671e-03 9.35e-05 1.16e-09  -5.7 1.20e-04    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.2004668e-03 3.40e-05 1.37e-09  -8.6 7.61e-05    -  1.00e+00 1.00e+00h  1
      11  4.2004668e-03 1.02e-10 2.19e-14  -8.6 1.73e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 11
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.2004668177067269e-03    4.2004668177067269e-03
    Dual infeasibility......:   2.1936684218932273e-14    2.1936684218932273e-14
    Constraint violation....:   3.7671747534961315e-13    1.0186340659856798e-10
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 12
    Number of objective gradient evaluations             = 12
    Number of equality constraint evaluations            = 12
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 12
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 11
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.013
    Total CPU secs in NLP function evaluations           =      0.046
    
    EXIT: Optimal Solution Found.
    2020-07-21 06:17:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-21 06:17:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-21 06:17:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-21 06:17:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-21 06:17:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-21 06:17:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-21 06:17:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-21 06:17:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  5.8069564e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  7.6962465e-02 1.42e+03 4.70e-01  -1.0 1.37e+04    -  9.87e-01 1.00e+00f  1
       2  6.6404860e-03 9.01e+00 1.80e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  6.6497947e-03 4.22e-01 7.60e-05  -2.5 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  4.9276656e-03 6.63e+03 5.16e-02  -3.8 1.94e+00    -  7.56e-01 5.00e-01h  2
       5  5.6985028e-03 3.75e+03 2.95e-02  -3.8 5.31e-01    -  9.76e-01 5.00e-01h  2
       6  6.0126872e-03 1.03e+03 5.78e-03  -3.8 3.83e-01    -  1.00e+00 1.00e+00h  1
       7  5.3357816e-03 7.07e+02 9.67e-03  -3.8 3.92e-01    -  1.00e+00 5.00e-01h  2
       8  4.8827587e-03 1.31e+02 1.95e-03  -3.8 1.39e-01    -  1.00e+00 1.00e+00h  1
       9  4.8473810e-03 6.08e+00 7.39e-05  -3.8 2.93e-02    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.8464633e-03 8.57e-03 6.76e-08  -3.8 1.06e-03    -  1.00e+00 1.00e+00h  1
      11  4.8451349e-03 1.60e-01 7.33e-06  -5.7 5.26e-03    -  1.00e+00 1.00e+00h  1
      12  4.8451155e-03 3.89e-05 5.83e-10  -5.7 7.78e-05    -  1.00e+00 1.00e+00h  1
      13  4.8451153e-03 2.32e-05 1.08e-09  -8.6 6.33e-05    -  1.00e+00 1.00e+00h  1
      14  4.8451153e-03 1.46e-10 2.16e-14  -8.6 1.14e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 14
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.8451153122791735e-03    4.8451153122791735e-03
    Dual infeasibility......:   2.1556129727055630e-14    2.1556129727055630e-14
    Constraint violation....:   5.3816782192801870e-13    1.4551915228366849e-10
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 19
    Number of objective gradient evaluations             = 15
    Number of equality constraint evaluations            = 19
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 15
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 14
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.109
    Total CPU secs in NLP function evaluations           =      0.021
    
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
          <td>0.498685</td>
          <td>-0.424704</td>
        </tr>
        <tr>
          <th>1</th>
          <td>0.443437</td>
          <td>-0.382180</td>
        </tr>
        <tr>
          <th>2</th>
          <td>0.493806</td>
          <td>-0.421352</td>
        </tr>
        <tr>
          <th>3</th>
          <td>0.436215</td>
          <td>-0.376787</td>
        </tr>
        <tr>
          <th>4</th>
          <td>0.433841</td>
          <td>-0.374938</td>
        </tr>
        <tr>
          <th>5</th>
          <td>0.493270</td>
          <td>-0.420706</td>
        </tr>
        <tr>
          <th>6</th>
          <td>0.494588</td>
          <td>-0.421912</td>
        </tr>
        <tr>
          <th>7</th>
          <td>0.471588</td>
          <td>-0.404111</td>
        </tr>
        <tr>
          <th>8</th>
          <td>0.460986</td>
          <td>-0.395813</td>
        </tr>
        <tr>
          <th>9</th>
          <td>0.490593</td>
          <td>-0.418624</td>
        </tr>
      </tbody>
    </table>
    </div>

