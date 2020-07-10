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

    2020-07-10 21:50:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:50:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:50:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:50:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:50:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:50:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:50:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:50:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:50:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:50:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:50:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:50:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:50:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:50:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:50:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:50:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:50:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.014
    Total CPU secs in NLP function evaluations           =      0.013
    
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

    2020-07-10 21:51:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  6.9788893e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  1.1683758e-01 1.40e+03 5.27e-01  -1.0 1.37e+04    -  9.89e-01 1.00e+00f  1
       2  6.3189660e-03 9.68e+01 3.50e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  6.4162349e-03 3.43e+01 2.73e-04  -1.7 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  6.4041252e-03 1.36e+00 1.17e-04  -2.5 1.43e-02    -  1.00e+00 1.00e+00h  1
       5  6.2778266e-03 3.89e+01 4.13e-04  -3.8 7.15e-02    -  1.00e+00 1.00e+00h  1
       6  5.8288504e-03 4.90e+02 3.83e-03  -3.8 2.59e-01  -4.0 1.00e+00 1.00e+00h  1
       7  4.8912569e-03 6.44e+02 1.24e-02  -3.8 3.33e-01  -3.6 1.00e+00 1.00e+00h  1
       8  4.8184370e-03 1.63e+02 1.72e-03  -3.8 1.60e-01    -  1.00e+00 1.00e+00h  1
       9  4.7708867e-03 7.62e+00 4.54e-05  -3.8 3.21e-02    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.7687732e-03 7.84e-05 8.47e-08  -3.8 1.77e-04    -  1.00e+00 1.00e+00h  1
      11  4.7672447e-03 2.34e-01 7.49e-06  -5.7 6.33e-03    -  1.00e+00 1.00e+00h  1
      12  4.7672152e-03 8.99e-05 1.24e-09  -5.7 1.18e-04    -  1.00e+00 1.00e+00h  1
      13  4.7672150e-03 3.32e-05 1.10e-09  -8.6 7.56e-05    -  1.00e+00 1.00e+00h  1
      14  4.7672150e-03 5.82e-11 2.13e-14  -8.6 1.70e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 14
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.7672149859795858e-03    4.7672149859795858e-03
    Dual infeasibility......:   2.1323829129745596e-14    2.1323829129745596e-14
    Constraint violation....:   2.1526712877120750e-13    5.8207660913467407e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 15
    Number of objective gradient evaluations             = 15
    Number of equality constraint evaluations            = 15
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 15
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 14
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.010
    Total CPU secs in NLP function evaluations           =      0.012
    
    EXIT: Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:12 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:13 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:14 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:15 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:16 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:17 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  6.2014158e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  9.3288263e-02 1.41e+03 4.96e-01  -1.0 1.37e+04    -  9.88e-01 1.00e+00f  1
       2  6.5487948e-03 2.32e+01 1.77e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  6.5731669e-03 2.28e-01 5.32e-05  -2.5 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  5.3778282e-03 2.02e+03 1.77e-02  -3.8 5.30e-01    -  9.52e-01 1.00e+00h  1
       5  5.3599335e-03 3.97e+02 2.73e-03  -3.8 2.29e-01    -  1.00e+00 1.00e+00h  1
       6  5.0962186e-03 1.71e+02 6.05e-03  -3.8 1.79e-01    -  1.00e+00 1.00e+00h  1
       7  4.8783595e-03 4.97e+01 1.17e-03  -3.8 9.08e-02    -  1.00e+00 1.00e+00h  1
       8  4.8665368e-03 3.47e+00 5.76e-05  -3.8 2.33e-02    -  1.00e+00 1.00e+00h  1
       9  4.8662732e-03 1.14e-02 1.61e-07  -3.8 1.31e-03    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.8648792e-03 1.82e-01 7.22e-06  -5.7 5.60e-03    -  1.00e+00 1.00e+00h  1
      11  4.8648564e-03 5.16e-05 7.32e-10  -5.7 8.96e-05    -  1.00e+00 1.00e+00h  1
      12  4.8648562e-03 2.62e-05 1.07e-09  -8.6 6.72e-05    -  1.00e+00 1.00e+00h  1
      13  4.8648562e-03 2.91e-11 2.19e-14  -8.6 1.30e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.8648561856265166e-03    4.8648561856265166e-03
    Dual infeasibility......:   2.1876396249383213e-14    2.1876396249383213e-14
    Constraint violation....:   5.6843418860808015e-14    2.9103830456733704e-11
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
    Total CPU secs in NLP function evaluations           =      0.013
    
    EXIT: Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:18 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:19 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:20 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:21 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:22 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:23 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:24 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  6.5073316e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  1.0258592e-01 1.41e+03 4.87e-01  -1.0 1.37e+04    -  9.88e-01 1.00e+00f  1
       2  7.1304950e-03 1.84e+01 1.73e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  7.1407299e-03 8.09e-01 1.01e-04  -2.5 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  5.2274801e-03 6.60e+03 5.26e-02  -3.8 1.92e+00    -  7.57e-01 5.00e-01h  2
       5  6.1565705e-03 3.73e+03 3.01e-02  -3.8 5.27e-01    -  9.78e-01 5.00e-01h  2
       6  6.5005428e-03 1.03e+03 5.82e-03  -3.8 3.80e-01    -  1.00e+00 1.00e+00h  1
       7  5.6571879e-03 8.23e+02 1.18e-02  -3.8 4.89e-01  -4.0 1.00e+00 5.00e-01h  2
       8  5.1606478e-03 1.67e+02 2.26e-03  -3.8 1.58e-01    -  1.00e+00 1.00e+00h  1
       9  5.1070388e-03 8.65e+00 8.10e-05  -3.8 3.48e-02    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  5.1053091e-03 1.11e-02 4.99e-08  -3.8 1.19e-03    -  1.00e+00 1.00e+00h  1
      11  5.1040265e-03 1.44e-01 6.45e-06  -5.7 5.00e-03    -  1.00e+00 1.00e+00h  1
      12  5.1040068e-03 3.10e-05 4.59e-10  -5.7 6.97e-05    -  1.00e+00 1.00e+00h  1
      13  5.1040066e-03 2.09e-05 9.56e-10  -8.6 6.03e-05    -  1.00e+00 1.00e+00h  1
      14  5.1040066e-03 2.91e-11 2.03e-14  -8.6 1.02e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 14
    
                                       (scaled)                 (unscaled)
    Objective...............:   5.1040065528444832e-03    5.1040065528444832e-03
    Dual infeasibility......:   2.0276268216203606e-14    2.0276268216203606e-14
    Constraint violation....:   5.6843418860808015e-14    2.9103830456733704e-11
    Complementarity.........:   2.5059035596800626e-09    2.5059035596800626e-09
    Overall NLP error.......:   2.5059035596800626e-09    2.5059035596800626e-09
    
    
    Number of objective function evaluations             = 20
    Number of objective gradient evaluations             = 15
    Number of equality constraint evaluations            = 20
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 15
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 14
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.018
    Total CPU secs in NLP function evaluations           =      0.022
    
    EXIT: Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:25 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:26 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:27 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:28 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:29 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:30 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  4.6133399e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  5.8143920e-02 1.42e+03 4.78e-01  -1.0 1.37e+04    -  9.88e-01 1.00e+00f  1
       2  5.1197209e-03 8.31e+00 1.75e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  5.1600986e-03 4.58e-02 2.19e-05  -2.5 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  4.5748547e-03 5.60e+02 4.98e-03  -3.8 2.81e-01    -  9.98e-01 1.00e+00h  1
       5  4.5291982e-03 1.67e+02 8.38e-03  -3.8 1.82e-01    -  1.00e+00 1.00e+00h  1
       6  4.1344120e-03 6.92e+01 1.96e-03  -3.8 1.08e-01    -  1.00e+00 1.00e+00h  1
       7  4.1112454e-03 9.63e+00 1.71e-04  -3.8 3.90e-02    -  1.00e+00 1.00e+00h  1
       8  4.1109492e-03 1.32e-01 1.72e-06  -3.8 4.48e-03    -  1.00e+00 1.00e+00h  1
       9  4.1092901e-03 2.98e-01 1.05e-05  -5.7 7.10e-03    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.1092623e-03 1.54e-04 1.81e-09  -5.7 1.54e-04    -  1.00e+00 1.00e+00h  1
      11  4.1092620e-03 4.25e-05 1.54e-09  -8.6 8.49e-05    -  1.00e+00 1.00e+00h  1
      12  4.1092620e-03 1.31e-10 2.20e-14  -8.6 2.21e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 12
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.1092620100202084e-03    4.1092620100202084e-03
    Dual infeasibility......:   2.2014777328921516e-14    2.2014777328921516e-14
    Constraint violation....:   4.8435103973521685e-13    1.3096723705530167e-10
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 13
    Number of objective gradient evaluations             = 13
    Number of equality constraint evaluations            = 13
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 13
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 12
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.011
    Total CPU secs in NLP function evaluations           =      0.010
    
    EXIT: Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:31 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:32 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:33 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:34 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:35 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:36 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:37 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  5.1397157e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  7.1050094e-02 1.42e+03 4.60e-01  -1.0 1.37e+04    -  9.87e-01 1.00e+00f  1
       2  5.8858099e-03 7.08e+00 1.88e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  5.8297497e-03 6.60e+00 3.60e-04  -3.8 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  5.8099600e-03 3.96e-02 2.59e-05  -3.8 2.59e-03  -2.0 1.00e+00 1.00e+00h  1
       5  5.7794454e-03 4.42e-01 2.74e-05  -5.7 8.02e-03  -2.5 1.00e+00 1.00e+00h  1
       6  5.6728785e-03 5.12e+00 6.01e-05  -5.7 2.74e-02  -3.0 1.00e+00 1.00e+00h  1
       7  5.1659565e-03 1.11e+02 1.21e-03  -5.7 1.28e-01  -3.4 1.00e+00 1.00e+00h  1
       8  4.3222399e-03 5.30e+02 9.17e-03  -5.7 4.44e-01  -3.9 1.00e+00 1.00e+00H  1
       9  4.7997763e-03 3.56e+02 8.38e-03  -5.7 3.07e-01    -  1.00e+00 1.00e+00H  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.4990235e-03 1.03e+02 3.61e-03  -5.7 1.21e-01    -  1.00e+00 1.00e+00h  1
      11  4.4607724e-03 8.86e+00 2.27e-04  -5.7 3.69e-02    -  1.00e+00 1.00e+00h  1
      12  4.4595334e-03 6.94e-02 1.05e-06  -5.7 3.22e-03    -  1.00e+00 1.00e+00h  1
      13  4.4595247e-03 2.42e-06 2.22e-11  -5.7 1.86e-05    -  1.00e+00 1.00e+00h  1
      14  4.4595245e-03 2.92e-05 1.30e-09  -8.6 7.07e-05    -  1.00e+00 1.00e+00h  1
      15  4.4595245e-03 8.73e-11 2.22e-14  -8.6 1.46e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 15
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.4595244815720187e-03    4.4595244815720187e-03
    Dual infeasibility......:   2.2173832739832745e-14    2.2173832739832745e-14
    Constraint violation....:   3.2290069315681125e-13    8.7311491370201111e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 19
    Number of objective gradient evaluations             = 16
    Number of equality constraint evaluations            = 19
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 16
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 15
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.012
    Total CPU secs in NLP function evaluations           =      0.011
    
    EXIT: Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:38 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:39 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:40 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:41 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:42 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:43 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:44 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  5.9421757e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  9.4159818e-02 1.42e+03 4.76e-01  -1.0 1.37e+04    -  9.88e-01 1.00e+00f  1
       2  6.6178921e-03 1.81e+01 1.74e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  6.6442300e-03 4.65e-01 8.13e-05  -2.5 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  4.9524114e-03 5.61e+03 4.37e-02  -3.8 8.86e-01    -  8.94e-01 1.00e+00h  1
       5  6.4480637e-03 1.50e+03 1.27e-02  -3.8 4.78e-01    -  1.00e+00 1.00e+00h  1
       6  5.0043833e-03 1.41e+02 2.91e-03  -3.8 3.93e-01    -  1.00e+00 1.00e+00h  1
       7  4.8492917e-03 5.32e+01 1.24e-03  -3.8 2.12e-01    -  1.00e+00 1.00e+00h  1
       8  4.8417837e-03 4.95e+00 8.31e-05  -3.8 2.80e-02    -  1.00e+00 1.00e+00h  1
       9  4.8414552e-03 2.44e-02 3.42e-07  -3.8 1.92e-03    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.8401080e-03 1.66e-01 7.32e-06  -5.7 5.34e-03    -  1.00e+00 1.00e+00h  1
      11  4.8400898e-03 4.18e-05 6.02e-10  -5.7 8.06e-05    -  1.00e+00 1.00e+00h  1
      12  4.8400896e-03 2.40e-05 1.08e-09  -8.6 6.44e-05    -  1.00e+00 1.00e+00h  1
      13  4.8400895e-03 5.82e-11 2.21e-14  -8.6 1.18e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 13
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.8400895485921375e-03    4.8400895485921375e-03
    Dual infeasibility......:   2.2084761545179689e-14    2.2084761545179689e-14
    Constraint violation....:   2.1526712877120750e-13    5.8207660913467407e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 14
    Number of objective gradient evaluations             = 14
    Number of equality constraint evaluations            = 14
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 14
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 13
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.012
    Total CPU secs in NLP function evaluations           =      0.010
    
    EXIT: Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:45 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:46 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:47 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:48 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:49 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:50 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:51 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  5.5744277e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  7.5513345e-02 1.42e+03 4.72e-01  -1.0 1.37e+04    -  9.87e-01 1.00e+00f  1
       2  6.3386122e-03 8.05e+00 1.78e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  6.3592343e-03 2.01e-01 5.16e-05  -2.5 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  4.8541715e-03 8.20e+03 6.11e-02  -3.8 1.08e+00    -  8.65e-01 1.00e+00h  1
       5  5.2744200e-03 6.29e+03 4.71e-02  -3.8 6.13e-01    -  1.00e+00 2.50e-01h  3
       6  5.8864333e-03 3.58e+03 2.70e-02  -3.8 5.26e-01    -  1.00e+00 5.00e-01h  2
       7  6.0418625e-03 1.04e+03 5.85e-03  -3.8 3.82e-01    -  1.00e+00 1.00e+00h  1
       8  5.3650957e-03 1.51e+03 1.62e-02  -3.8 2.57e+00  -4.0 1.00e+00 1.25e-01h  4
       9  4.8344614e-03 3.29e+02 3.17e-03  -3.8 2.18e-01    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.7239369e-03 1.63e+01 8.38e-05  -3.8 4.54e-02    -  1.00e+00 1.00e+00h  1
      11  4.7191363e-03 1.33e-03 4.24e-07  -3.8 9.92e-04    -  1.00e+00 1.00e+00h  1
      12  4.7177593e-03 1.78e-01 7.77e-06  -5.7 5.53e-03    -  1.00e+00 1.00e+00h  1
      13  4.7177377e-03 4.91e-05 6.97e-10  -5.7 8.73e-05    -  1.00e+00 1.00e+00h  1
      14  4.7177374e-03 2.58e-05 1.15e-09  -8.6 6.66e-05    -  1.00e+00 1.00e+00h  1
      15  4.7177374e-03 5.82e-11 2.22e-14  -8.6 1.28e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 15
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.7177374445802282e-03    4.7177374445802282e-03
    Dual infeasibility......:   2.2204291028011037e-14    2.2204291028011037e-14
    Constraint violation....:   2.1526712877120750e-13    5.8207660913467407e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 23
    Number of objective gradient evaluations             = 16
    Number of equality constraint evaluations            = 23
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 16
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 15
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.014
    Total CPU secs in NLP function evaluations           =      0.017
    
    EXIT: Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:52 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:53 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:54 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:55 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:56 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:57 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  5.1448320e+01 2.79e+00 3.68e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  5.9095889e-02 1.12e+03 3.38e-01  -1.0 1.22e+04    -  1.00e+00 1.00e+00f  1
       2  6.1278643e-03 2.36e+01 1.63e-02  -1.7 3.73e+02    -  1.00e+00 1.00e+00h  1
       3  6.1000723e-03 1.52e-01 4.58e-05  -2.5 4.14e-01    -  1.00e+00 1.00e+00h  1
       4  5.1298434e-03 1.27e+03 1.03e-02  -3.8 4.19e-01    -  9.72e-01 1.00e+00h  1
       5  4.6735003e-03 1.15e+02 2.38e-03  -3.8 1.26e-01    -  1.00e+00 1.00e+00h  1
       6  4.5398033e-03 2.62e+01 7.74e-04  -3.8 6.91e-02    -  1.00e+00 1.00e+00h  1
       7  4.5371874e-03 2.17e+00 3.88e-05  -3.8 1.87e-02    -  1.00e+00 1.00e+00h  1
       8  4.5370913e-03 6.51e-03 9.00e-08  -3.8 9.99e-04    -  1.00e+00 1.00e+00h  1
       9  4.5356595e-03 2.00e-01 7.66e-06  -5.7 5.85e-03    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.5356347e-03 6.33e-05 7.86e-10  -5.7 9.90e-05    -  1.00e+00 1.00e+00h  1
      11  4.5356344e-03 2.86e-05 1.13e-09  -8.6 7.01e-05    -  1.00e+00 1.00e+00h  1
      12  4.5356344e-03 5.82e-11 2.21e-14  -8.6 1.44e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 12
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.5356344278953357e-03    4.5356344278953357e-03
    Dual infeasibility......:   2.2074644433110812e-14    2.2074644433110812e-14
    Constraint violation....:   2.1526712877120750e-13    5.8207660913467407e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 13
    Number of objective gradient evaluations             = 13
    Number of equality constraint evaluations            = 13
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 13
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 12
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.011
    Total CPU secs in NLP function evaluations           =      0.011
    
    EXIT: Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:58 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:51:59 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:00 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:01 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:02 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:03 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:04 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  5.3999422e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  8.5102697e-02 1.41e+03 4.95e-01  -1.0 1.37e+04    -  9.88e-01 1.00e+00f  1
       2  5.5627990e-03 1.96e+01 1.75e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  5.6105739e-03 6.70e-04 1.86e-05  -2.5 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  5.0450376e-03 5.16e+02 4.86e-03  -3.8 2.68e-01    -  1.00e+00 1.00e+00h  1
       5  4.7127357e-03 1.28e+02 4.62e-03  -3.8 3.16e-01  -4.0 1.00e+00 1.00e+00H  1
       6  4.7319940e-03 1.26e+02 4.36e-03  -3.8 4.47e-01  -4.5 1.00e+00 6.25e-02h  5
       7  4.9909027e-03 7.73e+01 2.40e-03  -3.8 1.97e-01  -4.1 1.00e+00 1.00e+00H  1
       8  4.4719913e-03 6.15e+02 9.65e-03  -3.8 1.47e+01    -  3.63e-01 1.99e-02h  5
       9  4.4119484e-03 1.30e+02 1.31e-03  -3.8 1.41e-01    -  1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  4.3851158e-03 5.62e+00 2.66e-05  -3.8 2.77e-02    -  1.00e+00 1.00e+00h  1
      11  4.3840994e-03 1.65e-03 2.80e-08  -3.8 5.01e-04    -  1.00e+00 1.00e+00h  1
      12  4.3825085e-03 2.67e-01 9.23e-06  -5.7 6.73e-03    -  1.00e+00 1.00e+00h  1
      13  4.3824820e-03 1.19e-04 1.59e-09  -5.7 1.36e-04    -  1.00e+00 1.00e+00h  1
      14  4.3824817e-03 3.77e-05 1.35e-09  -8.6 8.02e-05    -  1.00e+00 1.00e+00h  1
      15  4.3824817e-03 8.73e-11 2.09e-14  -8.6 1.94e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 15
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.3824817357682620e-03    4.3824817357682620e-03
    Dual infeasibility......:   2.0926587372469633e-14    2.0926587372469633e-14
    Constraint violation....:   3.2290069315681125e-13    8.7311491370201111e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 28
    Number of objective gradient evaluations             = 16
    Number of equality constraint evaluations            = 28
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 16
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 15
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.023
    Total CPU secs in NLP function evaluations           =      0.027
    
    EXIT: Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:05 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:06 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:07 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:08 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:09 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 1 optimal - Optimal Solution Found.
    2020-07-10 21:52:10 [INFO] idaes.init.fs.state_block: Initialization Step 2 optimal - Optimal Solution Found.
    2020-07-10 21:52:11 [INFO] idaes.init.fs.state_block: Initialization Step 3 optimal - Optimal Solution Found.
    2020-07-10 21:52:11 [INFO] idaes.init.fs.state_block: Initialization Step 4 optimal - Optimal Solution Found.
    2020-07-10 21:52:11 [INFO] idaes.init.fs.state_block: Initialization Step 5 optimal - Optimal Solution Found.
    2020-07-10 21:52:11 [INFO] idaes.init.fs.state_block: State Released.
    2020-07-10 21:52:11 [INFO] idaes.init.fs.state_block: Initialization Complete: optimal - Optimal Solution Found
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
       0  6.9891101e+01 3.15e+00 4.16e+01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
       1  1.0914129e-01 1.40e+03 5.19e-01  -1.0 1.37e+04    -  9.88e-01 1.00e+00f  1
       2  6.7330302e-03 6.40e+01 3.20e-02  -1.7 4.74e+02    -  1.00e+00 1.00e+00h  1
       3  6.9472606e-03 7.38e+01 6.19e-04  -1.7 6.62e-01    -  1.00e+00 1.00e+00h  1
       4  6.8677701e-03 1.27e+01 1.17e-04  -1.7 4.25e-02    -  1.00e+00 1.00e+00h  1
       5  6.8336554e-03 8.60e-01 1.04e-04  -2.5 1.11e-02    -  1.00e+00 1.00e+00h  1
       6  6.6410505e-03 7.44e+01 7.81e-04  -3.8 9.94e-02    -  1.00e+00 1.00e+00h  1
       7  6.5876448e-03 7.22e-04 3.49e-05  -3.8 3.47e-03  -2.0 1.00e+00 1.00e+00h  1
       8  6.5682888e-03 2.69e-01 2.11e-05  -5.7 6.18e-03  -2.5 1.00e+00 1.00e+00h  1
       9  6.5030444e-03 3.06e+00 3.93e-05  -5.7 2.09e-02  -3.0 1.00e+00 1.00e+00h  1
    iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
      10  6.1747901e-03 6.69e+01 7.61e-04  -5.7 9.82e-02  -3.4 1.00e+00 1.00e+00h  1
      11  5.8339283e-03 1.31e+01 2.85e-04  -5.7 4.66e-02  -3.0 1.00e+00 1.00e+00h  1
      12  5.1066858e-03 1.84e+02 3.40e-03  -5.7 1.70e-01  -3.5 1.00e+00 1.00e+00h  1
      13  4.9662816e-03 2.46e+00 3.35e-04  -5.7 2.57e-02    -  1.00e+00 1.00e+00h  1
      14  4.9630170e-03 4.45e-01 1.22e-05  -5.7 8.73e-03    -  1.00e+00 1.00e+00h  1
      15  4.9629477e-03 5.88e-04 1.17e-08  -5.7 3.07e-04    -  1.00e+00 1.00e+00h  1
      16  4.9629474e-03 2.71e-05 9.97e-10  -8.6 6.85e-05    -  1.00e+00 1.00e+00h  1
      17  4.9629474e-03 5.82e-11 2.12e-14  -8.6 1.36e-08    -  1.00e+00 1.00e+00h  1
    
    Number of Iterations....: 17
    
                                       (scaled)                 (unscaled)
    Objective...............:   4.9629473998646671e-03    4.9629473998646671e-03
    Dual infeasibility......:   2.1168218711018691e-14    2.1168218711018691e-14
    Constraint violation....:   2.1526712877120750e-13    5.8207660913467407e-11
    Complementarity.........:   2.5059035596800622e-09    2.5059035596800622e-09
    Overall NLP error.......:   2.5059035596800622e-09    2.5059035596800622e-09
    
    
    Number of objective function evaluations             = 18
    Number of objective gradient evaluations             = 18
    Number of equality constraint evaluations            = 18
    Number of inequality constraint evaluations          = 0
    Number of equality constraint Jacobian evaluations   = 18
    Number of inequality constraint Jacobian evaluations = 0
    Number of Lagrangian Hessian evaluations             = 17
    Total CPU secs in IPOPT (w/o function evaluations)   =      0.020
    Total CPU secs in NLP function evaluations           =      0.022
    
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
          <td>0.473655</td>
          <td>-0.406260</td>
        </tr>
        <tr>
          <th>1</th>
          <td>0.485703</td>
          <td>-0.415158</td>
        </tr>
        <tr>
          <th>2</th>
          <td>0.500225</td>
          <td>-0.426065</td>
        </tr>
        <tr>
          <th>3</th>
          <td>0.448204</td>
          <td>-0.386185</td>
        </tr>
        <tr>
          <th>4</th>
          <td>0.472713</td>
          <td>-0.404900</td>
        </tr>
        <tr>
          <th>5</th>
          <td>0.488732</td>
          <td>-0.417268</td>
        </tr>
        <tr>
          <th>6</th>
          <td>0.483488</td>
          <td>-0.413248</td>
        </tr>
        <tr>
          <th>7</th>
          <td>0.475450</td>
          <td>-0.407060</td>
        </tr>
        <tr>
          <th>8</th>
          <td>0.459248</td>
          <td>-0.394895</td>
        </tr>
        <tr>
          <th>9</th>
          <td>0.486293</td>
          <td>-0.415815</td>
        </tr>
      </tbody>
    </table>
    </div>

